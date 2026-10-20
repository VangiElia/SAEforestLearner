// Fast, memory-light prediction for ranger regression forests.
//
// ranger::predict() sorts every covariate column before predicting (save.memory
// is hard-coded to FALSE) and converts the whole forest again on each call. With
// ~150 covariates and millions of population units this dominates run time.
// Here the forest is flattened once in R (flatten_forest()), packed into 16-byte
// nodes and traversed branch-free, several rows in lockstep, on std::threads.
// Trees are processed in groups small enough to stay in the L3 cache while all
// rows stream past them. Results are bit-identical to ranger (same comparisons,
// same summation order over trees).

#include <Rcpp.h>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <thread>
#include <vector>

using namespace Rcpp;

namespace {

struct Node {
  double split;  // +Inf for terminal nodes (x > +Inf is never true)
  int var;       // covariate index (0-based)
  int left;      // left child; right child is left + 1 (ranger layout); self for terminal
};

template <class F>
void run_threads(int nt, F f) {
  if (nt <= 1) {
    f();
    return;
  }
  std::vector<std::thread> pool;
  pool.reserve(nt);
  for (int k = 0; k < nt; ++k) pool.emplace_back(f);
  for (auto& th : pool) th.join();
}

}  // namespace

// Maximal depth of every tree (children always have larger ids than their parent in ranger).
// [[Rcpp::export(rng = false)]]
IntegerVector cpp_tree_depth(const IntegerVector& left, const IntegerVector& right,
                             const IntegerVector& tree_start) {
  const int ntree = tree_start.size();
  const int nnode = left.size();
  IntegerVector out(ntree);
  std::vector<int> depth;
  for (int t = 0; t < ntree; ++t) {
    const int off = tree_start[t];
    const int end = (t + 1 < ntree) ? tree_start[t + 1] : nnode;
    depth.assign(end - off, 0);
    int maxd = 0;
    for (int i = 0; i < end - off; ++i) {
      const int l = left[off + i], r = right[off + i];
      if (l == 0 && r == 0) {
        maxd = std::max(maxd, depth[i]);
        continue;
      }
      if (l <= i || r <= i || l >= end - off || r >= end - off)
        stop("Unexpected ranger forest layout (child id not larger than parent id).");
      depth[l] = depth[i] + 1;
      depth[r] = depth[i] + 1;
    }
    out[t] = maxd;
  }
  return out;
}

// Mean prediction over trees for every row of x (numeric matrix, columns in the
// order of the forest's independent variables).
// [[Rcpp::export(rng = false)]]
NumericVector cpp_predict_forest(const NumericMatrix& x,
                                 const IntegerVector& var,
                                 const NumericVector& val,
                                 const IntegerVector& left,
                                 const IntegerVector& right,
                                 const IntegerVector& tree_start,
                                 const IntegerVector& tree_depth,
                                 int num_threads,
                                 int block_rows = 256,
                                 double group_mb = 3.0) {
  const std::size_t n = x.nrow();
  const std::size_t p = x.ncol();
  const std::size_t ntree = tree_start.size();
  const std::size_t nnode = var.size();
  NumericVector out(n);
  if (n == 0 || ntree == 0) return out;

  // pack nodes (tree-local child ids)
  std::vector<Node> nodes(nnode);
  for (std::size_t t = 0; t < ntree; ++t) {
    const std::size_t off = tree_start[t];
    const std::size_t end = (t + 1 < ntree) ? (std::size_t)tree_start[t + 1] : nnode;
    for (std::size_t i = off; i < end; ++i) {
      const int loc = (int)(i - off);
      if (left[i] == 0 && right[i] == 0) {
        nodes[i].split = std::numeric_limits<double>::infinity();
        nodes[i].var = 0;
        nodes[i].left = loc;
      } else {
        if (right[i] != left[i] + 1)
          stop("Unexpected ranger forest layout (right child != left child + 1).");
        if (var[i] < 0 || (std::size_t)var[i] >= p)
          stop("Split variable index out of range: covariate columns do not match the forest.");
        nodes[i].split = val[i];
        nodes[i].var = var[i];
        nodes[i].left = left[i];
      }
    }
  }

  // tree groups that fit in cache
  std::vector<std::size_t> gstart(1, 0);
  {
    const double max_nodes = std::max(1.0, group_mb * 1048576.0 / sizeof(Node));
    std::size_t acc_nodes = 0;
    for (std::size_t t = 0; t < ntree; ++t) {
      const std::size_t sz = ((t + 1 < ntree) ? (std::size_t)tree_start[t + 1] : nnode) - tree_start[t];
      if (acc_nodes > 0 && acc_nodes + sz > max_nodes) {
        gstart.push_back(t);
        acc_nodes = 0;
      }
      acc_nodes += sz;
    }
    gstart.push_back(ntree);
  }

  const int nt = std::max(1, num_threads);
  const std::size_t block = std::max(16, block_rows);
  const std::size_t nblocks = (n + block - 1) / block;
  const double* X = x.begin();

  // column-major -> row-major copy (each row's covariates contiguous)
  std::vector<double> rows(n * p);
  {
    std::atomic<std::size_t> next(0);
    run_threads(std::min<int>(nt, (int)nblocks), [&]() {
      for (;;) {
        const std::size_t b = next.fetch_add(1);
        if (b >= nblocks) break;
        const std::size_t r0 = b * block;
        const std::size_t r1 = std::min(n, r0 + block);
        for (std::size_t j = 0; j < p; ++j) {
          const double* col = X + j * n;
          for (std::size_t i = r0; i < r1; ++i) rows[i * p + j] = col[i];
        }
      }
    });
  }

  std::vector<double> acc(n, 0.0);
  const double* PV = val.begin();
  const int* T = tree_start.begin();
  const int* D = tree_depth.begin();
  const Node* NB = nodes.data();
  const double INF = std::numeric_limits<double>::infinity();
  const std::size_t lane = 16;  // rows traversed in lockstep

  for (std::size_t g = 0; g + 1 < gstart.size(); ++g) {
    const std::size_t t0 = gstart[g], t1 = gstart[g + 1];
    std::atomic<std::size_t> next(0);
    run_threads(std::min<int>(nt, (int)nblocks), [&]() {
      int cur[lane];
      for (;;) {
        const std::size_t b = next.fetch_add(1);
        if (b >= nblocks) break;
        const std::size_t r0 = b * block;
        const std::size_t nb = std::min(n, r0 + block) - r0;
        for (std::size_t t = t0; t < t1; ++t) {
          const Node* tr = NB + T[t];
          const double* pv = PV + T[t];
          const int depth = D[t];
          for (std::size_t g0 = 0; g0 < nb; g0 += lane) {
            const std::size_t ng = std::min(lane, nb - g0);
            const double* xr = &rows[(r0 + g0) * p];
            for (std::size_t k = 0; k < ng; ++k) cur[k] = 0;
            for (int d = 0; d < depth; ++d) {
              bool done = true;
              for (std::size_t k = 0; k < ng; ++k) {
                const Node& nd = tr[cur[k]];
                cur[k] = nd.left + (xr[k * p + nd.var] > nd.split);
                done &= (nd.split == INF);
              }
              if (done) break;
            }
            for (std::size_t k = 0; k < ng; ++k) acc[r0 + g0 + k] += pv[cur[k]];
          }
        }
      }
    });
  }

  double* O = out.begin();
  for (std::size_t i = 0; i < n; ++i) O[i] = acc[i] / (double)ntree;
  return out;
}

// Sum of x by integer group g (1-based), accumulated into a vector of length ng.
// [[Rcpp::export(rng = false)]]
NumericVector cpp_group_sum(const NumericVector& x, const IntegerVector& g, int ng) {
  NumericVector out(ng);
  const std::size_t n = x.size();
  if ((std::size_t)g.size() != n) stop("x and g must have the same length.");
  for (std::size_t i = 0; i < n; ++i) {
    const int k = g[i];
    if (k < 1 || k > ng) stop("Group index out of range.");
    out[k - 1] += x[i];
  }
  return out;
}

// One replicate of block-bootstrapped unit-level errors for the whole population,
// without materialising them: population unit i (domain dom[i]) draws an error
// uniformly from pool[offset[d] .. offset[d] + K[d] - 1] (domain-specific residuals
// for sampled domains, all residuals otherwise). Returns the error sums per cell
// and the errors of the units listed in sel (1-based, sorted ascending).
// Optional heteroscedastic errors: if scale has length n, the error of unit i is
// pool[...] * scale[i] (pool then holds standardised residuals).
// Optional Box-Cox back-transformation (bc_lambda >= 0 and base of length n): instead of the
// error sums, E holds the cell sums of Tinv(base[i] + u_dom[d] + e), the bootstrap values on the
// original scale, with Tinv(t) = max(lambda t + 1, 0)^(1/lambda) - shift (exp(t) - shift for
// lambda = 0). e_sel always holds the errors (transformed scale). bc_lambda < 0: no transformation.
static inline double bc_inv(double t, double lambda, double shift) {
  if (lambda == 0.0) return std::exp(t) - shift;
  const double b = lambda * t + 1.0;
  return (b > 0.0 ? std::pow(b, 1.0 / lambda) : 0.0) - shift;
}

// Box-Cox back-transformation of predictions summed per cell: unit i gets
// h_k(eta[i] + u_dom[d]), with d = dom[i], k = k_dom[d] and h_k the k-th column of the table
// tab (regular grid starting at t0 with step dt), interpolated linearly; values outside the
// grid are clamped. Returns the cell sums and the number of clamped units.
// [[Rcpp::export]]
List cpp_bt_sum(const NumericVector& eta, const IntegerVector& dom, const NumericVector& u_dom,
                const IntegerVector& k_dom, const NumericMatrix& tab, double t0, double dt,
                const IntegerVector& cell, int ncell) {
  const std::size_t n = eta.size();
  if ((std::size_t)dom.size() != n || (std::size_t)cell.size() != n) stop("eta, dom and cell must have the same length.");
  const int ng = tab.nrow();
  const double xmax = (double)(ng - 1) - 1e-9;
  NumericVector out(ncell);
  double clamped = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const int d = dom[i] - 1;
    const int k = k_dom[d] - 1;
    double x = (eta[i] + u_dom[d] - t0) / dt;
    if (x < 0.0) { x = 0.0; clamped += 1; } else if (x > xmax) { x = xmax; clamped += 1; }
    const int j = (int)x;
    const double w = x - j;
    out[cell[i] - 1] += tab(j, k) * (1.0 - w) + tab(j + 1, k) * w;
  }
  return List::create(Named("sum") = out, Named("clamped") = clamped);
}

// [[Rcpp::export]]
List cpp_boot_errors(const IntegerVector& dom, const IntegerVector& cell, int ncell,
                     const IntegerVector& offset, const IntegerVector& K,
                     const NumericVector& pool, const IntegerVector& sel,
                     const NumericVector& scale = NumericVector(0),
                     const NumericVector& base = NumericVector(0),
                     const NumericVector& u_dom = NumericVector(0),
                     double bc_lambda = -1.0, double bc_shift = 0.0) {
  const std::size_t n = dom.size();
  const std::size_t ns = sel.size();
  if ((std::size_t)cell.size() != n) stop("dom and cell must have the same length.");
  const bool scaled = scale.size() > 0;
  if (scaled && (std::size_t)scale.size() != n) stop("scale must be empty or have the length of dom.");
  const bool transform = bc_lambda >= 0.0;
  if (transform && (std::size_t)base.size() != n) stop("base must have the length of dom for the back-transformation.");
  if (transform && u_dom.size() < K.size()) stop("u_dom must have one value per domain.");
  NumericVector E(ncell);
  NumericVector esel(ns);
  std::size_t j = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const int d = dom[i] - 1;
    const int k = (int)R_unif_index((double)K[d]);
    const double e = scaled ? pool[offset[d] + k] * scale[i] : pool[offset[d] + k];
    E[cell[i] - 1] += transform ? bc_inv(base[i] + u_dom[d] + e, bc_lambda, bc_shift) : e;
    while (j < ns && (std::size_t)(sel[j] - 1) == i) {
      esel[j] = e;
      ++j;
    }
  }
  if (j != ns) stop("sel must be sorted and within 1..length(dom).");
  return List::create(Named("E") = E, Named("e_sel") = esel);
}
