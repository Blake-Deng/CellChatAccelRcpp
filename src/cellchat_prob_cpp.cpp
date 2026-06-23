#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
List cellchat_prob_boot_cpp(
    NumericMatrix dataLavg,
    NumericMatrix dataRavg,
    NumericMatrix dataLavgBoot,
    NumericMatrix dataRavgBoot,
    NumericMatrix agonist,
    NumericMatrix antagonist,
    NumericMatrix agonistBoot,
    NumericMatrix antagonistBoot,
    LogicalVector hasAgonist,
    LogicalVector hasAntagonist,
    double Kh,
    double n_power,
    int nboot) {

  const int nLR = dataLavg.nrow();
  const int K = dataLavg.ncol();
  const int KK = K * K;
  NumericVector prob(Dimension(K, K, nLR));
  NumericVector pval(Dimension(K, K, nLR));

  const double Kh_n = std::pow(Kh, n_power);

  for (int lr = 0; lr < nLR; ++lr) {
    std::vector<double> pnull(KK);
    bool all_zero = true;

    for (int s = 0; s < K; ++s) {
      const double L = dataLavg(lr, s);
      for (int r = 0; r < K; ++r) {
        const double R = dataRavg(lr, r);
        const double dataLR = L * R;
        double p1 = 0.0;
        if (dataLR > 0.0) {
          const double dataLR_n = std::pow(dataLR, n_power);
          p1 = dataLR_n / (Kh_n + dataLR_n);
        }
        double val = p1;
        if (hasAgonist[lr]) val *= agonist(lr, s) * agonist(lr, r);
        if (hasAntagonist[lr]) val *= antagonist(lr, s) * antagonist(lr, r);

        const int idx = s + K * r;
        pnull[idx] = val;
        prob[s + K * r + KK * lr] = val;
        if (val != 0.0) all_zero = false;
      }
    }

    if (all_zero) {
      for (int idx = 0; idx < KK; ++idx) {
        pval[idx + KK * lr] = 1.0;
      }
      continue;
    }

    std::vector<int> reject(KK, 0);
    for (int b = 0; b < nboot; ++b) {
      const int row = lr + nLR * b;
      for (int s = 0; s < K; ++s) {
        const double L = dataLavgBoot(row, s);
        for (int r = 0; r < K; ++r) {
          const double R = dataRavgBoot(row, r);
          const double dataLR = L * R;
          double p1 = 0.0;
          if (dataLR > 0.0) {
            const double dataLR_n = std::pow(dataLR, n_power);
            p1 = dataLR_n / (Kh_n + dataLR_n);
          }
          double val = p1;
          if (hasAgonist[lr]) val *= agonistBoot(row, s) * agonistBoot(row, r);
          if (hasAntagonist[lr]) val *= antagonistBoot(row, s) * antagonistBoot(row, r);
          const int idx = s + K * r;
          if (val - pnull[idx] > 0.0) reject[idx] += 1;
        }
      }
    }

    for (int idx = 0; idx < KK; ++idx) {
      double pv = static_cast<double>(reject[idx]) / static_cast<double>(nboot);
      if (pnull[idx] == 0.0) pv = 1.0;
      pval[idx + KK * lr] = pv;
    }
  }

  return List::create(
    Named("prob") = prob,
    Named("pval") = pval
  );
}

static inline double quantile_type7_sorted(const std::vector<double>& x, double p) {
  const int n = x.size();
  if (n == 0) return NA_REAL;
  if (n == 1) return x[0];
  const double h = 1.0 + (n - 1.0) * p;
  const int hf = static_cast<int>(std::floor(h));
  const double frac = h - hf;
  const int i0 = hf - 1;
  if (i0 < 0) return x[0];
  if (i0 >= n - 1) return x[n - 1];
  return x[i0] + frac * (x[i0 + 1] - x[i0]);
}

static inline double kth_value(std::vector<double> vals, int k) {
  std::nth_element(vals.begin(), vals.begin() + k, vals.end());
  return vals[k];
}

static inline double quantile_type7_unsorted(const std::vector<double>& x, double p) {
  const int n = x.size();
  if (n == 0) return NA_REAL;
  if (n == 1) return x[0];
  const double h = 1.0 + (n - 1.0) * p;
  const int hf = static_cast<int>(std::floor(h));
  const double frac = h - hf;
  const int i0 = hf - 1;
  if (i0 < 0) return kth_value(x, 0);
  if (i0 >= n - 1) return kth_value(x, n - 1);
  const double lo = kth_value(x, i0);
  const double hi = kth_value(x, i0 + 1);
  return lo + frac * (hi - lo);
}

static inline double tri_mean_sorted(const std::vector<double>& x) {
  if (x.empty()) return NA_REAL;
  const double q1 = quantile_type7_sorted(x, 0.25);
  const double q2 = quantile_type7_sorted(x, 0.50);
  const double q3 = quantile_type7_sorted(x, 0.75);
  return (q1 + 2.0 * q2 + q3) / 4.0;
}

static inline double tri_mean_unsorted(const std::vector<double>& x) {
  if (x.empty()) return NA_REAL;
  const double q1 = quantile_type7_unsorted(x, 0.25);
  const double q2 = quantile_type7_unsorted(x, 0.50);
  const double q3 = quantile_type7_unsorted(x, 0.75);
  return (q1 + 2.0 * q2 + q3) / 4.0;
}

// [[Rcpp::export]]
NumericMatrix group_tri_mean_cpp(NumericMatrix data, IntegerVector group, int K) {
  const int G = data.nrow();
  const int C = data.ncol();
  NumericMatrix out(G, K);

  std::vector< std::vector<int> > cells(K);
  for (int c = 0; c < C; ++c) {
    int g = group[c] - 1;
    if (g >= 0 && g < K) cells[g].push_back(c);
  }

  std::vector<double> vals;
  for (int gene = 0; gene < G; ++gene) {
    for (int k = 0; k < K; ++k) {
      vals.clear();
      vals.reserve(cells[k].size());
      for (int idx : cells[k]) {
        double v = data(gene, idx);
        if (!NumericVector::is_na(v)) vals.push_back(v);
      }
      out(gene, k) = tri_mean_unsorted(vals);
    }
  }
  return out;
}

// [[Rcpp::export]]
List group_tri_mean_boot_cpp(NumericMatrix data, IntegerMatrix groupBoot, int K) {
  const int G = data.nrow();
  const int C = data.ncol();
  const int B = groupBoot.ncol();
  NumericVector out(Dimension(G, K, B));

  std::vector< std::vector<int> > cells(K);
  std::vector<double> vals;
  for (int b = 0; b < B; ++b) {
    for (int k = 0; k < K; ++k) cells[k].clear();
    for (int c = 0; c < C; ++c) {
      int g = groupBoot(c, b) - 1;
      if (g >= 0 && g < K) cells[g].push_back(c);
    }
    const int offset = G * K * b;
    for (int gene = 0; gene < G; ++gene) {
      for (int k = 0; k < K; ++k) {
        vals.clear();
        vals.reserve(cells[k].size());
        for (int idx : cells[k]) {
          double v = data(gene, idx);
          if (!NumericVector::is_na(v)) vals.push_back(v);
        }
        out[gene + G * k + offset] = tri_mean_unsorted(vals);
      }
    }
  }
  return List::create(Named("avg_boot") = out);
}

static inline double avg_expr_matrix(const NumericMatrix& avg, const IntegerMatrix& idx, int row, int k) {
  double sum_log = 0.0;
  int n = 0;
  for (int j = 0; j < idx.ncol(); ++j) {
    const int g = idx(row, j) - 1;
    if (g < 0) continue;
    const double x = avg(g, k);
    if (NumericVector::is_na(x)) continue;
    sum_log += std::log(x);
    ++n;
  }
  if (n == 0) return NA_REAL;
  if (n == 1) return std::exp(sum_log);
  return std::exp(sum_log / static_cast<double>(n));
}

static inline double avg_expr_boot(const NumericVector& avgBoot, int G, int K,
                                   const IntegerMatrix& idx, int row, int k, int b) {
  double sum_log = 0.0;
  int n = 0;
  const int boot_offset = G * K * b;
  for (int j = 0; j < idx.ncol(); ++j) {
    const int g = idx(row, j) - 1;
    if (g < 0) continue;
    const double x = avgBoot[g + G * k + boot_offset];
    if (NumericVector::is_na(x)) continue;
    sum_log += std::log(x);
    ++n;
  }
  if (n == 0) return NA_REAL;
  if (n == 1) return std::exp(sum_log);
  return std::exp(sum_log / static_cast<double>(n));
}

static inline double coreceptor_factor_matrix(const NumericMatrix& avg, const IntegerMatrix& idx, int row, int k) {
  double out = 1.0;
  for (int j = 0; j < idx.ncol(); ++j) {
    const int g = idx(row, j) - 1;
    if (g < 0) continue;
    const double x = avg(g, k);
    if (NumericVector::is_na(x)) continue;
    out *= 1.0 + x;
  }
  return out;
}

static inline double coreceptor_factor_boot(const NumericVector& avgBoot, int G, int K,
                                            const IntegerMatrix& idx, int row, int k, int b) {
  double out = 1.0;
  const int boot_offset = G * K * b;
  for (int j = 0; j < idx.ncol(); ++j) {
    const int g = idx(row, j) - 1;
    if (g < 0) continue;
    const double x = avgBoot[g + G * k + boot_offset];
    if (NumericVector::is_na(x)) continue;
    out *= 1.0 + x;
  }
  return out;
}

static inline double agonist_factor_value(double x, double Kh_n, double n_power) {
  const double xn = std::pow(x, n_power);
  return 1.0 + xn / (Kh_n + xn);
}

static inline double antagonist_factor_value(double x, double Kh_n, double n_power) {
  const double xn = std::pow(x, n_power);
  return Kh_n / (Kh_n + xn);
}

static inline double cofactor_factor_matrix(const NumericMatrix& avg, const IntegerMatrix& idx, int row, int k,
                                            bool agonist_mode, double Kh_n, double n_power) {
  double out = 1.0;
  for (int j = 0; j < idx.ncol(); ++j) {
    const int g = idx(row, j) - 1;
    if (g < 0) continue;
    const double x = avg(g, k);
    if (NumericVector::is_na(x)) continue;
    out *= agonist_mode ? agonist_factor_value(x, Kh_n, n_power) : antagonist_factor_value(x, Kh_n, n_power);
  }
  return out;
}

static inline double cofactor_factor_boot(const NumericVector& avgBoot, int G, int K,
                                          const IntegerMatrix& idx, int row, int k, int b,
                                          bool agonist_mode, double Kh_n, double n_power) {
  double out = 1.0;
  const int boot_offset = G * K * b;
  for (int j = 0; j < idx.ncol(); ++j) {
    const int g = idx(row, j) - 1;
    if (g < 0) continue;
    const double x = avgBoot[g + G * k + boot_offset];
    if (NumericVector::is_na(x)) continue;
    out *= agonist_mode ? agonist_factor_value(x, Kh_n, n_power) : antagonist_factor_value(x, Kh_n, n_power);
  }
  return out;
}

// [[Rcpp::export]]
List cellchat_prob_from_avg_cpp(
    NumericMatrix avg,
    NumericVector avgBoot,
    IntegerMatrix ligandIdx,
    IntegerMatrix receptorIdx,
    IntegerMatrix coAIdx,
    IntegerMatrix coIIdx,
    IntegerMatrix agonistIdx,
    IntegerMatrix antagonistIdx,
    LogicalVector hasAgonist,
    LogicalVector hasAntagonist,
    double Kh,
    double n_power) {

  const IntegerVector dims = avgBoot.attr("dim");
  const int G = dims[0];
  const int K = dims[1];
  const int nboot = dims[2];
  const int nLR = ligandIdx.nrow();
  const int KK = K * K;
  const double Kh_n = std::pow(Kh, n_power);

  NumericVector prob(Dimension(K, K, nLR));
  NumericVector pval(Dimension(K, K, nLR));

  std::vector<double> L(K), R(K), ago(K), ant(K);
  std::vector<double> Lb(K), Rb(K), agoB(K), antB(K);

  for (int lr = 0; lr < nLR; ++lr) {
    for (int k = 0; k < K; ++k) {
      L[k] = avg_expr_matrix(avg, ligandIdx, lr, k);
      R[k] = avg_expr_matrix(avg, receptorIdx, lr, k);
      const double coA = coreceptor_factor_matrix(avg, coAIdx, lr, k);
      const double coI = coreceptor_factor_matrix(avg, coIIdx, lr, k);
      R[k] = R[k] * coA / coI;
      ago[k] = hasAgonist[lr] ? cofactor_factor_matrix(avg, agonistIdx, lr, k, true, Kh_n, n_power) : 1.0;
      ant[k] = hasAntagonist[lr] ? cofactor_factor_matrix(avg, antagonistIdx, lr, k, false, Kh_n, n_power) : 1.0;
    }

    std::vector<double> pnull(KK);
    bool all_zero = true;
    for (int s = 0; s < K; ++s) {
      for (int r = 0; r < K; ++r) {
        const double dataLR = L[s] * R[r];
        double p1 = 0.0;
        if (dataLR > 0.0) {
          const double dataLR_n = std::pow(dataLR, n_power);
          p1 = dataLR_n / (Kh_n + dataLR_n);
        }
        double val = p1 * ago[s] * ago[r] * ant[s] * ant[r];
        const int idx = s + K * r;
        pnull[idx] = val;
        prob[idx + KK * lr] = val;
        if (val != 0.0) all_zero = false;
      }
    }

    if (all_zero) {
      for (int idx = 0; idx < KK; ++idx) pval[idx + KK * lr] = 1.0;
      continue;
    }

    std::vector<int> reject(KK, 0);
    for (int b = 0; b < nboot; ++b) {
      for (int k = 0; k < K; ++k) {
        Lb[k] = avg_expr_boot(avgBoot, G, K, ligandIdx, lr, k, b);
        Rb[k] = avg_expr_boot(avgBoot, G, K, receptorIdx, lr, k, b);
        const double coA = coreceptor_factor_boot(avgBoot, G, K, coAIdx, lr, k, b);
        const double coI = coreceptor_factor_boot(avgBoot, G, K, coIIdx, lr, k, b);
        Rb[k] = Rb[k] * coA / coI;
        agoB[k] = hasAgonist[lr] ? cofactor_factor_boot(avgBoot, G, K, agonistIdx, lr, k, b, true, Kh_n, n_power) : 1.0;
        antB[k] = hasAntagonist[lr] ? cofactor_factor_boot(avgBoot, G, K, antagonistIdx, lr, k, b, false, Kh_n, n_power) : 1.0;
      }
      for (int s = 0; s < K; ++s) {
        for (int r = 0; r < K; ++r) {
          const double dataLR = Lb[s] * Rb[r];
          double p1 = 0.0;
          if (dataLR > 0.0) {
            const double dataLR_n = std::pow(dataLR, n_power);
            p1 = dataLR_n / (Kh_n + dataLR_n);
          }
          const double val = p1 * agoB[s] * agoB[r] * antB[s] * antB[r];
          const int idx = s + K * r;
          if (val - pnull[idx] > 0.0) reject[idx] += 1;
        }
      }
    }

    for (int idx = 0; idx < KK; ++idx) {
      double pv = static_cast<double>(reject[idx]) / static_cast<double>(nboot);
      if (pnull[idx] == 0.0) pv = 1.0;
      pval[idx + KK * lr] = pv;
    }
  }

  return List::create(Named("prob") = prob, Named("pval") = pval);
}

// [[Rcpp::export]]
List pathway_sum_cpp(NumericVector prob, NumericVector pval, IntegerVector pathway, int n_pathways, double thresh) {
  const IntegerVector dims = prob.attr("dim");
  const int K = dims[0];
  const int nLR = dims[2];
  const int KK = K * K;
  NumericVector prob_pathway(Dimension(K, K, n_pathways));
  NumericVector lr_sum(nLR);
  NumericVector pathway_sum(n_pathways);

  for (int lr = 0; lr < nLR; ++lr) {
    const int p = pathway[lr] - 1;
    if (p < 0 || p >= n_pathways) continue;
    double lr_total = 0.0;
    for (int idx = 0; idx < KK; ++idx) {
      const int src_tgt = idx;
      const int pos = src_tgt + KK * lr;
      double val = prob[pos];
      if (pval[pos] > thresh) val = 0.0;
      prob_pathway[src_tgt + KK * p] += val;
      lr_total += val;
      pathway_sum[p] += val;
    }
    lr_sum[lr] = lr_total;
  }

  return List::create(
    Named("prob_pathway") = prob_pathway,
    Named("lr_sum") = lr_sum,
    Named("pathway_sum") = pathway_sum
  );
}

// [[Rcpp::export]]
List aggregate_net_cpp(NumericVector prob, NumericVector pval, double thresh) {
  const IntegerVector dims = prob.attr("dim");
  const int K = dims[0];
  const int nLR = dims[2];
  const int KK = K * K;
  NumericMatrix count(K, K);
  NumericMatrix weight(K, K);

  for (int lr = 0; lr < nLR; ++lr) {
    for (int idx = 0; idx < KK; ++idx) {
      const int pos = idx + KK * lr;
      const double val = prob[pos];
      double pv = pval[pos];
      if (val == 0.0) pv = 1.0;
      if (pv < thresh && val > 0.0) {
        const int s = idx % K;
        const int r = idx / K;
        count(s, r) += 1.0;
        weight(s, r) += val;
      }
    }
  }

  return List::create(Named("count") = count, Named("weight") = weight);
}
