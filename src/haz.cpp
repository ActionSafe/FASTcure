#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector haz_(NumericVector Time, NumericVector event, NumericVector death_point, NumericVector wexpb) {

  int n = death_point.size();
  NumericVector lambda(n);

  int j = 0;

  for (int i = 0; i < n; i++) {
      while (j < Time.size() && Time[j] < death_point[i]) {
        j++;
      }
      double sum = 0.0;
      for(int k = j; k < wexpb.size(); k++) {
        sum += wexpb[k];
      }
      // double temp = sum(temp3[Range(j, temp3.size() - 1)]);
      lambda[i] = event[i] / sum;
  }
  return lambda;
}
