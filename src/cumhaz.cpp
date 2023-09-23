#include <Rcpp.h>
using namespace Rcpp;

inline int binarySearch(NumericVector vec, double target) {
  int low = 0;
  int high = vec.size() - 1;
  int result = -1;

  while (low <= high) {
    int mid = (low + high) / 2;
    if (vec[mid] <= target) {
      result = mid;
      low = mid + 1;
    } else {
      high = mid - 1;
    }
  }
  return result;
}

// [[Rcpp::export]]
NumericVector cumhaz(NumericVector Time, NumericVector lambda, NumericVector death_point) {
  int n = Time.size();
  NumericVector HHazard(n);
  double max_death_point = max(death_point);
  double min_death_point = min(death_point);
  for (int i = 0; i < n; i++) {
    if(Time[i] > max_death_point){
      HHazard[i] = 996.0;
    }
    else if(Time[i] < min_death_point){
      HHazard[i] = 0.0;
    }
    else {
      double sum = 0;
      // for(int k = 0; k < death_point.size(); k++) {
      //   if(Time[i] >= death_point[k]) {
      //     temp += lambda[k];
      //   } else {
      //     break;
      //   }
      // }
      int index = binarySearch(death_point, Time[i]);
      // HHazard[i] = sum(death_point[Range(0, index)]);
      for (int k = 0; k <= index; k++) {
        sum += lambda[k];
      }
      HHazard[i] = sum;
    }
  }
  return HHazard;
}
