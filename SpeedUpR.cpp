#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
CharacterVector cppCheckFour(DataFrame x) {
  NumericVector col1 = x["col1"];
  NumericVector col2 = x["col2"];
  NumericVector col3 = x["col3"];
  NumericVector col4 = x["col4"];
  int n = col1.size();
  CharacterVector out(n);
  for (int i=0; i<n; i++){
  if ((col1[i] + col2[i] + col3[i] + col4[i]) > 4) {
    out[i] = "greater";
  } else {
    out[i] = "smaller_or_equal";
  }
}
return out;
}
