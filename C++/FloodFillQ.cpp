#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
IntegerMatrix patches(IntegerMatrix im_mat,
                      IntegerMatrix patch_mat,
                      int x, int y) {
  int new_queue[2][2];

  new_queue[0][0] = x;
  new_queue[0][1] = y;

  int old_len = 1;
  int new_len = 1;

  IntegerVector dx[] = {-1, 0, 1, -1, 1, -1, 0, 1};
  IntegerVector dy[] = {-1, -1, -1, 0, 0, 1, 1, 1};

  for(int qi = 0; qi < new_len; qi++) {
    for(int vi = 0; vi < 8; vi++) {
      old_len = vi * qi;
    }
  }

  return patch_mat;
}
