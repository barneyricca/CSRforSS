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

// The flood_fill implemented in package:spatialwarnings just passes
//  the matrices. It does not use globals. However, let's use the globals
//  here, as we have set the matrix size to 150 x 150.

// Without the global, the proliferation due to recursion was a killer.

IntegerMatrix patch_nums(150, 150);

void flood_fill(const IntegerMatrix seg1,
                int row_pos, int col_pos, int current_species, int patch_num)
{
  if(patch_nums(row_pos,col_pos) != 0) {  // Already been done
    return;
  }

  if(seg1(row_pos, col_pos) != current_species) {  // Not this species
    return;
  }

  patch_nums(row_pos, col_pos) = patch_num;

  // No wraps!
  int row_north = row_pos - 1;
  int row_south = row_pos + 1;
  int col_east = col_pos + 1;
  int col_west = col_pos - 1;
  flood_fill(seg1, row_south, col_pos, current_species, patch_num);   // then search can go south...
  flood_fill(seg1, row_north, col_pos, current_species, patch_num);   // or north...
  flood_fill(seg1, row_pos, col_east, current_species, patch_num);    // or east...
  flood_fill(seg1, row_pos, col_west, current_species, patch_num);    // or west...
  flood_fill(seg1, row_south, col_east, current_species, patch_num);  // or southeast...
  flood_fill(seg1, row_north, col_west, current_species, patch_num);  // or northwest...
  flood_fill(seg1, row_north, col_east, current_species, patch_num);  // or northeast...
  flood_fill(seg1, row_south, col_west, current_species, patch_num);  // or southwest...

  return;
}


// [[Rcpp::export]]
IntegerMatrix patchsizes(IntegerMatrix seg) {
  //
  // Very important: seg must have a white border!!!!
  //

  int world_size = seg.nrow();   // Get the number of rows

  // patch_nums contains the patch number for each cell
  for(int i = 0; i < world_size; i++) {
    for(int j = 0; j < world_size; j++) {
      patch_nums(i,j) = 0;
    }
  }

  // seg is an n x n matrix with 0 (blank), 1, and 2.
  // I need to return a matrix that is 3 columns x however many
  //  rows: The rows are the number of patches of species s with
  //  size sz. We know the population of each species coming in,
  //  so we should pass the maximum of that number (pop_size_max)
  //  as a way of defining the matrix. Then we call flood_fill (modified
  //  to wrap)
  int patch_num = 1;

  // Guaranteed a white border
  for(int row = 1; row < world_size-1; row++) {
    for(int column = 1; column < world_size-1; column++) {
      if(patch_nums(row, column) == 0) {
        flood_fill(seg, row, column, seg(row,column), patch_num);
        patch_num++;
      }
    }
  }

//  return patch_size_counts;
  return patch_nums;
}


// // [[Rcpp::export]]
// IntegerMatrix interior(IntegerMatrix seg) {
//   IntegerMatrix interior_counts(seg.nrow(), seg.ncol());
//   int row, column;
//   int row_north, row_south;
//   int col_east, col_west;

//   for(row = 0; row < seg.nrow(); row++) {
//     for(column = 0; column < seg.ncol(); column++) {
//       interior_counts(row, column) = -1;
//     }
//   }

// Count the number of interior points, and return a vector of them
// There can be, so far, no more than 121, so just return that
//  number for each
//   for(row = 0; row < 11; row++) {
//     for(column = 0; column < 11; column++){
//       row_north = row - 1;
//       if(row_north < 0) { row_north = 10;}
//       row_south = row + 1;
//       if(row_south > 10) {row_south = 0;}
//       col_east = column + 1;
//       if(col_east > 10) {col_east = 0;}
//       col_west = column - 1;
//       if(col_west < 0) {col_west = 10;}

//       int check = seg(row, column);

//       if( seg(row_north, column) == check) {
//         if( seg(row, col_east) == check) {
//           if( seg(row, col_west) == check) {
//             if( seg(row_south, column) == check) {
//               if( seg(row_north, col_east) == check) {
//                 if( seg(row_north, col_west) == check) {
//                   if( seg(row_south, col_east) == check) {
//                     if( seg(row_south, col_west) == check) {
//                       interior_counts(row, column) = check;
//                     }
//                   }
//                 }
//               }
//             }
//           }
//         }
//       }
//     }
//   }
//   return interior_counts;
// }
