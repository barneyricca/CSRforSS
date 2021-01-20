{         # Put the cursor on this line, and run it; because
  #  of the curly braces, all the code within
  #  the braces will be run.
  c("bmp",          # Read .BMP files
    "data.table",   # fread() for fast data input
    "dplyr",        # Data wrangling
    "dtplyr",       # For faster tidyverse data wrangling
    "dtw",          # Dynamic time warping
    "gridExtra",    # For better table printing
    "here",         # To assist with folder structure
    "igraph",       #
    "raster",       # For image processing and coverings
    "R.matlab",     # To read matlab data from Zhang et al. (2020)
    "Rcpp",         # For faster processing of Betti numbers
    "TDA") ->
    package_names

  for (package_name in package_names) {
    if (!is.element(package_name, installed.packages()[, 1])) {
      install.packages(package_name, repos = "http://lib.stat.cmu.edu/R/CRAN")
    }
    library(package_name,
            character.only = TRUE,
            quietly = TRUE,
            verbose = FALSE
    )
  }
  rm(list = c("package_names", "package_name"))
}

fread(file = "/Users/Barney/Dropbox/Dynamic analysis of gait and kicking/Beginning Walk Full Term.csv",
      header = TRUE) -> walk1_df
fread(file = "/Users/Barney/Dropbox/Dynamic analysis of gait and kicking/Beginning Walk Pre Term.csv",
      header = TRUE) -> walk2_df
fread(file = "/Users/Barney/Dropbox/Dynamic analysis of gait and kicking/Experienced Walk Full Term.csv",
      header = TRUE) -> walk3_df
fread(file = "/Users/Barney/Dropbox/Dynamic analysis of gait and kicking/Experienced Walk Pre Term.csv",
      header = TRUE) -> walk4_df

names(walk1_df)
# I need some help on a few of the names, perhaps, just for my purposes:
#  MetheadVel
#  lateralMalleoleiVel
#  footRotationAngle
#  hipAbdAdd
#  hipFlex
# And...
#  What are the directions for X, Y, and Z? (E.g., front-back,
#   left-right, and up-down?)
#  What are the reference origins? (E.g., 0 degrees is...)

# As an example
plot(walk1_df$Time, walk1_df$ankleAngle,
     type = 'l')

# OK, can do some DTW on these things.
dtw(walk1_df$ankleAngle, walk3_df$ankleAngle) -> dtw1
