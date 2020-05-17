prepwork <- function(package_names) {
# Still to do:
#  Check on those packages needing devtools for install
  if(class(packages_names) == "character") {
    if( length(package_names) != 0) {
      c(package_names,
        "data.table",     # Faster data access
        "here",           # To find/store files w/o setwd() and getwd()
        "Rcpp",           # Will move the functions to C++, for speed
        "tidyverse",      # For data manipulation
        "tseriesChaos"    # Some useful routines
      ) -> package_names
      for(package_name in package_names) {
        if(!is.element(package_name, installed.packages()[,1])) {
          install.packages(package_name,
                           repos = "http://cran.mtu.edu/")
        }
        library(package_name, character.only=TRUE,
                quietly=TRUE,verbose=FALSE)
      }
    } else { # If package names is empty
      c("data.table",     # Faster data access
        "here",           # To find/store files w/o setwd() and getwd()
        "Rcpp",           # Will move the functions to C++, for speed
        "tidyverse",      # For data manipulation
        "tseriesChaos"
      ) -> package_names
      for(package_name in package_names) {
        if(!is.element(package_name, installed.packages()[,1])) {
          install.packages(package_name,
                           repos = "http://cran.mtu.edu/")
        }
        library(package_name, character.only=TRUE,
                quietly=TRUE,verbose=FALSE)
      }
    }
  } else { # if package_names is not a character
    c("data.table",     # Faster data access
      "here",           # To find/store files w/o setwd() and getwd()
      "Rcpp",           # Will move the functions to C++, for speed
      "tidyverse",      # For data manipulation
      "tseriesChaos"
    ) -> package_names
    for(package_name in package_names) {
      if(!is.element(package_name, installed.packages()[,1])) {
        install.packages(package_name,
                         repos = "http://cran.mtu.edu/")
      }
      library(package_name, character.only=TRUE,
              quietly=TRUE,verbose=FALSE)
    }
  }
}
