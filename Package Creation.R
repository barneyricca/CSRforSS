# Package Creation Script ####
# Just copy, change the names, etc., and use it!

# This follows from:
#  https://r-mageddon.netlify.app/post/writing-an-r-package-from-scratch/

# For now, this is just a "hold your nose and go" thing,
#  but I might document it better later. See the above for
#  some documentation, or look the following for a more detailed story.
#  https://www.r-bloggers.com/easymts-my-first-r-package-story-and-results/
# Some of the stuff below was taken from this site and from the sites that
#  it references.

# Additional help at:
#  https://www.hvitfeldt.me/blog/usethis-workflow-for-package-development/


##  Set-up  ####
#  - You need a free GitHub account
#  - Do NOT create a GitHub repository before starting this
#  - Log out of GitHub
#  - Must have Git/GitHub installed on local computer
#  - If necessary, install some packages.
# Do NOT install roxygen if it is already installed; you'll probably get
#  a never-ending set of error messages.
{
  c("devtools", "usethis") -> package_names

  for (package_name in package_names) {
    if (!is.element(package_name, installed.packages()[, 1])) {
      install.packages(package_name, repos = "http://lib.stat.cmu.edu/R/CRAN")
    }
    library(package_name, character.only = TRUE,
            quietly = TRUE, verbose = FALSE)
  }
  rm(list = c("package_names", "package_name"))

  if(!is.element("available", installed.packages()[,1])) {
    devtools::install_github("r-lib/available")
  }

  if(!is.element("roxygen2", installed.packages()[,1])) {
    devtools::install_github("klutometis/roxygen")
  }
}
#  - In a clean RStudio session, do the remainder of this script. (I.e., restart RStudio at this point.)

## Step 0 ####
# Set up: You'll need to change these to appropriate things.
#  You'll also need to change the token in the use_github() command below,
#  but you can't get that yet.
"Bernard Ricca" -> licenseName
"barneyricca" -> gitHubUserName
"barneyricca@gmail.com" -> gitHubEmail
"/Users/Barney/Google Drive/Ricca - Complexity Research Methods/" -> pathName # Where to create the package
                                     # This is the OSX (Mac) version
"complexR" -> packageName            # Name to create


## Step 0.1 ####
available::available(packageName)    # Check to see if the packageName is
                                     #  available! This can take some time,
                                     #  and the opened web-pages must be examined
                                     #  manually.
# If all is okay, proceed to Step 1.

## Step 1 ####
usethis::create_package(paste0(pathName,packageName))
# This will open up a new R session, so you'll have to reopen this file to continue,
#  and re-run Step 0. (This file, if you change and save it, will appear in the
#  recent files list.)

## Step 2 ####
# Add functions
"complexR" -> initialScript           # Initial R file to create
usethis::use_git_config(user.name = gitHubUserName,
                        user.email = gitHubEmail)
usethis::use_r(initialScript)
#
# Now, add a function to the file above.
# Functions to add (check package:tseriesChaos, too!):
# bootNetworkSig()
# butterfly()
# cusp()
# fold()
# makeMarkov()
# networkEntropy()
# orbde()
# swallowtail()
# transitionNetwork()
# windowedEntropy()
#
# For intialization purposes, each function has this form:
#  test_add <- function (a,b) {
#    return(a + b)
#  }
#
devtools::load_all()       # Load the functions
# Try them! All should return 5
butterfly(2,3)
cusp(2,3)
fold(2,3)
networkEntropy(2,3)
orbde(2,3)
swallowtail(2,3)
transitionNetwork(2,3)
windowedEntropy(2,3)

## Step 3 ####
# See steps 4, 5, and 6 at:
#  https://corinnovation.com/2019/04/26/creating-an-r-package/
devtools::document()
# Now, edit the DESCRIPTION (see the Files tab in the Files, Plots, Packaages, Help, Viewer pane)
# Now, edit the documentation in the initialScript file. For example, I did the
#  following. It is important to include the #' at the beginning of each
#  line, and not just the '

#' test_add
#'
#' This function returns the sum of two numbers.
#' @param <a> first addend
#' @param <b> second addend
#' @keywords add
#' @export
#' @examples
#' test_add(2,3)

# Every time you add documentation to the script, you must rerun
#  devtools::document() to get the results
devtools::document()

# Now, install the package we're working on
devtools::install()

# Voila! A help file exists for the function.
?cusp

# You'll have to do Step 3 for every function. (Well, you could do
#  all the edits at once in the script, and just go from there.)

## Step 4 ####
usethis::use_mit_license(licenseName)
#
# Edit the Description file (see the Files tab)
#

## Step 5 ####

usethis::use_package_doc()
devtools::document()
#
# Edit documentation
# See: http://r-pkgs.had.co.nz/description.html
#
# Then:
devtools::install()

## Step 6 ####
# Git repo!

# This gives you a chance to commit things. You should use it whenever there needs to
#  be something committed.
usethis::use_git()
# Accept the defaults on the next one. When you submit this,
#  you will get a Git token. Copy that for use in the use_github() command
usethis::browse_github_pat()  # Get a Git token (for the next command)

# Be sure to change the auth_token in the next one!
usethis::use_github(protocol = "https",
                    auth_token = "4971d0c7f53afb58a868a90bc1bd68a4d2d9a982")

## Step 7 ####
# Add a Readme. By using .Rmd rather than .md, R code can be
#  embedded, and it can be knit to .md. If you don't need to embed R code,
#  then just use_readme_md() will work.
usethis::use_readme_rmd()
#usethis::use_readme_md()


# Now, anyone can add your package by doing this:
#  devtools::install_github("barneyricca/complexR")


# Now, do all the scripting, etc. and then do the next stuff.

## Going further ####
#  Build | Clean and Rebuild  (from menu)

usethis::use_testthat()
# Call 'use_test()' to initialize a basic test file and open
#  it for editing

usethis::use_vignette('catastRophe')
# Modify 'vignettes/CSRforSS.Rmd'
# Some thoughts on this are at:
# https://qualityandinnovation.com/2019/10/13/easymts-r-package-quick-solver-for-mahalanobis-taguchi-system-mts-problems/

usethis::use_citation()
# Modify 'inst/CITATION'

# Add a logo
# Something like this:
usethis::use_logo("C:/Users/Tom/Desktop/dogs_hex.png")

## Add Dependencies  ####
# This is an example; add each dependency this way:
use_package('data.table')
# Refer to functions with data.table::function_name()


## Add data ####
# The next stuff is from:
# https://qualityandinnovation.com/2019/10/13/my-first-r-package-part-3/
#
# I want to include two files, one data frame containing 50 observations
# of a healthy group with 5 predictors each, and another data frame
# containing 15 observations from an abnormal or unhealthy group (also
# with 5 predictors). I made sure the two CSV files I wanted to add to
# the package were in my working directory first by using dir().

# > use_data_raw()
# ✔ Creating 'data-raw/'
# ✔ Adding '^data-raw$' to '.Rbuildignore'
# ✔ Writing 'data-raw/DATASET.R'
# ● Modify 'data-raw/DATASET.R'
# ● Finish the data preparation script in 'data-raw/DATASET.R'
# ● Use `usethis::use_data()` to add prepared data to package

# > mtsdata1 <- read.csv("MTS-Abnormal.csv") %>% mutate(abnormal=1)
# > usethis::use_data(mtsdata1)
# ✔ Creating 'data/'
# ✔ Saving 'mtsdata1' to 'data/mtsdata1.rda'

# > mtsdata2 <- read.csv("MTS-Normal.csv") %>% mutate(normal=1)
# > usethis::use_data(mtsdata2)
# ✔ Saving 'mtsdata2' to 'data/mtsdata2.rda'
# Magically, this added my two files (in .rds format) into my /data
# directory. (Now, though, I don’t know why the /data-raw directory is
# there… maybe we’ll figure that out later.) I decided it was time to
# commit these to my repository again.



## Checking for CRAN ####
# Release Issue Checklist #
use_release_issue(version = NULL)

# More stuff for pre-check
# https://www.r-bloggers.com/everything-you-should-know-about-winbuilder/

# More testing?
# https://cloud.r-project.org/web/packages/tinytest/index.html

# More help
# https://www.r-bloggers.com/%f0%9f%93%a6-managing-dependencies-in-packages/
# (also at:)
# https://irudnyts.github.io//managing-dependencies-in-packages/

# Free checks at R-hub:
# https://r-hub.github.io/rhub/index.html

