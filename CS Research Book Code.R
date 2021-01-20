# ***************************************************************
# This file contains all of the R code for the examples and figures
#  in the book.
#
# ABM models are found elsewhere.
#
#  Make this into a package:
#   https://corinnovation.com/2019/04/26/creating-an-r-package/
#
# What about Xquartz and/or Rtools?
# ***************************************************************


#
# Change this to Rmd files with mother & child:
# The script header for the mother file should be something like this:
# ```{r child = here::here("Script", "Rmds", "Chapter-1.Rmd")}
# ```
#
#

# ***************************************************************
## Setup for all code ####
#
# Run the next set of commands before doing anything else. This set up must be
#  repeated each time you re-open this file.
#
# You do NOT need to understand any of this stuff; just run
#  the commands.
# ****************************

# Install and/or load all necessary packages:
{         # Put the cursor on this line, and run it; because
          #  of the curly braces, all the code within
          #  the braces will be run.
  c("conflicted",
    "cusp",         # Does cusp catastrophe
    "data.table",   # fread() for FAST and flexible data input and shift()
    "dtplyr",       # For faster data wrangling with tidyverse structure
    "here",         # To assist with folder structure
    "igraph",       # For creating and plotting graphs;
                    #   ggraph & tidygraph are probably better.
                    #   If one is not displaying graphs, then
                    #   package:igraph is not necessary.
    "MASS",         # For fitting distributions
    "poweRlaw",     # To fit power law data; MASS::fitdistr() can't
    "pracma",       # package:pracma is used for demonstration
                    #  only; it is not needed for the analysis
                    #  of nonlinear time series.
    "Rcpp",         # For some faster routines
#    "rgl",          # 3-d plotting with rotation
    "TDA",          # Topological Data Analysis
    "tseries",      # Basic time series routines
    "tseriesChaos") -> # Useful for mutual information calculations
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

  conflict_prefer("decompose", "stats")
}


# ***************************************************************
# Figure 1.9
#
rep(c(-1, 0, 1),          # 3 values to repeat: -1, 0, & 1
    each = 16) -> accel   # Repeat each one 16 times; put into variable accel
seq(from = 0.25,          # Lowest number in the sequence
    to = 12,              # Highest number in the sequence
    by = 0.25) -> gap     # Interval between entries; put into variable gap
plot(gap, accel,          # Plot with gap on horizontal & accel on vertical
     pch = 20,            # Use smaller, solid dots
     main = "Car Acceleration vs. Gap",        # Main title
     xlab = "Gap (arbitrary units)",           # Title on x-axis
     ylab = "Acceleration (arbitrary units)")  # Title on y-axis


# ***************************************************************
# Figure 3.1
#

# Simulate the Lorenz system
sim.cont(syst = lorenz.syst,
         start.time = 0,
         end.time = 8000,
         dt = 0.05,
         start.x=c(5,5,5),
         parms=c(10, 28, -8/3),
         obs.fun = function(x) list(x)) -> lorenz_ts

# Convert from lorenz_ts list to 3 time series. There are
#  other ways to do this, but this works
unlist(lorenz_ts) -> temp_ts    # Change from list to named vector
unname(temp_ts[which(names(temp_ts) == "1")]) -> x_ts
unname(temp_ts[which(names(temp_ts) == "2")]) -> y_ts
unname(temp_ts[which(names(temp_ts) == "3")]) -> z_ts

# library(plot3D)
50 -> n
lines3D(x_ts[1:n], y_ts[1:n], z_ts[1:n])

# Figure 3.1
plot(z_ts,
     type = 'l',
     main = "z-component of Lorenz map",
     xlab = "time index",
     ylab = "z")

# Figure 3.2
data.frame("Time" = 1:(length(z_ts)),
           "Z" = z_ts) -> lorenz_z
lm(Z ~ Time, data = lorenz_z) -> z_lm
{
  plot(z_ts,
       type = 'l',
       main = "z-component of Lorenz map",
       xlab = "time index",
       ylab = "z")
  abline(a = z_lm$coef[1],
         b = z_lm$coef[2],
         col = "green")
}

# Figure 3.3
hist(z_lm$residuals,
     breaks = 50,
     main = "Histogram of residuals",
     xlab = "Residual")



# ***************************************************************
# Figure 3.4
#
c(2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6,
  7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9,
  10, 10, 10, 11, 11, 12) -> dice_rolls
hist(dice_rolls,
     main = "Frequency of Results from Two Dice",
     xlab = "total",
     ylab = "frequency",
     xlim = c(1.5, 12.5),
     breaks = seq(from = 1.5,
                  to = 12.5,
                  by = 1))
#
# ***************************************************************



# ***************************************************************
# Figure 3.5
#
R = 1e6
replicate(R,
          sum(sample(1:6, 2, replace = TRUE))) ->
  two_dice_results
ks.test(results,
        "pnorm")
#
# ***************************************************************



# ***************************************************************
# Figure 3.7
# ****************************
set.seed(42)
# First, randomly set up the cups:
sample(rep(c("M","T"), each = 4), 8, replace = FALSE) -> treatment

# Determine how many times to try guessing
1e6 -> N

# The guessing function
guess1 <- function(treat) {
  # Make the guess
  sample(rep(c("M","T"), each = 4), 8, replace = FALSE) -> attempt
  # Compare to the treatment
  if(all(attempt == treat)) {  # Compare the guesses to the treatments
    return(TRUE)   # The guesses were all correct!
  } else {
    return(FALSE)  # At least one guess was wrong
  }
}

replicate(N, guess1(treatment)) -> results
sum(results) / N


# ****************************
# End of Figure 3.7
# ***************************************************************


# Spam data (Poisson)
# From my SJFC spam filter: Number of spam messages caught daily, from 1/23/20
#  through about the next month
c(4, 3, 2, 1, 3, 1, 3, 2, 2, 0, 0, 0, 1, 1, 6, 5, 1, 2, 1, 2, 3, 1, 2,
  3, 1, 0, 1, 2, 1, 3, 3, 2, 2) -> spam
ks.test(spam,         # data to test
        ppois,        # hypothesized distribution
        mean(spam))   # parameter for ppois

library(MASS)
fitdistr(spam, "Poisson")  -> pois_est

pois_est$estimate
pois_est$sd

{
  hist(spam,
       main = "Daily spam probabilities",
       xlab = "Number of spam messages caught",
       ylab = "Probability",
       breaks = seq(       # "breaks =" makes a better looking graph
         from = -0.5,
         to = max(spam) + 0.5,
         by = 1),
       prob = TRUE)
  lines(0:max(spam),
        dpois((0:max(spam)),
              lambda = pois_est$estimate),
        lwd = 2,
        col = "mediumblue")
}

# ***************************************************************
# Chapter 5 ####
#  TraMineR: http://traminer.unige.ch/
# ****************************

# Residual analysis and Figures 5.1-5.4


1:15 -> x
3*x + 4 -> y
set.seed(42)
runif(15, -5, 5) + y -> yunif
rnorm(15, 0, 3) + y -> ynorm
lm(yunif ~ x) -> yunif.lm
lm(ynorm ~ x) -> ynorm.lm
summary(yunif.lm)
summary(ynorm.lm)
hist(yunif.lm$residuals,
     breaks = 10,
     main = "Histogram of residuals",
     xlab = "residual")

summary(yunif.lm)

{
  plot(x, yunif,              # Plot the data
       main = "Dummy Data",   # Title above the plot
       xlab = "X",            # x-axis label
       ylab = "Y",            # y-axis label
       xlim = c(0,16),        # x-axis range
       ylim = c(0,50),        # y-axis range
       pch = 16,              # Use a solid dot
       cex = 0.8)             # make the dots a little smaller
  abline(a = yunif.lm$coef[1],    # Add a line, with intercept a
         b = yunif.lm$coef[2],    #  and slope b
         col = "darkgreen")   #  in a nice dark green color
  text(1, 40,                 # Put some text info on the plot
       labels = paste0("slope = ",
                       round(yunif.lm$coefficients[2], digits = 2)),
       adj = c(0,0))          # Position the lower left at (1,40)
  text(1, 35,                 # Put some more text on the plot
       labels = paste0("intercept = ",
                       round(yunif.lm$coefficients[1], digits = 2)),
       adj = c(0,0))
}



# Figure 5.2
{
  plot(x, yunif.lm$residuals,              # Plot the data
       main = "Residuals",   # Title above the plot
       xlab = "X",            # x-axis label
       ylab = "Residuals",            # y-axis label
       xlim = c(0,16),        # x-axis range
       ylim = c(-5,5),        # y-axis range
       pch = 16,              # Use a solid dot
       cex = 0.8)             # make the dots a little smaller

}

library(lmPerm)
lmp(yunif ~ x, center = FALSE) -> y_perm.lm
summary(y_perm.lm)

# THis is taken from Coghlan, 2014. Get a link to the online version
births <- scan("http://robjhyndman.com/tsdldata/data/nybirths.dat")
birthstimeseries <- ts(births,          # Make a time series of the data
                       frequency=12,    # 12 months per year
                       start=c(1946,1)) # Start with January 1946.
plot.ts(birthstimeseries)
# library(TTR)
decompose(birthstimeseries) -> births_dec
plot(births_dec)

# The decompose() function is similar to the function tslm(), but decompose() has slightly nicer plotting functions for our purposes.

# package:forecast
# auto.arima()
# forecast()
# package:arfima
# Some nice functions (and AIC, too)

tslm(birthstimeseries ~ trend + season) -> births_tslm
plot(forecast(births_tslm, h = 10))


# ***************************************************************
# ***************************************************************
# ***************************************************************
#
##           CHAPTER 6 ####
#
# ***************************************************************
# ***************************************************************
# ***************************************************************




# ***************************************************************
# Figure 6.x
#  Here is a modified logistic that has a significant non-zero slope when plotted
#   as x vs. n, but which largely retains its parabolic shape when plotted x_{n+1} vs. x_n.
# ****************************
3.57 -> R
1e3 -> N
rep(0, N) -> x
1:N -> n
0.4 -> x[1]

1.9 * (1 - R/4) / N -> B

for(i in 1:(N-1)) {
  R*x[i]*(1-x[i]) + B*x[i]*i -> x[i+1]
}

plot(n,x,
     xlim = c(100, 200),
     pch = 16,
     cex = 0.4)

summary(lm(x ~ n))

x[2:N] -> xnp1
x[1:(N-1)] -> xn
plot(xn, xnp1,
     pch = 16,
     cex = .6)

library(pracma)
# The polyfit() yields -3.49 x^2 + 3.59 x + 0.02, rather nicely.
# The generation was done by something like -3.57 x^2 + 3.59 * (1 + epsilon) x
polyfit(xn, xnp1, 2) -> pfit
# ****************************
# End of Figure 6.x
# ***************************************************************


# ***************************************************************
# Figure 6.x
# ****************************
plot(x = 1:(length(logistic_series)),
     y = logistic_series,
     xlab = "n",
     ylab = expression('x'[n]),
     xlim = c(500,700),  # Just look at part of it
     pch = 16,           # Use dots
     cex = 0.5)          # Use smaller dots

shift(logistic_series,
      1) -> xnp1   # "x n plus 1"

plot(xnp1, logistic_series,
     xlab = expression('x'[n]),   # [] puts in a subscript
     ylab = expression('x'[n+1]))

# ****************************
# End of Figure 6.x
# ***************************************************************


# ***************************************************************
# Figure 6.x
# ****************************
length(logistic_series) -> n
logistic_series[1:n-1] -> xn
logistic_series[2:n] -> xnp1
polyfit(xn, xnp1, 2) -> logistic_poly   # Fit a 2nd order polynomial

{
  plot(xn, xnp1,
       col = "green",
       xlab = expression('x'[n]),   # [] puts in a subscript
       ylab = expression('x'[n+1]))
  points(xn,
         logistic_poly[1]*xn*xn + logistic_poly[2]*xn + logistic_poly[3],
         pch = 21,
         cex = 0.4,
         col = "blue")
}
# ****************************
# End of Figure 6.x
# ***************************************************************





# ***************************************************************
# Figure 6.x
# ****************************
library(devtools)
devtools::install_github("ropensci/workloopR")

library(workloopR)

# Example partially taken from:
#  https://ropensci.org/technotes/2019/11/14/workloopr-release/

## import the workloop.ddf file included within workloopR
wl_dat <-read_ddf(system.file("extdata", "workloop.ddf",
                              package = 'workloopR'),
                  phase_from_peak = TRUE)

  names(attributes(wl_dat))

  attr(wl_dat,"stim_table")
  attr(wl_dat,"total_cycles")


  library(ggplot2)

  ## select cycles 3 through 5 using a l0-to-l0 definition
  # There are 6 cycles
  wl_selected <- select_cycles(wl_dat, cycle_def = "lo", keep_cycles = 1:6)

  ## apply a gear ratio correction, run the analysis function,
  ## and then get the full object
  wl_analyzed <- analyze_workloop(wl_selected, GR = 2)

  ## base R work loop plot for the second retained cycle (cycle "b")
  plot(wl_analyzed$cycle_b$Position,
       wl_analyzed$cycle_b$Force,
       xlab = "Position (mm)",
       ylab = "Force (mN)",
       main = "Work loop plot via base R",
       bty = "n",
       tck = 0.02)

  ## now via ggplot
  ggplot(wl_analyzed$cycle_b, aes(x = Position, y = Force)) +
    geom_path(lwd = 2) +
    labs(y = "Force (mN)", x = "Position (mm)") +
    ggtitle("Work loop plot via ggplot2") +
    theme_minimal()

# ****************************
# End of Figure 6.x
# ***************************************************************



# ***************************************************************
# ***************************************************************
# ***************************************************************
#
#           CHAPTER 7
#
# ***************************************************************
# ***************************************************************
# ***************************************************************





# ***************************************************************
# ***************************************************************
# ***************************************************************
#
#           CHAPTER 8
#
# ***************************************************************
# ***************************************************************
# ***************************************************************


# ***************************************************************
# Figure 8.7
# ****************************

graph_from_adjacency_matrix(adj_mat,
                            weighted = TRUE) -> g


# This was set for zooming in on the graph.
# Making good looking graph plots is, as a rule, very hard!
{
  set.seed(42)
  plot(g,
       layout = layout_with_fr,
       vertex.size = counts[V(g)$name]/50,
       vertex.label = V(g)$name,
       edge.width = sqrt(E(g)$weight)/4,
       edge.arrow.size = 0.4,
       edge.curved = 0.5,
       vertex.label.cex = 2,
       vertex.label.dist = 2,
       vertex.label.degree = c(-9*pi/10, 3*pi/4, -pi, 0, pi/2, -2*pi/5, -pi),
       margin = c(0,0,0,0))
}
# ****************************
# End of Figure 8.7
# ***************************************************************
