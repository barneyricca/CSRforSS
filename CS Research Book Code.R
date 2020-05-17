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
# Analyses included here:
# Chapter 1:
# Chapter 2:
# Chapter 3:
#   Fitting distributions
#   ks.test() and the like
# Chapter 4:
#   inference via resampling
# Chapter 5:
#   ARIMA models
#   windowed entropy approaches
# Chapter 6:
#   orbital decomposition
#   Markov models
#   recurrence quantification analysis
#   cross-recurrence quantification analysis
# Chapter 7:
#   cusp catastrophe
# Chapter 8:
#   persistent homology
#   clustering
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

# Reprodcibility
set_here()      # This helps with folder structure when shifting computers

# Usually, seeds for random number generation would be set here. They are, for pedagogical
#  reasons, set with each snippet of code. The general process would be this:
#  as.integer(Sys.time()) %% 100000 -> time_seed    # Get a seed based on the current time/date
#  paste("Seed:", time_seed)                        # display the seed
#  set.seed(time_seed)                              # Set the seed

# Stuff I prefer
options(show.signif.stars = FALSE) # Do not show significance stars! Showing
                                   #  significance stars encourages the
                                   #  conflation of significant and effect
                                   #  size.
options(digits = 5)                # Round to 5 digits


# ****************************
# End of Setup
# ***************************************************************



# ***************************************************************
## Function Definitions ####
# For ease of use, some code is collected in functions.
#  Those functions(and expected parameters) are:
#   Chapter 5:
#    census() - a helper function for win_ent()
#    win_ent() - returns the windowed entropy of a data sequence
#    OD(sequence of data) - performs orbital decomposition
# ****************************
# ****************************

OD <- function(data_seq) {
  # Returns a table  with the final orbital decomposition analysis:
  #  - String length, C
  #  - Trace of the transition matrix for the string length C, M^C
  #  - Topological entropy, H_T (this is better than H_S for these purposes; use minimum of H_T)
  #  - Shannon entropy, H_S
  #  - D_L - dimensionality via Lyapunov
  #  - chi-square
  #  - degrees of freedom, df
  #  - p-value, p
  #  - Number of possible strengths of length C, N*
  #  - phi-square (phi-sq = chi-sq/N*)
  # This all follows Guastello, Chapter 21, in Guastello & Gregson, 2011, _Nonlinear
  #   Dynamical Systems Analysis for the Behavioral Sciences Using Real Data_.

  # First, recast everything into single character codes:
  unique(data_seq) -> uni_seq
  length(uni_seq) -> n_codes
  if(n_codes > 26) {
    cat("Cannot do more than 26 codes\n")
    return(NULL)
  }
  LETTERS[1:n_codes] -> uni_rep
  uni_seq -> names(uni_rep)
  uni_rep[data_seq] -> data_seq

  # Begin processing data
  length(data_seq) -> seq_len
  1 -> C
  length(data_seq) -> N_star

  TRUE -> continue
  0 -> recurs
  data_seq -> keep1
  data_seq -> seq1

  shift(seq1, 1) -> seq2
  length(which(seq1 == seq2)) -> recurs

  unique(seq1[which(seq1==seq2)]) -> repeats
  length(repeats) -> trMC
  log2(trMC)/C -> H_T
  exp(H_T) -> D_L
  length(unique(data_seq)) - 1 -> dof # For the singlets, this is true;
                                      #  after this, it is the number of unique
                                      #  repeated codes

  table(seq1) -> F_obs_tab
  F_obs_tab / seq_len -> p_obs_keep   # This gets used for calculations later
  rep(seq_len / length(unique(data_seq)),
      length(unique(data_seq))) -> F_exp
  F_obs_tab / seq_len -> p_obs_tab
  -1 * sum(p_obs_tab * log(p_obs_tab)) -> H_S
  sum(F_obs_tab * log(F_obs_tab / F_exp)) * 2 -> chi_sq
  dchisq(chi_sq, dof) -> p
  chi_sq / N_star[C] -> phi_sq

  while(continue) {
    # Identify all recurrences of length C
    N_star[C] - 1 -> N_star[C + 1]
    for(i in 1:N_star[C + 1]) {
      paste(seq1[i], seq1[i+1], sep = "") -> seq1[i]
    }
    seq1[-N_star[C]] -> seq1
    shift(seq1, (C+1)) -> seq2
    length(which(seq1 == seq2)) -> recurs[C + 1]

    C + 1 -> C
    unique(seq1[which(seq1==seq2)]) -> repeats   # Which codes are repeated?
    length(repeats) -> trMC[C]                   # How many repeated codes are there?
    log2(trMC[C])/C -> H_T[C]                    # Topological entropy
    exp(H_T[C]) -> D_L[C]                        # Lyapunov dimension

    # The number of unique repeated codes:
    table(seq1) -> repeated_codes -> F_obs_tab
    length(repeated_codes[which(repeated_codes > 1)]) -> dof[C]

    F_obs_tab[which(F_obs_tab > 1)] -> F_rep_obs_tab  # Repeated codes

    # Calculate H_S:
    length(F_obs_tab) -
      length(F_rep_obs_tab) -> singlets               # Number of singlets; for H_S
    F_rep_obs_tab / N_star[C] -> p_obs_tab            # For H_S
    -1 * sum(p_obs_tab * log(p_obs_tab)) -> H_S[C]    # Multiple repeats
    H_S[C] + singlets *
      (log(N_star[C]) / N_star[C]) -> H_S[C]          # Singlets

    # Frequency expected
    length(F_rep_obs_tab) -> n_rep
    rep(1, n_rep + 1) -> F_exp                        # For all the repeats and a
    0 -> F_exp[n_rep+1]                               #  slot for the singlets

    for(i in 1:n_rep) {
      for(j in 1:C) {
        substr(names(F_rep_obs_tab)[i], j, j) -> fn
        p_obs_keep[fn] * F_exp[i] -> F_exp[i]
      }
      F_exp[i] * N_star[C] -> F_exp[i]
    }
    N_star[C] - sum(F_exp) -> F_exp[n_rep+1]

    #    F_reps_obs_tab / N_star[C] -> F_obs
    N_star[C] - sum(F_rep_obs_tab) -> F_rep_obs_tab[length(F_rep_obs_tab) + 1]

    sum(F_rep_obs_tab * log(F_rep_obs_tab / F_exp)) * 2 -> chi_sq[C]
    dchisq(chi_sq[C], dof[C]) -> p[C]
    chi_sq[C] / N_star[C] -> phi_sq[C]

    if(recurs[C] == 0) {   # Keep going or not?
      FALSE -> continue
    }
  }

  data.table("C" = 1:C,
             "trMC" = trMC,
             "H_T" = H_T,
             "D_L" = D_L,
             "chi-squared" = chi_sq,
             "df" = dof,
             "p-value" = p,
             "Nstar" = N_star,
             "phi-squared" = phi_sq,
             "H_S" = H_S) -> OD_tab
  OD_tab[-nrow(OD_tab),] -> OD_tab     # Fix this some day: I dunno why this stops
  #  after the second zero rather than the first
  return(OD_tab)
}

census <- function(data_vec, codes) {
  rep(0, length(codes)) -> census_vec
  codes -> names(census_vec)

  for(i in seq_along(data_vec)) {
    census_vec[data_vec[i]] + 1 -> census_vec[data_vec[i]]
  }
  return(census_vec)
}

win_ent <- function(data_seq, parts = 16, lag_max = 50) {
  # Follows Wiltshire, Fiore, & Butner (2017)
  # Returns a data frame (or NULL, if failed) of two time series:
  #  1. A (numeric) windowed entropy
  #  2. A marker of maxima (TRUE for a maxima, FALSE for not)
  #
  if(mode(data_seq) == "integer") {
    as.character(data_seq) -> data_seq
  }

  if(mode(data_seq) == "numeric") {
    max(data_seq) -> max_data
    min(data_seq) -> min_data
    as.character(1:parts) -> codes
    ceiling( parts *
               (data_seq - min_data) /
               (max_data - min_data)) -> partition
    1 -> partition[which(partition == 0)]  # partition goes 0:parts, not 1:parts
    codes[partition] -> data_seq
  }

  if(mode(data_seq) != "character") {
    return(NULL)
  }

  require(data.table)
  unique(data_seq) -> codes
  length(codes) -> n_codes
  table(data_seq) -> code_counts
  1:n_codes -> nums
  codes -> names(nums)

  code_counts / sum(code_counts) -> code_probs

  # The appropriate lag, according to Fraser & Swinney, 1986, is at
  #  the first minimum of mutual information.
  which.min(mutual(unname(nums[data_seq]),
                   partitions = n_codes,       # Arbitrary partitions needed                                                          #  for continuous data; any
                   #  partition >= n_codes works
                   #  for categorical data.
                   lag.max = lag_max,
                   plot = FALSE)) -> wind_sz

  # Calculate the windowed entropy, ent[]:
  # This is different than what Wiltshire, Fiore, & Butner (2017) do. However,
  #  I think that they overestimate the issue (and they probably would think
  #  that I underestimate the issue.) They use the natural log, and only use
  #  the probabilities within the window.
  length(data_seq) -> seq_len
  rep(0.0, length = (seq_len - wind_sz + 1)) -> ent

  for(i in 1:(length(ent))) {
    for(j in 0:(wind_sz-1)) {
      ent[i] - (code_probs[data_seq[i+j]] *
                  log2(code_probs[data_seq[i+j]])) -> ent[i]
    }
  }

  # Find the appropriate maxima:
  # Pass 1 - find all maxima:
  vector("logical", length(ent)) -> maxes  # All entries are FALSE by default
  for (i in 2:(length(ent)-1)) {
    if(ent[i] >= ent[i-1] &&
       ent[i] >= ent[i+1]) {
      TRUE -> maxes[i]
    }
  }

  # Pass 2 - which maxima indicate code probability changes from before to after a
  #  maximum located in Pass 1
  1 -> i_old
  length(data_seq) -> n_max
  for(i in 1:length(maxes)) {
    if(maxes[i]) {
      0 -> cont_tab
      if((i+1) < (n_max - n_codes - wind_sz)) {       # Skip if too close to the end
        census(data_seq[i_old:i], codes) -> before
        census(data_seq[(i+1):n_max], codes) -> after
        rbind(before,after) -> cont_tab

        # I prefer the fisher.test() for the next one, but it sometimes
        #  runs out of memory. The chi-squared sometimes fails for really
        #  small cell values (and zeroes); hence the is.na().
        #  The chisq.test() is also very biased, so it will generate warnings.
        if(sum(cont_tab) > 0 ) {
          chisq.test(cont_tab) -> p_temp
          if(is.na(p_temp$p.value)) {
            FALSE -> maxes[i]
          } else {
            if( (chisq.test(cont_tab))$p.value > 0.10) {  # 0.1 is arbitrary
              FALSE -> maxes[i]  # No change
            } else {
              i -> i_old         # Only include info from this sub-segment
            }
          }
        } else {
          FALSE -> maxes[i]
        }
      }
    }
  }
  # In principle, this refinement should be iterated until it converges. However,
  #  as we know this approach is biased, we won't try to refine something that
  #  will still be a bit biased.

  data.frame("Windowed_Entropy" = ent,    # Entropy
             "Maxima" = maxes) -> ent_df  # True if a local maximum of entropy
  return(ent_df)
}



# ****************************
# End of Function definitions
# ***************************************************************


# ***************************************************************
# ***************************************************************
# ***************************************************************
#
##           CHAPTER 1 ####
#
# Traffic (NetLogo) simulation
# Conway's Game of Life (There's one in NetLogo)
# Schelling's (Netlogo) model
#
# Langton's Ant
#
#
# ***************************************************************
# ***************************************************************
# ***************************************************************


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
# ***************************************************************
# ***************************************************************
#
##           CHAPTER 2 ####
#
# Rebellion (NetLogo) simulation
# Conway's Game of Life
# Schelling's (Netlogo) model
#
#
# ***************************************************************
# ***************************************************************
# ***************************************************************




# ***************************************************************
# ***************************************************************
# ***************************************************************
#
##           CHAPTER 3 ####
#
# ***************************************************************
# ***************************************************************
# ***************************************************************

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
# Figure 3.2
#


# ***************************************************************
# Figure 3.3
#



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
# Figure 3.8
#

# ***************************************************************
# Figure 3.9
#

# ***************************************************************
# Figure 3.10
#

# ***************************************************************
# Figure 3.11
#

# ***************************************************************
# Figure 3.12
#

# ***************************************************************
# Figure 3.13
#

# ***************************************************************
# Figure 3.14
#




# ***************************************************************
# ***************************************************************
# ***************************************************************
#
##           CHAPTER 4 ####
#
# ***************************************************************
# ***************************************************************
# ***************************************************************




# ***************************************************************
# Figure 4.1
# This replicates Figure 1 from Hoffman et al., (2018).
#  The exact values will be different here from Hoffman et al., because of the differences
#  in pseudo-random number generation (e.g., different "seeds").
# ****************************
set.seed(630)

generate_scores <- function(x, y) {
  # Here are the parameters used by Hoffman:
  0.5 -> mu_int_young
  2.5 -> mu_int_old
  0.1 -> mu_slope
  0 -> mu_error

  0.1 -> sd_int_young
  0.1 -> sd_int_old
  0.02 -> sd_slope
  0.2 -> sd_error

  ifelse(y < x,
         rnorm(1, mu_int_young, sd_int_young) + rnorm(1, mu_slope, sd_slope) * y,
         rnorm(1, mu_int_old, sd_int_old) + rnorm(1, mu_slope, sd_slope) * y) -> scores
  scores + rnorm(length(scores), mu_error, sd_error) -> scores
  return(scores)
}

# Generate data
30 -> N                                           # 30 students
6 -> mu_mature                                    # For age at transition
0.5 -> sd_mature                                  # For age at transition
rnorm(N, mu_mature, sd_mature) -> transitions     # Age of transition
seq(from = 4, to = 8, by = 0.02) -> ages

sapply(transitions,            # sapply() takes the first argument (variables)
       generate_scores,        #  and uses each on in sequence in the second argument
                               #  (the function)
       y = ages) -> score      # and adds additional arguments. Then put the result
                               #  for each student into the score variable
rowSums(score) / N -> mean_score

# Make the data plot
{
  plot(ages, score[,1],
       type = 'l',
       col = "grey",
       xlab = "Age",
       ylab = "Score")
  for(i in 2:N) {
    lines(ages, score[,i],
          col = "grey",
          lwd = 0.6)
  }
  lines(ages, mean_score,   # Add the mean of the data at each age
        col = "black",
        lwd = 0.9)
}


# ****************************
# End of Figure 4.1
# ***************************************************************


# ***************************************************************
# Figure 4.2
#
# This is identical to Figure 4.1, except for the use of the median_score
#  in place of the mean_score
# ****************************

set.seed(630)

generate_scores <- function(x, y) {
  # Here are the parameters used by Hoffman:
  0.5 -> mu_int_young
  2.5 -> mu_int_old
  0.1 -> mu_slope
  0 -> mu_error

  0.1 -> sd_int_young
  0.1 -> sd_int_old
  0.02 -> sd_slope
  0.2 -> sd_error

  ifelse(y < x,
         rnorm(1, mu_int_young, sd_int_young) + rnorm(1, mu_slope, sd_slope) * y,
         rnorm(1, mu_int_old, sd_int_old) + rnorm(1, mu_slope, sd_slope) * y) -> scores
  scores + rnorm(length(scores), mu_error, sd_error) -> scores
  return(scores)
}

# Generate data
30 -> N                                           # 30 students
6 -> mu_mature                                    # For age at transition
0.5 -> sd_mature                                  # For age at transition
rnorm(N, mu_mature, sd_mature) -> transitions     # Age of transition
seq(from = 4, to = 8, by = 0.02) -> ages

sapply(transitions,            # sapply() takes the first argument (variables)
       generate_scores,        #  and uses each on in sequence in the second argument
       #  (the function)
       y = ages) -> score      # and adds additional arguments. Then put the result
#  for each student into the score variable

apply(score, 1, median) -> median_score

# Make the data plot
{
  plot(ages, score[,1],
       type = 'l',
       col = "grey",
       xlab = "Age",
       ylab = "Score")
  for(i in 2:N) {
    lines(ages, score[,i],
          col = "grey",
          lwd = 0.6)
  }
  lines(ages, median_score,   # Add the mean of the data at each age
        col = "black",
        lwd = 0.9)
}



# ****************************
# End of Figure 4.2
# ***************************************************************

# ***************************************************************
# ***************************************************************
# ***************************************************************
#
##           CHAPTER 5
#
# ***************************************************************
# ***************************************************************
# ***************************************************************


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

# Figure 5.6
1:10 -> x
0.5 * x * x + 1 -> y
set.seed(42)
rnorm(10, 0, 2) + y -> y
plot(x, y,
     pch = 16,
     cex = 0.7,
     main = "Fictitious Quadratic Data",
     xlim = c(0, 11),
     ylim = c(0, 50))
lm(y ~ x) -> yquad.lm
plot(x, yquad.lm$residuals,
     main = "Residual Plot",
     ylab = "residual",
     pch = 16,
     cex = 0.7)

x * x -> x2
lm(y ~ x2) -> yquad2_lm

# Figure 5.8
{
plot(x, y,
     pch = 16,
     cex = 0.7,
     main = "Fictitious Quadratic Data",
     xlim = c(0, 11),
     ylim = c(0, 50))
lines(x, yquad2_lm$fitted.values,
      col = "blue")
lines(x, yquad.lm$fitted.values,
      col = "darkgreen")
}

# Figure 5.9
plot(x, yquad2_lm$residuals,
     main = "Residual Plot",
     ylab = "residual",
     pch = 16,
     cex = 0.7)

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
