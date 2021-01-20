#### Resampling
#
# All resampling techniques take original data, and randomly generate
#  a population of samples from the original data, to get some results.
#  How the population of samples are generated is where most of the
#  differences lie.

############## Permutation test ##############
# The main idea:
#  1. Put the (N1 + N2) results from two groups into a single big bin
#  2. Randomly sample N1 from the big bin into a first group.
#  3. Put the rest into a second group.
#  4. Compare means
#  5. Repeat 1-4 many, many times to get a normal distribution
#  6. Compare the original data to the resulting distribution to get
#     a p-value
# (For the mathematically particular, this should actually be a 
#   "combination" test rather than a "permutation" test)

# Example from PLOS One: "Beer consumption increases human attractiveness
#  to malaria mosquitoes" (DOI: 10.1371/journal.pone.0009546). First seen
#  at:
#  http://sas-and-r.blogspot.com/2014/11/example-201413-statistics-doesnt-have.html

beer = c(27, 19, 20, 20, 23, 17, 21, 24, 31, 26, 28, 20, 27, 19, 25, 31, 24, 28, 24, 29, 21, 21, 18, 27, 20)
water = c(21, 19, 13, 22, 15, 22, 15, 22, 20, 12, 24, 24, 21, 19, 18, 16, 23, 20)

ds = data.frame(y = c(beer, water), 
                x = c(rep("beer", length(beer)), rep("water", length(water))))

# With mosaic...skip for now
# require(mosaic)
# obsdiff = compareMean(y ~ x, data=ds)
# nulldist = do(999)*compareMean(y ~ shuffle(x), data=ds)
# histogram(~ result, xlab="permutation differences", data=nulldist)
# ladd(panel.abline(v=obsdiff, col="red", lwd=2))

# > obsdiff
# [1] -4.377778
# > tally(~ abs(result) > abs(obsdiff), format="percent", data=nulldist)

# TRUE FALSE 
# 0.1  99.9 

# Set up the big bin
alldata <- c(beer, water)
labels <- c(rep("beer", length(beer)), rep("water", length(water)))
obsdiff <- mean(alldata[labels=="beer"]) - mean(alldata[labels=="water"])

# Do the permutation tests: Version 1
N <- 100000   # Number of permutations to test
resample_diff <- vector("numeric", N)

# This function only works for the "beer" and "Water" data.
#  It may have to change a bit for general use
resample1 <- function(dat, lab, N) {
  resample_diff <- vector("numeric", N)
  for(i in 1:N) {
    resample_labels <- sample(lab)
    resample_diff[i] <- mean(dat[resample_labels=="beer"]) - 
      mean(dat[resample_labels=="water"])
  }
  return(resample_diff)
}

resample_diff <- resample1(alldata, labels, N)
hist(resample_diff)

# Do the permutation tests: Version 2
resamp_means = function(data, labs){
  resample_labels = sample(labs)
  resample_diff = mean(data[resample_labels=="beer"]) - 
    mean(data[resample_labels=="water"])
  return(resample_diff)
}
nulldist = replicate(N,resamp_means(alldata,labels))


hist(nulldist, col="cyan")
abline(v = obsdiff, col = "red")
# The p-value is obtained by counting the proportion of statistics 
#  (including the actual observed difference) among greater than or equal 
#  to the observed statistic:
  
alldiffs = c(obsdiff,nulldist)
p = sum(abs(alldiffs >= obsdiff)/ N)

############ Bootstrap ##############
#
# "[T]he essence of bootstrap: resampling the observed data with 
#  replacement and computing the statistic of interest...many times on the 
#  resampled data to get a distribution of the statistic of interest. 
#  This distribution of the statistic of interest can then be used to 
#  compute, for example, confidence intervals."
#
