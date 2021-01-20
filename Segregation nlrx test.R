library(nlrx)
library(future)   # For parallelizing

# Look at advantages of nlrx over RNetLogo. The biggest will be simply the
#  reproducibility and the "tidy" nature of the data that comes back.
# Generally, it appears that you will have to clean up your own temp
#  directory, lest it get really bloated. (I could put in an R script
#  to do that at the end of this script; file.remove() will do it.
#  Be sure to delete the file and not the directory! You'll need the
#   directory.

# RNetLogo is relatively straightforward, provided you have the correct
#  version of Java. look at the web-site for an example of it.

# There's a lot more here (see Salecki et al. 2019) but this isn't
#  a NetLogo tutorial nor is it a course on ABM. Hence, we'll just
#  do the basics.

# The next four commands may be necessary on a Mac. (I needed them!)
# THese solve issues with the temp directory and with finding Java.
unixtools::set.tempdir("/Users/barneyricca/temp")
# Can the previous line be replaced with:
# unixtools::set.tempdir("~/temp")
# so that it will be the same on all machines? (Or at least, all Macs....)

system("echo $JAVA_HOME") # Probably returns a blank. If so, run the next two.
# Go to the Terminal tab in the Console pane, and type:
# java -version
# to make sure you have Java (and find what version you have). Then:
# /usr/libexec/java_home
# to find the information needed for the next line. Use whatever is returned
#  in place of the quoted text in the next line.
Sys.setenv(JAVA_HOME = "/Library/Java/JavaVirtualMachines/jdk-15.0.1.jdk/Contents/Home")
system("echo $JAVA_HOME")  # Hopefully returns what you put in
# End of fixing Java and tempdir() problems.

# 1: Create an nl object

netlogopath <- file.path("/Applications/NetLogo 6.1.1")
modelspath <- file.path("models/Sample Models/Social Science")
modelpath <- file.path(netlogopath, "Segregation.nlogo")
outpath <- file.path("/Users/barneyricca/Google Drive")

nl <- nl(nlversion = "6.1.1",
         nlpath = netlogopath,
         modelpath = modelpath,
         jvmmem = 1024)

# 2. Attach an experiment
# Get metrics, metrics.patches, constants, metrics.turtles, and variables from
#  the NetLogo code. Show how each of these connects:
#  metrics - globals
#  metrics.patches - patches-own
#  constants - on interface
#  metrics.turtles - turtles own
#  metrics.links - links own
#  variables - on interface
# Look @ ?nlrx::experiment
# also: a slot "repetition" works, but must set seed manually. Probably it 
#  is best for replication to set repetition = 1 and increase nseeds.
#  It is implied, but never said, that nseeds =  will always send the
#   same seeds. Should run things twice to see what happens.
# repetition = 1 is the default. Use that (and use nseeds to set
#  the number of attempts.)
#
# See what the variables list really produces
# Are there any other defaults we can eliminate?
# idsetup and idgo are pretty standard, but not always. Look at the
#  interface.
nl@experiment <- experiment(expname="Segregation1",  # Experiment name
                            outpath=outpath,         # Output directory
                            repetition=4,            # How many repetitions
                            tickmetrics="true",      # record @ each tick?
                            idsetup="setup",         # Name of setup procedure
                            idgo="go",               # Name of go procedure
                            runtime=5000,            # Number of ticks
                            # metrics.patches = c()   #
                            # constants = list()      # If one wants to set these and leave them
                            metrics=c("percent-unhappy"),
                            metrics.turtles=list("turtles"= c("similar-nearby")),
                            variables = list('density' = list(min=30, max=60, 
                                                              qfun="qunif"),
                                             '%-similar-wanted' = list(min=30, 
                                                                       max=80, 
                                                                       qfun="qunif")))

# 3. Create a design
# A. May need to create an evaluation function.
# B. See https://cran.rstudio.com/web/packages/nlrx/vignettes/furthernotes.html
#    for design types
# We'll just do a simple design, though, so the next is sufficient.
nl@simdesign <- simdesign_lhs(nl=nl,
                              # evalcrit = function(x){},
                              samples=10,
                              nseeds=3,
                              precision=3)

# 4. Run simulations
# package:progress has stuff too

# package:future has the next function
# Just do this, with the split, just in case
#  Want nseeds * split = number of available cores for maximum 
#  efficiency.
plan(multisession)           # Run in parallel on local machine:

# Started at 10:59 a.m.
results <- run_nl_all(nl)
#  results %<-% run_nl_all(nl = nl, split = 4)  # Run on 4 times the number
#                                             # nseed cores.
# Also seed = nl@simdesign@simseeds[1] can be used
# Ended at 11:08 a.m. (with all 8 cores, I presume)

# stopcond(), etc.
# Suggestion: In NetLogo do this:
#to-report green.patches
#  report count patches with [pcolor = green]
#end
# and then "green.patches" can be used in the set up in place of 
#  "report count patches with [pcolor = green]"

save(results, file = "Segregation results 1.RData")
# The file is kinda big...60 MB


load("Segregation results 1.RData")
results -> setsim(nl, "simoutput")
#write_simoutput(nl)


# Investigate output
# Need a lot of work/guidance here.
unnest_simoutput(nl) -> dum1
analyze_nl(nl) -> dum2

# Also:
# tibble::enframe(results) -> nl@simdesign@simoutput

# getsim(nl, "simoutput") %>%
  # dplyr, tidyr stuff.
  # See:
  # https://cran.rstudio.com/web/packages/nlrx/vignettes/abc.html

# https://cran.rstudio.com/web/packages/nlrx/vignettes/sensitivity.html
