library(nlrx)
library(future)   # For parallelizing

# The next four commands may be necessary on a Mac. (I needed them!)
unixtools::set.tempdir("/Users/barneyricca/temp")

system("echo $JAVA_HOME") # Probably returns a blank. If so, run the next two.

# Go to the Terminal tab in the Console pane, and type:
# java -version
# to find what version of Java you have. Then:
# /usr/libexec/java_home
# to find the information needed for the next line
Sys.setenv(JAVA_HOME = "/Library/Java/JavaVirtualMachines/jdk-15.0.1.jdk/Contents/Home")
system("echo $JAVA_HOME")  # Hopefully returns what you put in

# End of fixing Java and tempdir() problems.

netlogopath <- file.path("/Applications/NetLogo 6.1.1")
modelpath <- file.path(netlogopath, "models/Sample Models/Social Science/Segregation.nlogo")
outpath <- file.path("/Users/barneyricca/Google Drive")

nl <- nl(nlversion = "6.1.1",
         nlpath = netlogopath,
         modelpath = modelpath,
         jvmmem = 1024)

nl@experiment <- experiment(expname="Segregation1",  # Experiment name
                            outpath=outpath,         # Output directory
                            repetition=1,            # How many repetitions
                            tickmetrics="true",      # record @ each tick?
                            idsetup="setup",         # Name of setup procedure
                            idgo="go",               # Name of go procedure
                            runtime=5000,            # Number of ticks
#                            evalticks=seq(40,50),
                            metrics=c("count sheep", "count wolves", "count patches with [pcolor = green]"),
                            metrics.turtles=list("turtles"= c("energy")),
                            variables = list('initial-number-sheep' = list(min=50, max=150, qfun="qunif"),
                                             'initial-number-wolves' = list(min=50, max=150, qfun="qunif")),
                            constants = list("model-version" = "\"sheep-wolves-grass\"",
                                             "grass-regrowth-time" = 30,
                                             "sheep-gain-from-food" = 4,
                                             "wolf-gain-from-food" = 20,
                                             "sheep-reproduce" = 4,
                                             "wolf-reproduce" = 5,
                                             "show-energy?" = "false"))

nl@simdesign <- simdesign_lhs(nl=nl,
                              samples=10,
                              nseeds=3,
                              precision=3)

# Run in parallel on local machine:
plan(multisession)
results <- run_nl_all(nl)
