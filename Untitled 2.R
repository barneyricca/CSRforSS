# Entropy and sequence length:
entr <- function(le, nu) {
  sample(1:nu, le, replace = TRUE) -> x
  table(x) -> freqs
  freqs / sum(unname(freqs)) -> probs
  -sum(freqs * probs * log(probs,2)) -> en
  return(en)
}

set.seed(42)
1e4 -> N
100 -> num_lens
10 -> spacing
20 -> size

matrix(0,
       ncol = 2,
       nrow = num_lens) -> ent_mat
for(i in 1:num_lens) {
  spacing * i -> len
  replicate(N, entr(len, size)) -> ent_dist
  mean(ent_dist) -> ent_mat[i,1]
  sd(ent_dist) -> ent_mat[i,2]
}

seq(from = spacing, to = spacing * num_lens, by = spacing) -> x
plot(x, ent_mat[,1], pch = 16, cex = 0.6)
plot(x, ent_mat[,2], pch = 16, cex = 0.6)
hist(ent_mat[,1], breaks = 50)

lm(ent_mat[,1] ~ x) -> ent_fit
summary(ent_fit)   # Slope = .3322

# OK, figure out where .3322 comes from in the size = 10 case
# Well:
log(10,2) / 10
# [1] 0.3321928
log(20,2) / 20
## [1] 0.2160964 vs. 0.2161 slope

# Slope is the entropy per sample, so we can just normalize 

# Just because...
library(poweRlaw)
replicate(N, entr(100, 20)) -> ent100
conlnorm$new(ent100) -> ent_ln
ent_ln$setPars(estimate_pars(ent_ln))

plot(ent_ln, pch = 16, cex = 0.6)
lines(ent_ln, col = 2)
# Well, it is a fatter tail than a log-normal, but the log-lot plot
#  doesn't bring about a power law. Maybe a stretched-normal? I dunno,
#  and I don't care so far.



# Now, for a non uniform probability

entr8 <- function(le, pr) {
  sample(1:8, le, replace = TRUE, prob = pr) -> x
  table(x) -> freqs
  freqs / sum(unname(freqs)) -> probs
  -sum(freqs * probs * log(probs,2)) -> en
  return(en)
}
  
set.seed(42)
1e4 -> N
10 -> num_lens
10 -> spacing

c(0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.69) -> prob_dist
  
matrix(0,
       ncol = 2,
       nrow = num_lens) -> ent_mat
for(i in 1:num_lens) {
  spacing * i -> len
  replicate(N, entr8(len, prob_dist)) -> ent_dist
  mean(ent_dist) -> ent_mat[i,1]
  sd(ent_dist) -> ent_mat[i,2]
}
  
seq(from = spacing, to = spacing * num_lens, by = spacing) -> x
plot(x, ent_mat[,1], pch = 16, cex = 0.6)
plot(x, ent_mat[,2], pch = 16, cex = 0.6)
hist(ent_mat[,1], breaks = 50)
  
lm(ent_mat[,1] ~ x) -> ent_fit
summary(ent_fit)   # Slope = .3322


# c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.55) -> prob_dist
# 
# slope .3560349
# log(8,2) / 8 = 0.375
#

# c(0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.69) -> prob_dist
#
# slope 0.2643826