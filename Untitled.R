library(moments)

set.seed(42)

hotstreak <- function(N, num = 1, p = 0) {
  matrix(0, num, N) -> flips
  if(abs(p) > 0.5) {
    cat("Illegal probability adjustment.\nMust be in [-0.5,0.5]")
    return(NULL)
  }
  for(j in 1:num) {
    sample(c(-1,1), 1, prob = c(0.5, 0.5)) -> flips[j, 1]
    for(i in 2:N) {
      sample(c(-1,1), 1, 
            prob = c(0.5 - flips[j, i-1] * p, 
                      0.5 + flips[j, i-1] * p)) -> flips[j, i]
    }
  }  
  return(flips)
}

1e4 -> N      # Number of flips in each replication
1e3 -> n      # Number of replications
0.2 -> pdev
hotstreak(N, n, pdev) -> results
rowSums(results) -> dist1
hist(dist1, breaks = min(n, 30))
mean(dist1)
sd(dist1)
skewness(dist1)




set.seed(20191107)
rcauchy(50, 1, 2) -> a
rcauchy(30, 1.1, 2) -> b
data.frame("Value" = c(a,b),
           "Group" = c(rep("A",50), 
                       rep("B",30))) -> df
library(lmPerm)
aovp(Value ~ Group, data = df) -> df_ap
summary(df_ap)
aov(Value ~ Group, data = df) -> df_a
summary(df_a)
# Is there a problem with resampling from a cauchy?
rcauchy(5000, 2, 59) -> a1
hist(a1, breaks = 100)
min(a1)
max(a1)

# Is a cauchy the ratio of two normals?
rnorm(5000) -> an
rnorm(5000) -> bn
an/bn -> cn
head(sort(cn),20)
fitdistr(cn, "cauchy")
rcauchy(5000, 0.004, 1.015) -> dn
head(sort(dn), 20)


# Comparison of plots
library(mlbench)
data(PimaIndiansDiabetes)

# Get a shorter name

pm <- PimaIndiansDiabetes

par(mfrow=c(2,4))

for (ii in 1:(nrow(pm)-1)) {
  form <- as.formula(paste(names(pm)[ii]," ~ diabetes",sep=""))
  boxplot(form,data=pm,main=names(pm)[ii])
  grid()
}
par(mfrow=c(1,1)) # Reset the plot window

# vs
gather(pm,key="variable",val="value",-diabetes) %>%
  ggplot(aes(x=diabetes,y=value)) + 
  geom_boxplot() + 
  facet_wrap(~variable,scales="free") + 
  theme_bw()

#
splitdf <- split(pm,pm$diabetes)
par(mfrow=c(2,4))

for (ii in 1:8) {
  a <- paste("splitdf[[1]]$",names(pm)[ii],sep="")
  plot(density(eval(parse(text=a))),main=names(pm)[ii],col="red")
  a <- paste("splitdf[[2]]$",names(pm)[ii],sep="")
  lines(density(eval(parse(text=a))),col="blue")
  legend("bottomright",
         legend=c("neg","pos"),
         col=c("red","blue"),lty=1:1,cex=0.8)
  grid()
  
}
par(mfrow=c(1,1))


install.packages("DescTools")
library(DescTools)

set.seed(42)
sample.int(6, size = 1e5, replace = TRUE) -> a
sample.int(6, size = 1e5, replace = TRUE) -> b
MutInf(a, b)   # Should be zero

as.matrix(c(49, 27652, 141, 774106),
          nrow = 2,
          ncol = 2,
          byrow = TRUE) -> mat1
c(2,2) -> dim(mat1)
t(mat1) -> mat1
MutInf(mat1) # 0.0001105, according to https://nlp.stanford.edu/IR-book/html/htmledition/mutual-information-1.html

c(1/8, 1/16, 1/32, 1/32, 1/16, 1/8, 1/32, 1/32,
  1/16, 1/16, 1/16, 1/16, 1/4, 0, 0, 0) -> bt
c(4,4) -> dim(bt)
t(bt) -> bt
MutInf(bt) # 0.375 by https://www2.isye.gatech.edu/~yxie77/ece587/Lecture2.pdf


sample.int(6, 1e6, replace = TRUE) -> am
sample.int(6, 1e6, replace = TRUE) -> bm
c(1e3,1e3) -> dim(am) -> dim(bm)
MutInf(am, bm) # Should be zero



# Making a transition matrix / network
library(foreach)

c(13, 13,  8, 6, 13, 8, 6, 8, 12, 14, 3, 2, 1, 3, 3, 13,
  13, 8, 3, 9, 4, 12, 12, 8, 3, 3, 13, 13, 13, 13, 9, 8,
  12, 3, 9, 8, 9, 13, 4, 3, 8, 12, 4, 3, 9) -> samp_seq
unname(table(samp_seq)) -> symb_freqs
sort(unique(samp_seq)) -> syms
length(syms) -> num_syms

foreach(i=1:length(samp_seq), .combine = c) %do% {
  which(syms == samp_seq[i])} -> columns

matrix(0, 
       nrow = num_syms, 
       ncol = num_syms,
       dimnames = list(syms, syms)) -> trans_mat

for(i in 1:(length(samp_seq) - 1)) {
  trans_mat[columns[i], columns[i+1]] + 1 -> 
    trans_mat[columns[i], columns[i+1]]
}
trans_mat

library(network)
network(trans_mat,
        vertex.attrnames = syms,
        loops = TRUE) -> trans_net
trans_mat -> edge_mat
network::set.edge.attribute(trans_net,
                            "Weight",
                            edge_mat)
network::set.vertex.attribute(trans_net,
                              "Frequency",
                              as.character(symb_freqs))
trans_net

set.seed(42)
plot.network(trans_net,
     displaylabels = TRUE,
     edge.curve = TRUE,
     vertex.cex = 25 * sqrt(network::get.vertex.attribute(trans_net,
                                                         "Frequency")),
     edge.lwd = 2 * network::get.edge.attribute(trans_net,
                                                "Weight"))

# End of making the transition matrix / network




# ```{r clusters}
c(
  "magrittr",
  "cluster",
  "cluster.datasets",
  "cowplot",
  "NbClust",
  "clValid",
  "ggfortify",
  "clustree",
  "dendextend",
  "factoextra",
  "FactoMineR",
  "corrplot",
  "GGally",
  "ggiraphExtra",
) -> package_names  

data("all.mammals.milk.1956")
#raw_mammals <- all.mammals.milk.1956
# subset dataset
all.mammals.milk.1956 %>%    
  select(., -name) %>%       # subset data
  as_tibble(.) -> mammals    # Make it a tibble
# mammals <- as_tibble(mammals)

summary(mammals) %>% kable() %>% kable_styling()

mammals %>% 
  gather(Attributes, value, 1:5) %>% 
  ggplot(aes(x=value)) +
  geom_histogram(fill = "lightblue2", color = "black") + 
  facet_wrap(~Attributes, scales = "free_x") +
  labs(x = "Value", y = "Frequency")

corrplot(cor(mammals), type = "upper", method = "ellipse", tl.cex = 0.9)

mammals_scaled <- scale(mammals)
rownames(mammals_scaled) <- all.mammals.milk.1956$name

res.pca <- PCA(mammals_scaled,  graph = FALSE)
# Visualize eigenvalues/variances
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))

# Extract the results for variables
var <- get_pca_var(res.pca)
# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
# Control variable colors using their contributions to the principle axis
fviz_pca_var(res.pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
) + theme_minimal() + ggtitle("Variables - PCA")

kmean_calc <- function(df, ...){
  kmeans(df, scaled = ..., nstart = 30)
}
km2 <- kmean_calc(mammals_scaled, 2)
km3 <- kmean_calc(mammals_scaled, 3)
km4 <- kmeans(mammals_scaled, 4)
km5 <- kmeans(mammals_scaled, 5)
km6 <- kmeans(mammals_scaled, 6)
km7 <- kmeans(mammals_scaled, 7)
km8 <- kmeans(mammals_scaled, 8)
km9 <- kmeans(mammals_scaled, 9)
km10 <- kmeans(mammals_scaled, 10)
km11 <- kmeans(mammals_scaled, 11)
p1 <- fviz_cluster(km2, data = mammals_scaled, frame.type = "convex") + theme_minimal() + ggtitle("k = 2") 
p2 <- fviz_cluster(km3, data = mammals_scaled, frame.type = "convex") + theme_minimal() + ggtitle("k = 3")
p3 <- fviz_cluster(km4, data = mammals_scaled, frame.type = "convex") + theme_minimal() + ggtitle("k = 4")
p4 <- fviz_cluster(km5, data = mammals_scaled, frame.type = "convex") + theme_minimal() + ggtitle("k = 5")
p5 <- fviz_cluster(km6, data = mammals_scaled, frame.type = "convex") + theme_minimal() + ggtitle("k = 6")
p6 <- fviz_cluster(km7, data = mammals_scaled, frame.type = "convex") + theme_minimal() + ggtitle("k = 7")
plot_grid(p1, p2, p3, p4, p5, p6, labels = c("k2", "k3", "k4", "k5", "k6", "k7"))

gap_stat <- clusGap(mammals_scaled, FUN = kmeans, nstart = 30, K.max = 24, B = 50)
fviz_gap_stat(gap_stat) + theme_minimal() + ggtitle("fviz_gap_stat: Gap Statistic")



# Remember: CFA is a subset of SEM

# Simple CFA Model
mydata.cov <- cov(mydata)
model.mydata <- specify.model()
F1 ->  X1, lam1, NA
F1 ->  X2, lam2, NA
F1 ->  X3, lam3, NA
F2 ->  X4, lam4, NA
F2 ->  X5, lam5, NA
F2 ->  X6, lam6, NA
X1 <-> X1, e1,   NA
X2 <-> X2, e2,   NA
X3 <-> X3, e3,   NA
X4 <-> X4, e4,   NA
X5 <-> X5, e5,   NA
X6 <-> X6, e6,   NA
F1 <-> F1, NA,    1
F2 <-> F2, NA,    1
F1 <-> F2, F1F2, NA

mydata.sem <- sem(model.mydata, mydata.cov, nrow(mydata))
# print results (fit indices, paramters, hypothesis tests)
summary(mydata.sem)
# print standardized coefficients (loadings)
std.coef(mydata.sem)




library(fractal)
library(fractaldim)
library(tseriesChaos)

ts(logistic_series) -> log_ts
embedd(log_ts, m = 2, d = 1) -> log_e
plot(log_e, cex = 0.3)

x <- window(rossler.ts, start=90)
xyz <- embedd(x, m=3, d=8)


library(devtools)
devtools::install_github("ropensci/workloopR")

library(workloopR)

## import the workloop.ddf file included within workloopR
wl_dat <-read_ddf(system.file("extdata", "workloop.ddf", 
                              package = 'workloopR'),
                  phase_from_peak = TRUE)
Important metadata from the file are stored as attributes of the object. Here’s what’s stored for the work loop file:
  
  names(attributes(wl_dat))
