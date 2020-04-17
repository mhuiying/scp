## code to prepare `DATASET` dataset goes here

set.seed(123)
N = 41; n = N^2
S = seq(0,1,length=N)
s = expand.grid(S,S)
d = as.matrix(dist(s))

theta        = c(0,3,0.1,0.7)
names(theta) = c("Nugget","PartialSill","Range","Smoothness")
C = mat_cov(d,theta)
X = t(chol(C))%*%rnorm(n)
Y = X^3 + rnorm(n)

sample_data = list(s=s, Y=Y)
usethis::use_data(sample_data, compress = "xz")



# usethis::use_data("DATASET")
