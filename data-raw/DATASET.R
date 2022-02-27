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
w = pnorm(s[,1],0.5,0.1)
Y = sqrt(w)*X/sqrt(3) + sqrt(1-w)*rnorm(n)

names(s) = c("dim1", "dim2")
sample_data = list(s=s, Y=Y)
usethis::use_data(sample_data, compress = "xz", overwrite = TRUE)



# usethis::use_data("DATASET")
