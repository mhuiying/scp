library(abind)
library(fields)
library(geoR)
library(viridis)
library(doParallel)
cl = makeCluster(4)
registerDoParallel(cl)

# Simulation parameters
eps     = 0.1          # Significance level
sim_cnt = 4            # No. of methods: global conformal, local conformal, Kriging, laGP
N_sim   = 100          # repetition for each simulation

# True values
theta        = c(0,3,0.1,0.7)
names(theta) = c("Nugget","PartialSill","Range","Smoothness")

for(N in c(20,40)){

  print(N)

  # Fake data
  S = seq(0,1,length=N)
  n = N^2
  s = expand.grid(S,S)
  d = as.matrix(dist(s))
  C = mat_cov(d,theta)

  # 3D array to store result
  cover    = array(NA, dim=c(sim_cnt, N_sim, n))
  width    = array(NA, dim=c(sim_cnt, N_sim, n))
  intScore = array(NA, dim=c(sim_cnt, N_sim, n))

  for(sim in 1:N_sim){
    flag = TRUE
    while(flag == TRUE){
      X0 = t(chol(C))%*%rnorm(n)
      X = X0^3
      Y = X + rnorm(n)

      # Estimate spatial covariance parameters
      bins     = seq(0.01,0.2,0.01)
      thetaHat = get_theta(s,Y,dists=bins,plot_fitted=FALSE)

      # Get prediction intervals at all n sites
      Q        = try( solve(mat_cov(d,thetaHat)) )  # avoid possible singularity, especially the gamma situation
      flag = any(is.nan(Q))
    }

    gamma = foreach(j = 1:n) %dopar% {
      s0 = as.numeric(s[j,])
      rbind(
        scp(s0,s[-j,],Y[-j],alpha=eps),         # GSCP
        scp(s0,s[-j,],Y[-j],eta=0.1,alpha=eps), # LSCP
        Krige_int(Q,Y,j,eps=eps),
        laGP_int(s0,s[-j,],Y[-j],eps=eps)
      )
    }
    gamma = do.call(abind, c(gamma, list(along=3)))

    # Summarize the output
    for(situ in 1:sim_cnt){
      width[situ,sim, ]    = gamma[situ,2,]-gamma[situ,1,]
      cover[situ,sim, ]    = (gamma[situ,1,]<=Y) & (Y<=gamma[situ,2,])
      intScore[situ,sim, ] = interval_score(l = gamma[situ,1,], u = gamma[situ,2,], Y)
    }
  }

  save(width, cover, intScore, file = paste("data/sim_results/cubic_N_", N, ".RData", sep = ""))

}
stopCluster(cl)
