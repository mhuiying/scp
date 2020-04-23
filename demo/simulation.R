library(scp)

## Simulation parameters
alpha   = 0.1          # Significance level
sim_cnt = 4            # No. of methods to compare: GSCP, LSCP, Kriging, laGP
N_sim   = 100          # repetition for each simulation
result_loc = "data/sim_results/"

## Simulation scenarios based on
## spatial observations G and white noise E
# scenario 1: stationary process
stat_sim = function(G,E,s){
  Y = G + E
  return(Y)
}
# scenario 2: cubic
cubic_sim = function(G,E,s){
  Y = G^3 + E
  return(Y)
}
# scenario 3: skewed distribution
skew_sim = function(G,E,s){
  X = qgamma(pnorm(G,sd=sqrt(3)),1,1/sqrt(3))
  Y = X + E
  return(Y)
}
# scenario 4: product of G and E
prod_sim = function(G,E,s){
  X = sqrt(3)*G
  Y = X * abs(E)
  return(Y)
}
# scenario 5: use g function to transform
gfun_sim = function(G,E,s){
  g = function(x, p){
    sign(x) * abs(x)^p
  }
  X = g(G, s[,1]+1)
  Y = X + E
  return(Y)
}
# scenario 6: process varies from west to east
eastwest_sim = function(G,E,s){
  w = pnorm(s[,1],0.5,0.1)
  Y = sqrt(w)*G/sqrt(3) + sqrt(1-w)*E
  return(Y)
}
# scenario 7: nugget change
nugget_sim = function(G,E,s){
  Y = G + s[,1]*E
  return(Y)
}
# scenario 8: adding a spike in the center
spike_sim = function(G,E,s){
  dist = sqrt( (0.5-s[,1])^2 + (0.5-s[,2])^2 )
  Y    = G + dnorm(dist,sd=0.1)/dnorm(0,sd=0.1)*10
  return(Y)
}

# Matern parameters
theta        = c(0,3,0.1,0.7)
names(theta) = c("Nugget","PartialSill","Range","Smoothness")

## Generate scenarios overview figure
## as Figure 1 in Mao. et al (2020)
# N = 40
# S = seq(0,1,length=N)
# n = N^2
# s = expand.grid(S,S)
# d = as.matrix(dist(s))
# C = mat_cov(d,theta)
# G = t(chol(C))%*%rnorm(n)
# E = rnorm(n)
#
# # simulation data
# Y1 = stat_sim(G,E,s)
# Y2 = cubic_sim(G,E,s)
# Y3 = skew_sim(G,E,s)
# Y4 = prod_sim(G,E,s)
# Y5 = gfun_sim(G,E,s)
# Y6 = eastwest_sim(G,E,s)
# Y7 = nugget_sim(G,E,s)
# Y8 = spike_sim(G,E,s)
#
# # visualize one realization
# pdf("man/figures/SIM-scenarios-overview.pdf", height = 10, width = 6)
# par(mfrow = c(4,2))
#
# plot_simY = function(Y,main){
#   image.plot(S,S,matrix(Y,N,N),col = viridis(50, option = "D"),
#              las=1, xlab = expression(s[x]), ylab = expression(s[y]),
#              main = main)
# }
#
# plot_simY(Y1, main = expression("Scenario 1. Y = G + E"))
# plot_simY(Y2, main = expression(paste("Scenario 2. Y = ", G^3, " + E")))
# plot_simY(Y3, main = expression(paste("Scenario 3. Y = q(", Phi, "(G/", sqrt(3),")) + E")))
# plot_simY(Y4, main = expression(paste("Scenario 4. Y = ", sqrt(3), G%.% "|", E, "|")))
# plot_simY(Y5, main = expression(paste("Scenario 5. Y = ", sign(G)%.%abs(G)^{s[x]+1}, " + E")))
# plot_simY(Y6, main = expression(paste("Scenario 6. Y = ", sqrt(w)%.%G/sqrt(3) + sqrt(1-w)%.%E)))
# plot_simY(Y7, main = expression(paste("Scenario 7. Y = G + Normal(0,",s[x],")")))
# plot_simY(Y8, main = expression("Scenario 8. Y = G + spike"))
#
# dev.off()

## Simulation starts here
for(sim_scenario in c("stat", "cubic", "skew", "prod", "gfun", "eastwest", "nugget", "spike")){

  if(sim_scenario == "stat")
    sim_fun = stat_sim
  if(sim_scenario == "cubic")
    sim_fun = cubic_sim
  if(sim_scenario == "skew")
    sim_fun = skew_sim
  if(sim_scenario == "prod")
    sim_fun = prod_sim
  if(sim_scenario == "gfun")
    sim_fun = gfun_sim
  if(sim_scenario == "eastwest")
    sim_fun = eastwest_sim
  if(sim_scenario == "nugget")
    sim_fun = nugget_sim
  if(sim_scenario == "spike")
    sim_fun = spike_sim

  if(sim_scenario == "spike"){
    N_list = c(20,40,60)
  } else{
    N_list = c(20,40)
  }

  for(N in N_list){

    # simulate on grid locations
    n = N^2
    S = seq(0,1,length=N)
    s = expand.grid(S,S)
    # # for simulation on uniformly sampled locations
    # # replace the above two lines with the following three commented code
    # Var1 = runif(n)
    # Var2 = runif(n)
    # s    = data.frame(Var1,Var2)

    d = as.matrix(dist(s))
    C = mat_cov(d,theta)

    # 3D array to store result
    cover    = array(NA, dim=c(sim_cnt, N_sim, n))
    width    = array(NA, dim=c(sim_cnt, N_sim, n))
    intScore = array(NA, dim=c(sim_cnt, N_sim, n))

    for(sim in 1:N_sim){
      flag = TRUE
      while( flag ){
        G = t(chol(C))%*%rnorm(n)
        E = rnorm(n)
        Y = sim_fun(G,E,s)

        # Estimate spatial covariance parameters
        bins     = seq(0.01,0.2,0.01)
        thetaHat = get_theta(s,Y,dists=bins,plot_fitted=FALSE)

        # Get prediction intervals at all n sites
        Q        = try( solve(mat_cov(d,thetaHat)) )  # avoid possible singularity, especially the gamma situation
        flag = any(is.nan(Q))
      }

      gamma = array(NA, dim=c(sim_cnt, 2, n))
      for(i in 1:n){
        s0 = as.numeric(s[i,])
        gamma[1,,i] = scp(s0=s0,s=s[-i,],Y=Y[-i],alpha=alpha)         # GSCP
        gamma[2,,i] = scp(s0=s0,s=s[-i,],Y=Y[-i],eta=0.1,alpha=alpha) # LSCP
        gamma[3,,i] = Krige_int(Q,Y,i,alpha=alpha)                    # Kriging
        gamma[4,,i] = laGP_int(s0,s[-i,],Y[-i],alpha=alpha)           # laGP
      }

      # Summarize the output
      for(situ in 1:sim_cnt){
        width[situ,sim, ]    = gamma[situ,2,]-gamma[situ,1,]
        cover[situ,sim, ]    = (gamma[situ,1,]<=Y) & (Y<=gamma[situ,2,])
        intScore[situ,sim, ] = interval_score(l = gamma[situ,1,], u = gamma[situ,2,], Y)
      }
    }

    save(s, width, cover, intScore,
         file = paste(result_loc, sim_scenario, "_N_", N, ".RData", sep = ""))

  }
}

## Analyze and present simulation results
# For scenarios stat, cubic, skew, and prod
# results are presented in tabular format
# as shown in Table 1 in Mao. et al (2020)

assemble_result_tabular = function(sim_scenario, N){

  result_file = paste(result_loc, sim_scenario, "_N_", N, ".RData", sep="")
  load(result_file)

  AveWidth    = apply(width, MARGIN = 1, FUN = mean, na.rm = TRUE)
  AveCover    = apply(cover, MARGIN = 1, FUN = mean, na.rm = TRUE)
  AveIntScore = apply(intScore, MARGIN = 1, FUN = mean, na.rm = TRUE)

  result          = as.data.frame(cbind(AveCover, AveWidth, AveIntScore))
  result$AveCover = paste(round(result$AveCover*100, digits = 1), "%", sep = "")
  return(result)
}

for(sim_scenario in c("stat", "cubic", "skew", "prod")){

  print(sim_scenario)
  result1 = assemble_result_tabular(sim_scenario = sim_scenario, N = 20)
  result2 = assemble_result_tabular(sim_scenario = sim_scenario, N = 40)
  print(cbind(result1, result2))

}

# For scenarios gfun, eastwest, and nugget
# results are presented by aggregating over sy and over datasets
# as shown in Figure 4, 8, and 9 in Mao. et al (2020)

assemble_result_sy = function(sim_scenario, N){

  result_file = paste(result_loc, sim_scenario, "_N_", N, ".RData", sep="")
  load(result_file)

  aggAveWidth = aggAveCover = aggIntScore = matrix(NA, N, sim_cnt)
  for(situ in 1:sim_cnt){
    AveWidth = apply(width[situ, , ], 2, mean, na.rm = TRUE)
    AveWidth = as.data.frame(cbind(AveWidth, s))
    aggAveWidth[,situ] = aggregate(AveWidth ~ Var1, data = AveWidth, mean)[,2]
    AveCover = apply(cover[situ, , ], 2, mean, na.rm = TRUE)
    AveCover = as.data.frame(cbind(AveCover, s))
    aggAveCover[,situ] = aggregate(AveCover ~ Var1, data = AveCover, mean)[,2]
    AveIntScore = apply(intScore[situ, , ], 2, mean, na.rm = TRUE)
    AveIntScore = as.data.frame(cbind(AveIntScore, s))
    aggIntScore[,situ] = aggregate(AveIntScore ~ Var1, data = AveIntScore, mean)[,2]
  }

  result = list(x = seq(0,1,length=N),
                aggAveWidth = aggAveWidth,
                aggAveCover = aggAveCover,
                aggIntScore = aggIntScore)
  return(result)

}

plot_result = function(x, y_metric, ylim, xlab = expression(s[x]), ylab, main = "", legend_loc = "top"){

  if(dim(y_metric)[1] != length(x))
    stop("Dimension not match")

  matplot(x, y_metric, type = "b", lwd = 1.5, ylim = ylim,
          pch = c(1,2,3,4), lty = c(1,2,3,4), col = c(1,2,3,4),
          xlab = xlab, ylab = ylab, main = main)
  legend(legend_loc, c("GSCP","LSCP", "Kriging", "laGP"),
         lty = c(1,2,3,4), pch = c(1,2,3,4), col = c(1,2,3,4))

}

for(sim_scenario in c("gfun", "eastwest", "nugget")){

  # N = 20 result
  result1 = assemble_result_sy(sim_scenario, N=20)
  x1           = result1$x
  aggAveWidth1 = result1$aggAveWidth
  aggAveCover1 = result1$aggAveCover
  aggIntScore1 = result1$aggIntScore

  # N = 40 result
  result2 = assemble_result_sy(sim_scenario, N=40)
  x2           = result2$x
  aggAveWidth2 = result2$aggAveWidth
  aggAveCover2 = result2$aggAveCover
  aggIntScore2 = result2$aggIntScore

  # plot
  par(mfrow = c(3,2))
  width_ylim = range(c(aggAveWidth1, aggAveWidth2))
  plot_result(x1, aggAveWidth1, ylim=width_ylim, ylab="Ave Width", main="N=20")
  plot_result(x2, aggAveWidth2, ylim=width_ylim, ylab="Ave Width", main="N=40")

  cover_ylim = range(c(aggAveCover1, aggAveCover2))
  plot_result(x1, aggAveCover1, ylim=cover_ylim, ylab="Ave Coverage")
  abline(h = 1-alpha, lty = 2)
  plot_result(x2, aggAveCover2, ylim=cover_ylim, ylab="Ave Coverage")
  abline(h = 1-alpha, lty = 2)

  score_ylim = range(c(aggIntScore1, aggIntScore2))
  plot_result(x1, aggIntScore1, ylim=score_ylim, ylab="Interval Score")
  plot_result(x2, aggIntScore2, ylim=score_ylim, ylab="Interval Score")

}

# For scenario spike, simulations have conducted for N = 20, 40, 60
# results are presented by aggregating over dist2center and over datasets
# as shown in Figure 10 in Mao. et al (2020)

assemble_dist2center_result = function(N, bin_width = 0.05, center = c(0.5, 0.5)){

  result_file = paste(result_loc, "spike_N_", N, ".RData", sep="")
  load(result_file)

  dist2center = sqrt( (center[1]-s[,1])^2 + (center[2]-s[,2])^2 )
  bins        = cut(dist2center, breaks = seq(from = 0, to = 1, by=bin_width))
  uniq_bins   = unique(bins)

  # sort uniq_bins
  labs      = as.vector(uniq_bins)
  upper     = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) )
  uniq_bins = uniq_bins[order(upper)]

  n_bins      = length(uniq_bins)
  aggAveWidth = aggAveCover = aggIntScore = matrix(NA, n_bins, sim_cnt)
  for(situ in 1:sim_cnt){
    for(i in 1:n_bins){
      dist = uniq_bins[i]
      aggAveWidth[i,situ] = mean(width[situ, , bins == dist])
      aggAveCover[i,situ] = mean(cover[situ, , bins == dist])
      aggIntScore[i,situ] = mean(intScore[situ, , bins == dist])
    }
  }

  result = list(x = sort(upper)-bin_width/2,
                aggAveWidth = aggAveWidth,
                aggAveCover = aggAveCover,
                aggIntScore = aggIntScore)
  return(result)

}

# N = 20 result
result1      = assemble_dist2center_result(N=20)
x1           = result1$x
aggAveWidth1 = result1$aggAveWidth
aggAveCover1 = result1$aggAveCover
aggIntScore1 = result1$aggIntScore

# N = 40 result
result2      = assemble_dist2center_result(N=40)
x2           = result2$x
aggAveWidth2 = result2$aggAveWidth
aggAveCover2 = result2$aggAveCover
aggIntScore2 = result2$aggIntScore

# N = 60 result
result3      = assemble_dist2center_result(N=60)
x3           = result3$x
aggAveWidth3 = result3$aggAveWidth
aggAveCover3 = result3$aggAveCover
aggIntScore3 = result3$aggIntScore

# plot
xlab = "distance to center"
par(mfrow = c(3,3))

width_ylim = range(c(aggAveWidth1, aggAveWidth2, aggAveWidth3))
plot_result(x1, aggAveWidth1, ylim=width_ylim, xlab = xlab, ylab="Ave Width", main="N=20")
plot_result(x2, aggAveWidth2, ylim=width_ylim, xlab = xlab, ylab="Ave Width", main="N=40")
plot_result(x3, aggAveWidth3, ylim=width_ylim, xlab = xlab, ylab="Ave Width", main="N=60")

legend_loc = "bottom"
cover_ylim = range(c(aggAveCover1, aggAveCover2, aggAveCover3))
plot_result(x1, aggAveCover1, ylim=cover_ylim, xlab = xlab, ylab="Ave Coverage", legend_loc = legend_loc)
abline(h = 1-alpha, lty = 2)
plot_result(x2, aggAveCover2, ylim=cover_ylim, xlab = xlab, ylab="Ave Coverage", legend_loc = legend_loc)
abline(h = 1-alpha, lty = 2)
plot_result(x3, aggAveCover3, ylim=cover_ylim, xlab = xlab, ylab="Ave Coverage", legend_loc = legend_loc)
abline(h = 1-alpha, lty = 2)

score_ylim = range(c(aggIntScore1, aggIntScore2, aggIntScore3))
plot_result(x1, aggIntScore1, ylim=score_ylim, xlab = xlab, ylab="Interval Score")
plot_result(x2, aggIntScore2, ylim=score_ylim, xlab = xlab, ylab="Interval Score")
plot_result(x3, aggIntScore3, ylim=score_ylim, xlab = xlab, ylab="Interval Score")
