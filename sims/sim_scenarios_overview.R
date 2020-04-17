library(scp)
library(fields)
library(viridis)

source("sims/sim_scenarios.R")

# setup
theta        = c(0,3,0.1,0.7)
names(theta) = c("Nugget","PartialSill","Range","Smoothness")

N = 40
S = seq(0,1,length=N)
n = N^2
s = expand.grid(S,S)
d = as.matrix(dist(s))
C = mat_cov(d,theta)
G = t(chol(C))%*%rnorm(n)
E = rnorm(n)

# simulation data
Y1 = stat_sim(G,E,s)
Y2 = cubic_sim(G,E,s)
Y3 = skew_sim(G,E,s)
Y4 = prod_sim(G,E,s)
Y5 = gfun_sim(G,E,s)
Y6 = eastwest_sim(G,E,s)
Y7 = nugget_sim(G,E,s)
Y8 = spike_sim(G,E,s)

# visualize one realization
pdf("man/figures/SIM-scenarios-overview.pdf", height = 10, width = 6)
par(mfrow = c(4,2))

plot_simY = function(Y,main){
  image.plot(S,S,matrix(Y,N,N),col = viridis(50, option = "D"),
             las=1, xlab = expression(s[x]), ylab = expression(s[y]),
             main = main)}

plot_simY(Y1, main = expression("Scenario 1. Y = G + E"))
plot_simY(Y2, main = expression(paste("Scenario 2. Y = ", G^3, " + E")))
plot_simY(Y3, main = expression(paste("Scenario 3. Y = q(", Phi, "(G/", sqrt(3),")) + E")))
plot_simY(Y4, main = expression(paste("Scenario 4. Y = ", sqrt(3), G%.% "|", E, "|")))
plot_simY(Y5, main = expression(paste("Scenario 5. Y = ", sign(G)%.%abs(G)^{s[x]+1}, " + E")))
plot_simY(Y6, main = expression(paste("Scenario 6. Y = ", sqrt(w)%.%G/sqrt(3) + sqrt(1-w)%.%E)))
plot_simY(Y7, main = expression(paste("Scenario 7. Y = G + Normal(0,",s[x],")")))
plot_simY(Y8, main = expression("Scenario 8. Y = G + spike"))

dev.off()

