## simulation scenarios

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

