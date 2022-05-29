MOVE <- function(){
  #trunc[1, 1] <- lat_min
  #trunc[1, 2] <- lat_max
  #trunc[2, 1] <- lon_min
  #trunc[2, 2] <- lon_max
  
  sd ~ dunif(0, 0.001)
  tau <- 1/(sd*sd)
  
  phi ~ dunif(0, 1)
  
  # variance-covariance matrix
  sigma_c[1, 1, 1] ~ dexp(1) # it cannot be negative
  sigma_c[1, 2, 1] <- 0
  sigma_c[2, 1, 1] <- 0
  sigma_c[2, 2, 1] ~ dexp(1)
  sigma_c[1, 1, 2] ~ dexp(1)
  sigma_c[1, 2, 2] <- 0
  sigma_c[2, 1, 2] <- 0
  sigma_c[2, 2, 2] ~ dexp(1)
  
  Rho[1, 1, 1] <- 1.0
  Rho[1, 2, 1] ~ dunif(-0.8, 0.8)
  Rho[2, 1, 1] <- Rho[1, 2, 1]
  Rho[2, 2, 1] <- 1.0
  Rho[1, 1, 2] <- 1.0
  Rho[1, 2, 2] ~ dunif(-0.8, 0.8)
  Rho[2, 1, 2] <- Rho[1, 2, 2]
  Rho[2, 2, 2] <- 1.0
  
  Sigma[1:2, 1:2, 1] <- sigma_c[,,1] %*% Rho[,,1] %*% sigma_c[,,1]
  Sigma[1:2, 1:2, 2] <- sigma_c[,,2] %*% Rho[,,2] %*% sigma_c[,,2]
  
  # priors on strength of autoregressive component for two states
  gamma[1] ~ dunif(0, 1)
  #u ~ dunif(0, 1)
  #gamma[2] <- gamma[1]*u
  gamma[2] ~ dunif(0, 1)
  
  # priors for drift component (one for x and one for y dim.)
  D[1, 1] ~ dnorm(0, 0.1)
  D[1, 2] ~ dnorm(0, 0.1)
  D[2, 1] ~ dnorm(0, 0.1)
  D[2, 2] ~ dnorm(0, 0.1)
  
  # N is the individual index
  for(k in 1:N){
    first.loc[k, 1] <- y[Yidx[k], 1]
    first.loc[k, 2] <- y[Yidx[k], 2]
    # mute measurement error components
    #logpsi[k] ~ dunif(-10, 10)
    #psi[k] <- exp(logpsi[k])
    b2[Xidx[k]] <- 1
    b[Xidx[k]] <- 2
    for(j in 1:2){
      # mute the measurement error component
      #x[Xidx[k], j] ~ dt(first.loc[k, j], itau2[Xidx[k], j]*psi[k], nu[Xidx[k], j])
      x[Xidx[k], j] <- first.loc[k, j]
    }
    # get x for the second recorded location, which does not use autoregressive component
    x[(Xidx[k]+1), 1:2] ~ dmnorm.vcov(x[Xidx[k],], Sigma[, , b[Xidx[k]]])
    
    # mute loglikelihood monitor for first two locations
    #log_lik[(Xidx[k])] <- logdensity.mnorm(x[(Xidx[k]), 1:2], x[(Xidx[k]), ], iSigma[, ]) 
    #log_lik[(Xidx[k]+1)] <- logdensity.mnorm(x[(Xidx[k]+1), 1:2], x[(Xidx[k]+1), ], iSigma[, ])
    # for regular time steps, t
    for(t in (Xidx[k]+1):(Xidx2[k]-1)){ 
      phi2[t] <- max(phi*b2[t-1], 0.001)
      b2[t] ~ dbern(phi2[t])
      b[t] <- b2[t] + 1
      # displacement at time t + 1, from t, includes an autoregressive component, drift, and a max hop distance
      # for x and y dims. 
      # combine x and y dims. and add displacement to the previous time step, t, to feed into dmnorm()
      displace[t, 1] <- max(-10, min((x[t, 1] - x[t-1, 1])*gamma[b[t]] + D[1, b[t]], 10))
      displace[t, 2] <- max(-10, min((x[t, 2] - x[t-1, 2])*gamma[b[t]] + D[2, b[t]], 10))
      # get x for t + 1
      x.mn[t, 1:2] <- x[t, 1:2] + displace[t, 1:2]
      #(x[t, 1:2] - x[t-1, 1:2])*gamma + D[1:2]
      x[t+1, 1:2] ~ dmnorm.vcov(x.mn[t, ], Sigma[, , b[t]])
    }

    for(i in (Yidx[k]+1):(Yidx[k+1]-1)) { 
      for(j in 1:2) {
        # mute the component that interpolates multiple detection within a single time step
        #yhat[i, j] <- w[i]*x[idx[i], j] + (1 - w[i])*x[idx[i+1], j] 
        yhat[i, j] <- x[idx[i], j]
        # mute the measurement error component
        #y[i, j] ~ dt(yhat[i, j], itau2[i, j]*psi[k], nu[i, j])
        y[i, j] ~ dnorm(yhat[i, j], tau)
        #;T(trunc[j, 1], trunc[j, 2])
      }
    }
  }
}

if (is.R()){
  filename <- file.path(tempdir(), "MOVE.bug")}
write.model(MOVE, filename)
inits <- list(list(gamma=c(0.5, 0.5), x=x, D=matrix(c(0,0,0,0), 2, 2), sigma_c=array(c(0.1,NA,NA,0.1,0.1,NA,NA,0.1), c(2, 2, 2)), sd=0.000001, phi=0.1, 
                   Rho=array(c(NA,NA,0.01,NA,NA,NA,0.01,NA), c(2, 2, 2))))
data <- list("Xidx","Xidx2", "Yidx", "y", "idx","N") 
             #"lat_min", "lat_max", "lon_min", "lon_max")
parameters <- c("gamma", "Sigma", "D", "b", "x", "phi")
MOVE <- jags(data=data, inits=inits, parameters.to.save=parameters, filename,
             n.chains=1, n.burnin=n.burnin, n.iter=n.iter+n.burnin, n.thin=1, DIC=TRUE)
MOVE.mcmc <- as.mcmc(MOVE)


b_post <- mat.or.vec(N*season_length, 100)
b_post[b_post==0] <- NA
for(q in 1:N){
  for(i in Xidx[q]:(Xidx2[q]-1)){
    b_post[i, ] <- MOVE$BUGSoutput$sims.array[, , paste("b[", i, "]", sep="")][1:100]
  }
}

# save MCMC output to objects to use for posterior checks
b_post <- mat.or.vec(N*season_length, 1)
b_post[b_post==0] <- NA
for(q in 1:N){
for(i in Xidx[q]:(Xidx2[q]-1)){
  b_post[i] <- mean(MOVE$BUGSoutput$sims.array[, , paste("b[", i, "]", sep="")])
}
}
b_post <- matrix(b_post, nrow=N, ncol=season_length, byrow=TRUE)
gamma_post <- mat.or.vec(2, length(MOVE$BUGSoutput$sims.array[, , "gamma[1]"]))
gamma_post[1,] <- MOVE$BUGSoutput$sims.array[, , "gamma[1]"]
gamma_post[2,] <- MOVE$BUGSoutput$sims.array[, , "gamma[2]"]
D_post <- array(NA, c(2, 2, length(MOVE$BUGSoutput$sims.array[, , "D[1,1]"])))
D_post[1, 1, ] <- MOVE$BUGSoutput$sims.array[, , "D[1,1]"]
D_post[2, 1, ] <- MOVE$BUGSoutput$sims.array[, , "D[2,1]"]
D_post[1, 2, ] <- MOVE$BUGSoutput$sims.array[, , "D[1,2]"]
D_post[2, 2, ] <- MOVE$BUGSoutput$sims.array[, , "D[2,2]"]
sigma_post1 <- array(NA, c(2, 2, length(MOVE$BUGSoutput$sims.array[, , "Sigma[1,1,1]"])))
sigma_post1[1,1, ] <- MOVE$BUGSoutput$sims.array[, , "Sigma[1,1,1]"]
sigma_post1[2,1, ] <- MOVE$BUGSoutput$sims.array[, , "Sigma[2,1,1]"]
sigma_post1[1,2, ] <- MOVE$BUGSoutput$sims.array[, , "Sigma[1,2,1]"]
sigma_post1[2,2, ] <- MOVE$BUGSoutput$sims.array[, , "Sigma[2,2,1]"]
sigma_post2 <- array(NA, c(2, 2, length(MOVE$BUGSoutput$sims.array[, , "Sigma[1,1,2]"])))
sigma_post2[1,1, ] <- MOVE$BUGSoutput$sims.array[, , "Sigma[1,1,2]"]
sigma_post2[2,1, ] <- MOVE$BUGSoutput$sims.array[, , "Sigma[2,1,2]"]
sigma_post2[1,2, ] <- MOVE$BUGSoutput$sims.array[, , "Sigma[1,2,2]"]
sigma_post2[2,2, ] <- MOVE$BUGSoutput$sims.array[, , "Sigma[2,2,2]"]
sigma_postl <- list(sigma_post1, sigma_post2)

# load MASS for random multivariate normal variables
library(MASS)

# we do not estimate a values for b in the last day of the season, so do not allow when getting posterior predictions
Xidx3[Xidx3==season_length] <- season_length - 1
# create plot comparing, for each species, observed locations to posterior predictions 
plot(1, 1, ylim=c(10, 60), xlim=c(-105, -45), col="white", xlab="Longitude", ylab="Latitude")
# posterior predictions for maximum longitude 
max_long <- mat.or.vec(100, N)
# posterior predictions for the mean longitude
mean_long <- mat.or.vec(100, N)
# posterior predictions for the max distance traveled (by latitude)
lat_dist <- mat.or.vec(100, N)
# create an array to store observed locations (individual x days since the start of the season x lat/long)
obs <- array(NA, c(N, season_length, 2))
# create a vector to store the number of recorded locations for each individual
num_recs <- mat.or.vec(N, 1)
cdf <- list()
# for N individuals
for(i in 1:N){
  # create an array to store posterior predictions (iterations x season length x lat/long)
  randloc <- array(NA, c(100, season_length, 2))
  # how many records for individual i
  num_recs[i] <- length(allindvs_ordered_filtered_nodups$days_since[Yidx[i]:(Yidx[i+1]-1)])
  # fill in array of NAs with recorded locations
  obs[i, allindvs_ordered_filtered_nodups$days_since[Yidx[i]:(Yidx[i+1]-1)], 1] <- y[Yidx[i]:(Yidx[i+1]-1), 1]
  obs[i, allindvs_ordered_filtered_nodups$days_since[Yidx[i]:(Yidx[i+1]-1)], 2] <- y[Yidx[i]:(Yidx[i+1]-1), 2]
  for(z in 1:100){
    # draw from iteration from posterior
    gamma <- mat.or.vec(2, 1)
    gamma[1] <- gamma_post[1, z]
    gamma[2] <- gamma_post[2, z]
    D <- mat.or.vec(2, 2)
    D[1,1] <- D_post[1,1,z]
    D[1,2] <- D_post[1,2,z]
    D[2,1] <- D_post[2,1,z]
    D[2,2] <- D_post[2,2,z]
    sigma1 <- sigma_post1[, , z]
    sigma2 <- sigma_post2[, , z]
    sigma <- list(sigma1, sigma2)
    #sigma <- solve(sigma)
    b <- mat.or.vec(season_length, 1)
    for(s in 1:season_length){
      b[s] <- b_post[i, s]
    }
    # draw a multivariate random variable using posterior values for parameters
    # start with the first recorded locations y[Yidx] at t = 0
    # t = 1 does not use the autoregressive components yet
    randloc[z, (Xidx[i]-season_length*(i-1)), ] <- as.numeric(y[Yidx[i], ])
    randloc[z, (Xidx[i]-season_length*(i-1)+1), ] <- mvrnorm(n=1, mu=as.numeric(y[Yidx[i], ]), Sigma=sigma[[b[(Xidx[i]-season_length*(i-1)+1)]]])
    # create an object to store displacement in the x and y dimension
    disp <- mat.or.vec(2, 1)
    disp[1] <- randloc[z, (Xidx[i]-season_length*(i-1)+1), 1] + max(-1, min((randloc[z, (Xidx[i]-season_length*(i-1)+1), 1] - y[Yidx[i], 1])*gamma[b[(Xidx[i]-season_length*(i-1)+2)]] + D[1,b[(Xidx[i]-season_length*(i-1)+2)]], 1))
    disp[2] <- randloc[z, (Xidx[i]-season_length*(i-1)+1), 2] + max(-1, min((randloc[z, (Xidx[i]-season_length*(i-1)+1), 2] - y[Yidx[i], 2])*gamma[b[(Xidx[i]-season_length*(i-1)+2)]] + D[2,b[(Xidx[i]-season_length*(i-1)+2)]], 2))
    # t = 2 uses the autoregressive component by referencing the first recorded location
    randloc[z, (Xidx[i]-season_length*(i-1)+2), ] <- mvrnorm(1, disp[], sigma[[b[(Xidx[i]-season_length*(i-1)+1)]]])
    for(t in (Xidx[i]-season_length*(i-1)+3):(season_length-1)){
      disp <- mat.or.vec(2, 1)
      disp[1] <- randloc[z, t-1, 1] +  max(-1, min((randloc[z, t-1, 1] - randloc[z, t-2, 1])*gamma[b[t]] + D[1, b[t]], 1))
      disp[2] <- randloc[z, t-1, 2] + max(-1, min((randloc[z, t-1, 2] - randloc[z, t-2, 2])*gamma[b[t]] + D[2, b[t]], 1))
      randloc[z, t, ] <- mvrnorm(1, disp[], sigma[[b[t]]])
      #if(randloc[z, t, 2] < lon_min){
      #  randloc[z, t, 2] <- lon_min
      #}
      #if(randloc[z, t, 2] > lon_max){
      #  randloc[z, t, 2] <- lon_max
      #}
      #if(randloc[z, t, 1] < lat_min){
      #  randloc[z, t, 1] <- lat_min
      #}
      #if(randloc[z, t, 1] > lat_max){
      #  randloc[z, t, 1] <- lat_max
      #}
      #allindvs_ordered_filtered_nodups$days_since  
    }
    points(randloc[z,allindvs_ordered_filtered_nodups$days_since[Yidx[i]:(Yidx[i+1]-1)] ,2], randloc[z,allindvs_ordered_filtered_nodups$days_since[Yidx[i]:(Yidx[i+1]-1)] ,1], col=rgb(0, 0, 0, 0.2), pch=16)
    # calculate a test statistic from posterior predictions, z, using only time steps that have recorded locations in the observed data araay
    max_long[z, i] <- mean(randloc[z, (Xidx[i]-season_length*(i-1)):(Xidx3[i]), 1], na.rm=TRUE)
    mean_long[z, i] <- mean(randloc[z, (Xidx[i]-season_length*(i-1)):(Xidx3[i]), 2], na.rm=TRUE)
    lat_dist[z, i] <- max(randloc[z,(Xidx[i]-season_length*(i-1)):(Xidx3[i]), 1]) - 
      min(randloc[z,(Xidx[i]-season_length*(i-1)):(Xidx3[i]), 1])
  }
  cdf[[i]] <- ecdf(randloc[ , allindvs_ordered_filtered_nodups$days_since[Yidx[i]:(Yidx[i+1]-1)] ,2])
  }
points(obs[2, , 2], obs[2, , 1], col=rgb(0, 0.5, 1, 1), pch=16)
lines(obs[2, , 2], obs[2, , 1], col=rgb(0, 0.5, 1, 1), pch=16)

# get a vector for the first recorded locations for each individual to sort posterior predictions by latitude
first_loc <- mat.or.vec(N, 1)
for(k in 1:N){
  first_loc[k] <- y[Yidx[k], 1]
}


# 3-panel plot for posterior checks using 3 test statistics 
plot_width <- 2
species_lab <- "Piping Plover"
quartz.options(width=8, height=5)
layout(matrix(c(1, 2, 3), 1, 3, byrow = TRUE))
par(mar=c(4.5, 3.5, 2, 1))
plot(1, 1, col="white", xlim=c(min(apply(obs[1:N,,1], MARGIN=1, FUN=mean, na.rm=TRUE)[num_recs>= 1]) - plot_width, max(apply(obs[1:N,,1], MARGIN=1, FUN=mean, na.rm=TRUE)[num_recs>= 1]) + plot_width), 
     ylim=c(0.5, N+0.5), ylab=" ", xlab=" ", main = species_lab, cex.axis=0.75, cex.main=0.85, bty="n")
for(i in 1:N){
  if(num_recs[i]>=1){
    segments(quantile(max_long[, order(first_loc)][, i], c(0.025)), i, quantile(max_long[, order(first_loc)][, i], c(0.975)), lend="butt", col="gray", lwd=4)
  }
}
points(apply(obs[1:N,,1], MARGIN=1, FUN=mean, na.rm=TRUE)[order(first_loc)][num_recs>= 1], (1:N)[num_recs>= 1], pch=16, cex=0.6)
abline(v = mean(apply(obs[1:N,,1], MARGIN=1, FUN=mean, na.rm=TRUE)[num_recs>= 1]))
mtext(side=1, line=2,"Mean latitude (dec. deg.)", cex=0.75)
mtext(side=2, line=2,"Individual (from lowest to highest latitude)", cex=0.75)


plot(1, 1, col="white", xlim=c(min(apply(obs[1:N,,2], MARGIN=1, FUN=mean, na.rm=TRUE)[num_recs>= 1]) - plot_width, max(apply(obs[1:N,,2], MARGIN=1, FUN=mean, na.rm=TRUE)[num_recs>= 1]) + plot_width), 
     ylim=c(0.5, N+0.5), ylab=" ", xlab=" ", main = species_lab, cex.axis=0.75, cex.main=0.85, bty="n")
for(i in 1:N){
  if(num_recs[i]>=1){
    segments(quantile(mean_long[, order(first_loc)][, i], c(0.025)), i, quantile(mean_long[, order(first_loc)][, i], c(0.975)), lend="butt", col="gray", lwd=4)
  }
}
points(apply(obs[1:N,,2], MARGIN=1, FUN=mean, na.rm=TRUE)[order(first_loc)][num_recs>= 1], (1:N)[num_recs>= 1], pch=16, cex=0.6)
mtext(side=1, line=2,"Mean longitude (dec. deg.)", cex=0.75)
mtext(side=2, line=2,"Individual (from lowest to highest latitude)", cex=0.75)
abline(v = mean(apply(obs[1:N,,2], MARGIN=1, FUN=mean, na.rm=TRUE)[num_recs>= 1]))


plot(1, 1, col="white", xlim=c(0, max(apply(obs[1:N,,1], MARGIN=1, FUN=max, na.rm=TRUE)[num_recs>= 1] - apply(obs[1:N,,1], MARGIN=1, FUN=min, na.rm=TRUE)[num_recs>= 1]) + plot_width), 
     ylim=c(0.5, N+0.5), ylab=" ", xlab=" ", main = species_lab, cex.axis=0.75, cex.main=0.85, bty="n")
for(i in 1:N){
  if(num_recs[i]>=1){
    segments(quantile(lat_dist[, order(first_loc)][, i], c(0.025)), i, quantile(lat_dist[, order(first_loc)][, i], c(0.975)), lend="butt", col="gray", lwd=4)
  }
}
points(apply(obs[1:N,,1], MARGIN=1, FUN=max, na.rm=TRUE)[order(first_loc)][num_recs>= 1] - apply(obs[1:N,,1], MARGIN=1, FUN=min, na.rm=TRUE)[order(first_loc)][num_recs>= 1], (1:N)[num_recs>= 1], pch=16, cex=0.6)
mtext(side=1, line=2,"Distance traveled", cex=0.75)
mtext(side=1, line=3.1,"(latitude in dec. deg.)", cex=0.75)
mtext(side=2, line=2,"Individual (from lowest to highest latitude)", cex=0.75)
abline(v = mean(apply(obs[1:N,,1], MARGIN=1, FUN=max, na.rm=TRUE)[num_recs>= 1] - apply(obs[1:N,,1], MARGIN=1, FUN=min, na.rm=TRUE)[num_recs>= 1]))


# create a plot that for each individual compares the cdf of the posterior to the observed point
# get a vector of observed test statistics for each individual; here shown for mean longitude 
obs_points <- apply(obs[1:N,,2], MARGIN=1, FUN=mean, na.rm=TRUE)
pit <- mat.or.vec(N, 1)
for(i in 1:N){
  # cdf for the posterior for each individual was estimated above in the test statistic code
  pit[i] <- cdf[[i]](obs_points[i])
}
hist(pit, breaks=50, bty="n", main=" ", xlab="Position of observed in posterior CDF", ylim=c(0, 6))
