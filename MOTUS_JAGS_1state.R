##################################################################################################################################
### JAGS model and associated posterior checks with plots
##################################################################################################################################

# general model structure modified from Baldwin et al. 2018
MOVE <- function(){
  #trunc[1, 1] <- lat_min
  #trunc[1, 2] <- lat_max
  #trunc[2, 1] <- lon_min
  #trunc[2, 2] <- lon_max
  
  sd ~ dunif(0, 0.0001)
  tau <- 1/(sd*sd)
  
  # variance-covariance matrix
  sigma_c[1, 1] ~ dexp(1) # it cannot be negative
  sigma_c[1, 2] <- 0
  sigma_c[2, 1] <- 0
  sigma_c[2, 2] ~ dexp(1)
  
  Rho[1, 1] <- 1.0
  Rho[1, 2] ~ dunif(-0.8, 0.8)
  Rho[2, 1] <- Rho[1, 2]
  #Rho[1, 2] <- 0
  #Rho[2, 1] <- 0
  Rho[2, 2] <- 1.0
  
  Sigma = sigma_c %*% Rho %*% sigma_c
  
  # get the first location for each indvidual
  for(k in 1:N){
    first.loc[k, 1] <- y[Yidx[k], 1]
    first.loc[k, 2] <- y[Yidx[k], 2]
  }
  
  # modify the model so it has only 1 state to choose from
  gamma ~ dunif(0, 1)
  
  # priors for drift component (1 for lat; 2 for lon)
  D[1] ~ dnorm(0, 0.1)
  D[2] ~ dnorm(0, 0.1)
  #D[1] <- 0
  #D[2] <- 0
  
  # N is the individual index
  for(k in 1:N){
    # mute measurement error components
    #logpsi[k] ~ dunif(-10, 10)
    #psi[k] <- exp(logpsi[k])
    
    # j is the index for lat and lon
    for(j in 1:2){
      # mute measurement error components
      #x[Xidx[k], j] ~ dt(first.loc[k, j], itau2[Xidx[k], j]*psi[k], nu[Xidx[k], j])
      # get x for the first location
      x[Xidx[k], j] <- first.loc[k, j]
    }
    # get x for the second recorded location, which does not use the autoregressive component
    x[(Xidx[k]+1), 1:2] ~ dmnorm.vcov(x[Xidx[k],], Sigma[, ])
    
    # mute loglikelihood monitor for first two locations
    #log_lik[(Xidx[k])] <- logdensity.mnorm(x[(Xidx[k]), 1:2], x[(Xidx[k]), ], iSigma[, ]) 
    #log_lik[(Xidx[k]+1)] <- logdensity.mnorm(x[(Xidx[k]+1), 1:2], x[(Xidx[k]+1), ], iSigma[, ])
    # for regular time steps, t
    for(t in (Xidx[k]+1):(Xidx2[k]-1)){ 
      # displacement at time t + 1, from t, includes an autoregressive component, drift, and a max hop distance
      # for lat and lon 
      displace[t, 1] <- max(-10, min((x[t, 1] - x[t-1, 1])*gamma + D[1], 10))
      displace[t, 2] <- max(-10, min((x[t, 2] - x[t-1, 2])*gamma + D[2], 10))
      # combine lat and lon and add displacement to the previous time step, t, to feed into dmnorm()
      x.mn[t, 1:2] <- x[t, 1:2] + displace[t, 1:2]
      # get x for t + 1
      x[t+1, 1:2] ~ dmnorm.vcov(x.mn[t, ], Sigma[, ])
    }
    
    for(i in (Yidx[k]+1):(Yidx[k+1]-1)) { 
      for(j in 1:2) {
        # mute the component that interpolates multiple detections within a single time step
        #yhat[i, j] <- w[i]*x[idx[i], j] + (1 - w[i])*x[idx[i+1], j] 
        yhat[i, j] <- x[idx[i], j]
        # mute the measurement error component
        #y[i, j] ~ dt(yhat[i, j], itau2[i, j]*psi[k], nu[i, j])
        y[i, j] ~ dnorm(yhat[i, j], tau);
        #T(trunc[j, 1], trunc[j, 2])
      }
    }
  }
}

if (is.R()){
  filename <- file.path(tempdir(), "MOVE.bug")}
write.model(MOVE, filename)
inits <- list(list(gamma=c(0.5), D=c(0, 0))) #x=x
data <- list("Xidx", "Xidx2", "Yidx", "y", "idx","N")
             #"lat_min", "lat_max", "lon_min", "lon_max")
parameters <- c("gamma", "Sigma", "D", "x") 
MOVE <- jags(data=data, inits=inits, parameters.to.save=parameters, filename,
             n.chains=1, n.burnin=n.burnin, n.iter=n.iter+n.burnin, n.thin=1, DIC=TRUE)
MOVE.mcmc <- as.mcmc(MOVE)

# create a simple plot to show predicted tracks alongside observed detections
# load shapefile of spatial units for plot
setwd("~/Documents/folders/URI/GIS/")
BOEM.sf <- st_read("BOEM_1deg_grid_latlon.shp")
bound <- getData('GADM', country='USA', level=1)
e <- c(min(BOEM$xmin)-1, max(BOEM$xmax)+1, min(BOEM$ymin)-1, max(BOEM$ymax)+1)
boundcrop <- crop(bound, e)
bound.simplify <- gSimplify(boundcrop, tol=0.01, topologyPreserve=TRUE)

par(mar=c(4, 4, 2, 1))
#(level=0 is country boundary, level=1 is state boundaries, level=2 is county boundaries)
plot(0, 0, xlim=c(e[1], e[2]), ylim=c(e[3], e[4]), xaxt="n", yaxt="n", yaxs="i", xaxs="i", xlab=" ", ylab=" ", xpd=FALSE, main=" ", cex.main=0.75)
plot(BOEM.sf, add=TRUE, border="black", col="white")
plot(bound.simplify, add=TRUE, col="white", border="dark gray", lwd=.9)
axis(side=1, line=-3, cex.axis=0.75)
axis(side=2, line=-3, cex.axis=0.75)
mtext(side=1, line=-1, "Longitude", cex=0.75)
mtext(side=2, line=-0.5, "Latitude", cex=0.75)

arrow_colors <- c(rgb(0, 0, 1, 0.5), rep(rgb(0, 0, 0, 0), 10))
lat_mu <- mat.or.vec(N*season_length, 1)
lon_mu <- mat.or.vec(N*season_length, 1)
lat_sd <- mat.or.vec(N*season_length, 1)
ind <- 50
for(i in ((ind-1)*season_length + 1):((ind-1)*season_length + season_length)){
  lat_mu[i] <- mean(MOVE$BUGSoutput$sims.array[, , paste("x[", i, ",1]", sep="")])
  lat_sd[i] <- sd(MOVE$BUGSoutput$sims.array[, , paste("x[", i, ",1]", sep="")])
  lon_mu[i] <- mean(MOVE$BUGSoutput$sims.array[, , paste("x[", i, ",2]", sep="")])
}
arrows(lon_mu[(Xidx[ind]):((ind-1)*season_length + season_length -1)], lat_mu[(Xidx[ind]):((ind-1)*season_length + season_length -1)], 
       lon_mu[(Xidx[ind] + 1):((ind-1)*season_length + season_length -1)+1], lat_mu[(Xidx[ind]):((ind-1)*season_length + season_length -1)+1], 
       lend="butt", length = 0.05, col = arrow_colors)
lines(lon_mu[(Xidx[ind]):((ind-1)*season_length + season_length)], lat_mu[Xidx[ind]:((ind-1)*season_length + season_length)], type="l", col="blue")
for(i in (Yidx[ind]):(Yidx[ind+1]-1)) { 
  points(y[i, 2], y[i, 1], col=rgb(0.1, 0.1, 0.1, 1), pch=16, cex=(0.9))
}
#y[(Yidx[ind]):(Yidx[ind+1]-2), ]

## save MCMC output to objects to use for posterior checks #######################################################################
gamma_post <- MOVE$BUGSoutput$sims.array[, , "gamma"]
D_post <- mat.or.vec(2, length(MOVE$BUGSoutput$sims.array[, , "D[1]"]))
D_post[1, ] <- MOVE$BUGSoutput$sims.array[, , "D[1]"]
D_post[2, ] <- MOVE$BUGSoutput$sims.array[, , "D[2]"]
sigma_post <- array(NA, c(2, 2, length(MOVE$BUGSoutput$sims.array[, , "Sigma[1,1]"])))
sigma_post[1,1, ] <- MOVE$BUGSoutput$sims.array[, , "Sigma[1,1]"]
sigma_post[2,1, ] <- MOVE$BUGSoutput$sims.array[, , "Sigma[2,1]"]
sigma_post[1,2, ] <- MOVE$BUGSoutput$sims.array[, , "Sigma[1,2]"]
sigma_post[2,2, ] <- MOVE$BUGSoutput$sims.array[, , "Sigma[2,2]"]


## create plot comparing, for each species, observed locations to posterior predictions ##########################################
# load MASS for random multivariate normal variables
library(MASS)

# we do not estimate a values for b in the last day of the season, so do not allow when getting posterior predictions
Xidx3[Xidx3==season_length] <- season_length - 1
# posterior predictions for maximum longitude 
max_long <- mat.or.vec(100, N)
# posterior predictions for the mean longitude
mean_long <- mat.or.vec(100, N)
# posterior predictions for the max distance traveled (by latitude)
lat_dist <- mat.or.vec(100, N)
# create an array to store observed locations (individual x days since the start of the season x lat/lon)
obs <- array(NA, c(N, season_length, 2))
# create a vector to store the number of recorded locations for each individual
num_recs <- mat.or.vec(N, 1)
cdf <- list()
# for N individuals
for(i in 1:N){
  # create an array to store posterior predictions (iterations x season length x lat/lon)
  randloc <- array(NA, c(100, season_length, 2))
  # how many records for individual i
  num_recs[i] <- length(allindvs_ordered_filtered_nodups$days_since[Yidx[i]:(Yidx[i+1]-1)])
  # fill in array of NAs with recorded locations
  obs[i, allindvs_ordered_filtered_nodups$days_since[Yidx[i]:(Yidx[i+1]-1)], 1] <- y[Yidx[i]:(Yidx[i+1]-1), 1]
  obs[i, allindvs_ordered_filtered_nodups$days_since[Yidx[i]:(Yidx[i+1]-1)], 2] <- y[Yidx[i]:(Yidx[i+1]-1), 2]
  for(z in 1:100){
    # draw from iteration from posterior
    gamma <- gamma_post[z]
    D <- mat.or.vec(2, 1)
    D[1] <- D_post[1, z]
    D[2] <- D_post[2, z]
    sigma <- sigma_post[,,z]
    # draw a multivariate random variable using posterior values for parameters
    # start with the first recorded locations y[Yidx] at t = 0
    # t = 1 does not use the autoregressive components yet
    randloc[z, (Xidx[i]-season_length*(i-1)), ] <- as.numeric(y[Yidx[i], ])
    randloc[z, (Xidx[i]-season_length*(i-1)+1), ] <- mvrnorm(n=1, mu=as.numeric(y[Yidx[i], ]), Sigma=sigma)
    # create an object to store displacement in the x and y dimension
    disp <- mat.or.vec(2, 1)
    disp[1] <- randloc[z, (Xidx[i]-season_length*(i-1)+1), 1] + max(-1, min((randloc[z, (Xidx[i]-season_length*(i-1)+1), 1] - y[Yidx[i], 1])*gamma + D[1], 1))
    disp[2] <- randloc[z, (Xidx[i]-season_length*(i-1)+1), 2] + max(-1, min((randloc[z, (Xidx[i]-season_length*(i-1)+1), 2] - y[Yidx[i], 2])*gamma + D[2], 2))
    # t = 2 uses the autoregressive component by referencing the first recorded location
    randloc[z, (Xidx[i]-season_length*(i-1)+2), ] <- mvrnorm(1, disp[], sigma)
    for(t in (Xidx[i]-season_length*(i-1)+3):(season_length-1)){
      disp <- mat.or.vec(2, 1)
      disp[1] <- randloc[z, t-1, 1] +  max(-1, min((randloc[z, t-1, 1] - randloc[z, t-2, 1])*gamma + D[1], 1))
      disp[2] <- randloc[z, t-1, 2] + max(-1, min((randloc[z, t-1, 2] - randloc[z, t-2, 2])*gamma + D[2], 1))
      randloc[z, t, ] <- mvrnorm(1, disp[], sigma)
    }
    #points(randloc[z,allindvs_ordered_filtered_nodups$days_since[Yidx[i]:(Yidx[i+1]-1)] ,], col=rgb(0, 0, 0, 0.2), pch=16)
    # calculate a test statistic from posterior predictions, z, using only time steps that have recorded locations in the observed data araay
    max_long[z,i] <- mean(randloc[z, (Xidx[i]-season_length*(i-1)):(Xidx3[i]), 1])
    mean_long[z,i] <- mean(randloc[z, (Xidx[i]-season_length*(i-1)):(Xidx3[i]), 2])
    lat_dist[z, i] <- max(randloc[z, (Xidx[i]-season_length*(i-1)):(Xidx3[i]), 1]) - 
      min(randloc[z, (Xidx[i]-season_length*(i-1)):(Xidx3[i]),1])
  }
  cdf[[i]] <- ecdf(randloc[ , allindvs_ordered_filtered_nodups$days_since[Yidx[i]:(Yidx[i+1]-1)] ,2])
}

# get a vector for the first recorded locations for each individual to sort posterior predictions by latitude
first_loc <- mat.or.vec(N, 1)
for(k in 1:N){
  first_loc[k] <- y[Yidx[k], 1]
}

# 3-panel plot for posterior checks using 3 test statistics 
plot_width <- 10
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
points((apply(obs[1:N,,1], MARGIN=1, FUN=max, na.rm=TRUE)[num_recs>= 1] - apply(obs[1:N,,1], MARGIN=1, FUN=min, na.rm=TRUE)[num_recs>= 1])[order(first_loc)], (1:N)[num_recs>= 1], pch=16, cex=0.6)
mtext(side=1, line=3.1,"(latitude in dec. deg.)", cex=0.75)
mtext(side=2, line=2,"Individual (from lowest to highest latitude)", cex=0.75)
abline(v = mean(apply(obs[1:N,,1], MARGIN=1, FUN=max, na.rm=TRUE)[num_recs>= 1] - apply(obs[1:N,,1], MARGIN=1, FUN=min, na.rm=TRUE)[num_recs>= 1]))


# create a plot that for each individual compares the cdf of the posterior to the observed point
# get a vector of observed test statistics for each individual; here shown for mean longitude 
quartz.options(width=5, height=4)
obs_points <- apply(obs[1:N,,2], MARGIN=1, FUN=mean, na.rm=TRUE)
pit <- mat.or.vec(N, 1)
for(i in 1:N){
  # cdf for the posterior for each individual was estimated above in the test statistic code
  pit[i] <- cdf[[i]](obs_points[i])
}
hist(pit, breaks=50, bty="n", main=" ", xlab="Position of observed in posterior CDF", ylim=c(0, 6))
