## load required packages
library(R2jags)
library(R2WinBUGS)
library(sp)
library(raster)
library(rgeos)
library(sf)
library(MASS)

# read a file that has the extent of user-defined spatial units as rows, labeled "xmin", "ymax", "xmax", "ymin"
# can use the field calculator in QGIS to get a data file with the extent of each cell and centroid coordinates
# e.g. x_min($geometry) on polygons and $x on a .shp for centroids
setwd("~/Documents/folders/GitProjects/URI_CRM/motus_data/")
BOEM <- read.csv('BOEM_1deg_XY.csv', header=TRUE)
# label columns
colnames(BOEM)[1] <- "xmin"
colnames(BOEM)[2] <- "ymax"
colnames(BOEM)[3] <- "xmax"
colnames(BOEM)[4] <- "ymin"

# load shapefile of spatial units for plot
setwd("~/Documents/folders/URI/GIS/")
BOEM.sf <- st_read("BOEM_1deg_grid_latlon.shp")

setwd("~/Documents/folders/URI/GIS/")
motus_rec <- st_read("receiver-deployments_BOEM_XY.shp")

# get centroids for polygons to calculate distance of each poly to nearest motus receiver
BOEM_center <- st_centroid(BOEM.sf)
BOEM_center_WGS84 <- st_transform(BOEM_center, st_crs(motus_rec))

# distance of each poly to each motus receiver
pt_dists <- st_distance(BOEM_center_WGS84, motus_rec)
min_motus <- mat.or.vec(length(pt_dists[,1]), 1)
# find closest receiver and covert distance to km
for(i in 1:length(pt_dists[,1])){
  min_motus[i] <- pt_dists[i, which(pt_dists[i,]==min(pt_dists[i,]))[1]]/1000
}

# initiate plot of randomly generated tracks
#quartz.options(width=7.48, height=5)
#layout(matrix(c(1, 2), 1, 2, byrow = TRUE))
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
#plot(BOEM_center, add=T)

# begin primary script to generate tracks and run simulation using known values
# simulation iterations
M <- 5
num_ind <- 10
params_out <- list()
out_bounds_array <- array(0, c(M, length(min_motus), 4))
post_pred_obs <- array(0, c(length(min_motus), M, 4))
post_pred_mod <- array(0, c(length(min_motus), M, 4))
post_pred_mod_check <- array(0, c(length(min_motus), M, 4))
post_pred_ecdf <- array(0, c(length(min_motus), M, 4))
params <- array(0, c(7, 4, M))
# mute1 determines whether the locations will be filtered by exisiting motus network (1 is filtered)
mute1 <- c(0, 0, 1, 1)
# mute2 determines if locations will be randomly assigned missing values (1 is random NA)
mute2 <- c(0, 1, 0, 0)
# mute3 determines whether individuals will take big flights into open ocean
mute3 <- c(0, 0, 0, 1)
# for p scenarios
for(p in 1:4){
  out_bounds <- mat.or.vec(M, length(min_motus))
  for(m in 1:M){
    rm(MOVE)
    num_day <- 100
    # values are divided by 10 to covert tracks that work for 10 time steps to 100 time steps
    XD <- array(0, c(num_day, 2, num_ind))
    #Sigma <- matrix(c(0.25, 0.03125, 0.03125, 1), nrow = 2)
    Sigma <- matrix(c(0.01, 0, 0, 0.01), nrow = 2)
    if(mute3[p] == 0){
    state <- rep(1, num_day)
    }else{
    state <- c(rep(1, num_day/2), rep(2, num_day/2))
    }
    # individual
    for(z in 1:num_ind){
      XD[1,,z] <- c(-70.242236, 41.242236)
      #XD[2,,z] <- c(rnorm(1, -70.242236, 0.1), runif(1, 41, 41.242236))
      XD[2,,z] <- mvrnorm(n = 1, XD[1,,z],  Sigma = Sigma)
      # time step
      for(i in 2:(num_day-1)){
        XD2x <- max(-100, min((XD[i, 1, z] - XD[i-1, 1, z])*0.2 - c(0.11, -0.3)[state[i]], 100)) + XD[i, 1, z]
        XD2y <- max(-100, min((XD[i, 2, z] - XD[i-1, 2, z])*0.2 - c(0.11, 0.05)[state[i]], 100)) + XD[i, 2, z] 
        XD[i+1,,z] <- mvrnorm(n = 1, c(XD2x, XD2y),  Sigma = Sigma)
      }
    }
    
    # plot randomly generated tracks
    #for(z in 1:num_ind){
    #for(i in 1:(num_day-1)){
    ##points(XD[i, 1, z], XD[i, 2, z], pch=16, col = rgb(0, 0, 1, 0.3))
    #arrows(XD[i, 1, z], XD[i, 2, z], XD[i+1, 1, z], XD[i+1, 2, z], lend="butt", length = 0, col = rgb(0, 0, 1, 0.4))
    #}
    ##points(XD[10, 1, z], XD[10, 2, z], pch = 16, col = rgb(0, 0, 1, 0.3))
    #}
    
    if(mute1[p] == 1){
      # match the locations for each track with any motus towers within 20 km
      loc_match <- array(0, c(num_day, 2, num_ind))
      # individual
      for(s in 1:num_ind){
        # location along track
        for(i in 1:num_day){
          loc_match_dists <- unlist(lapply(1:length(motus_rec$latitude), function(x) sqrt((motus_rec$latitude[x] -                                                                            XD[i, 2, s])^2 + (motus_rec$longitude[x] - XD[i, 1, s])^2)))
          which_min_dist <- which(unlist(loc_match_dists) == min(unlist(loc_match_dists)))[1]
          # convert dd to km
          min_dist <- loc_match_dists[which_min_dist]*111
          if(min_dist < 20){
            loc_match[i, ,s] <- XD[i, ,s]
          }else{
            loc_match[i, ,s] <- NA
          }
        }
      }
    }
    
    # mute next line when the NAs according to motus has already been run
    if(mute1[p]==0){
      loc_match <- XD
    }
    day_index <- 1:num_day
    #max_obs <- mat.or.vec(length(XD[1,1,]), 1)
    min_obs <- mat.or.vec(length(XD[1,1,]), 1)
    loc_match_indexed <- c(NA, NA)
    loc_match_indexed <- rbind(loc_match_indexed, c(NA, NA))
    Sind_obs <- 1
    for(i in 1:length(XD[1,1,])){
      tempsample <- sample(2:num_day, 50, replace=FALSE)
      if(mute2[p]==1){
      loc_match[tempsample, 1, i] <- NA
      loc_match[tempsample, 2, i] <- NA
      }
      #max_obs[i] <- max(which(!is.na(loc_match[, 1, i])))
      min_obs[i] <- min(which(!is.na(loc_match[, 1, i])))
      Sind_obs <- c(Sind_obs, (length(loc_match_indexed[,1]) - 1))
      loc_match_indexed <- rbind(loc_match_indexed, loc_match[!is.na(loc_match[, 1, i]), ,i])
    }
    loc_match_indexed <- loc_match_indexed[c(-1, -2),]
    Sind_obs <- Sind_obs[-1]
    
    Xidx <- mat.or.vec(length(XD[1,1,]), 1)
    for(i in 1:length(XD[1,1,])){
      Xidx[i] <- min_obs[i] + num_day*(i-1)
    }
    
    idx <- (1:(num_day*num_ind))[which(!is.na(as.vector(loc_match[,1,])))]
    # locations
    y <- loc_match_indexed
    # number of individuals
    N <- num_ind
    # there is an additional position added to the end of the vector to ensure that the index loops stay in bounds (as in Baldwin et al. 2018)
    Xidx <- c(Xidx, (num_day*N+1))
    # there is an additional position added to the end of the vector to ensure that the index loops stay in bounds (as in Baldwin et al. 2018)
    Yidx <- c(Sind_obs, (length(y[,1])+1))
    
    x1 <- rep(mean(y[,1]), num_day*N)
    x2 <- rep(mean(y[,2]), num_day*N)
    x <- cbind(x1, x2)
    for(i in 1:(length(Xidx)-1)){
      x[Xidx[i],] <-NA
    }
    
    # number of iterations for the JAGS model
    n.iter <- 5000
    n.burnin <- 10000
    
    sigma_c <- matrix(c(0.01, NA, NA, 0.01), 2, 2)
    
    if(mute3[p] != 9){
    MOVE <- function(){
      # get the first location for each indvidual
      for(k in 1:N){
        first.loc[k, 1] <- y[Yidx[k], 1]
        first.loc[k, 2] <- y[Yidx[k], 2]
      }
      
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
      
      # modify the model so it has only 1 state to choose from
      gamma ~ dunif(0, 1)
      
      # priors for drift component (1 for lat; 2 for lon)
      D[1] ~ dnorm(0, 0.1)
      D[2] ~ dnorm(0, 0.1)
      #D[1] <- 0
      #D[2] <- 0
      
      #for(r in 1:(max(Xidx)-1)){
      #lile[r, 1:2] ~ dmnorm.vcov(c(0, 0), Sigma[, ])
      #}
      
      # N is the individual index
      for(k in 1:N){
        # j is the index for lat and lon
        for(j in 1:2){
          x[Xidx[k], j] <- first.loc[k, j]
        }
        # get x for the second recorded location, which does not use the autoregressive component
        #x[(Xidx[k]+1), 1:2] <- x[Xidx[k],] + lile[Xidx[k]+1,]
        x[(Xidx[k]+1), 1:2] ~ dmnorm.vcov(x[Xidx[k],], Sigma[1:2, 1:2])
        # for regular time steps, t
        for(t in (Xidx[k]+1):(Xidx[k+1]-2)){ 
          # displacement at time t + 1, from t, includes an autoregressive component, drift, and a max hop distance
          # for lat and lon 
          displace[t, 1] <- max(-100, min((x[t, 1] - x[t-1, 1])*gamma + D[1], 100))
          displace[t, 2] <- max(-100, min((x[t, 2] - x[t-1, 2])*gamma + D[2], 100))
          # combine lat and lon and add displacement to the previous time step, t, to feed into dmnorm()
          x.mn[t, 1:2] <- x[t, 1:2] + displace[t, 1:2]
          # get x for t + 1
          #x[t+1, 1:2] <- x.mn[t, ] + lile[t+1, ]
          x[t+1, 1:2] ~ dmnorm.vcov(x.mn[t, ], Sigma[1:2, 1:2])
        }
        for(i in (Yidx[k]+1):(Yidx[k+1]-2)){ 
          for(j in 1:2) {
            yhat[i, j] <- x[idx[i], j]
            y[i, j] ~ dnorm(yhat[i, j], tau)
          }
        }
      }
    }
    
    if (is.R()){
      filename <- file.path(tempdir(), "MOVE.bug")}
    write.model(MOVE, filename)
    inits <- list(list(gamma=c(0.5), D=c(0, 0), sigma_c=sigma_c, x=x)) # x=x
    data <- list("Xidx", "Yidx", "y", "idx","N")
    parameters <- c("gamma", "Sigma", "D", "x", "Rho", "sigma_c", "sd") 
    MOVE <- jags(data=data, inits=inits, parameters.to.save=parameters, filename,
                 n.chains=1, n.burnin=n.burnin, n.iter=n.iter+n.burnin, n.thin=1, DIC=TRUE)
    }else{
      MOVE <- function(){
        sd ~ dunif(0, 0.0001)
        tau <- 1/(sd*sd)
        
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
        u ~ dunif(0, 1)
        gamma[2] <- gamma[1]*u
        
        # priors for drift component (one for x and one for y dim.)
        D[1, 1] ~ dnorm(0, 0.1)
        D[1, 2] ~ dnorm(0, 0.1)
        D[2, 1] ~ dnorm(0, 0.1)
        D[2, 2] ~ dnorm(0, 0.1)
        
        alpha[1] ~ dunif(0, 1)
        alpha[2] ~ dunif(0, 1)
        lambda[1] ~ dunif(0, 1)
        lambda[2] <- 1 - lambda[1]
        
        # get the first location for each indvidual
        for(k in 1:N){
          first.loc[k, 1] <- y[Yidx[k], 1]
          first.loc[k, 2] <- y[Yidx[k], 2]
        }
        
        # N is the individual index
        for(k in 1:N){
          # estimate first behavioural state for each individual
          b[Xidx[k]] ~ dcat(lambda[])
          for(j in 1:2){
            x[Xidx[k], j] <- first.loc[k, j]
          }
          # get x for the second recorded location, which does not use autoregressive component
          x[(Xidx[k]+1), 1:2] ~ dmnorm.vcov(x[Xidx[k],], Sigma[, , b[Xidx[k]]])
          
          # for regular time steps, t
          for(t in (Xidx[k]+1):(Xidx[k+1]-2)){ 
            # displacement at time t + 1, from t, includes an autoregressive component, drift, and a max hop distance
            # for x and y dims. 
            phi[t, 1] <- alpha[b[t-1]] 
            phi[t, 2] <- 1 - alpha[b[t-1]]
            b[t] ~ dcat(phi[t, ])
            # combine x and y dims. and add displacement to the previous time step, t, to feed into dmnorm()
            displace[t, 1] <- max(-1, min((x[t, 1] - x[t-1, 1])*gamma[b[t]] + D[1, b[t]], 1))
            displace[t, 2] <- max(-1, min((x[t, 2] - x[t-1, 2])*gamma[b[t]] + D[2, b[t]], 1))
            # get x for t + 1
            x.mn[t, 1:2] <- x[t, 1:2] + displace[t, 1:2]
            x[t+1, 1:2] ~ dmnorm.vcov(x.mn[t, ], Sigma[, , b[t]])
          }
          
          zeta[k, 1] <- alpha[b[Xidx[k+1]-2]]
          zeta[k, 2] <- 1 - zeta[k, 1]
          b[Xidx[k+1]-1] ~ dcat(zeta[k, ])
          
          for(i in (Yidx[k]+1):(Yidx[k+1]-1)) { 
            for(j in 1:2) {
              yhat[i, j] <- x[idx[i], j]
              y[i, j] ~ dnorm(yhat[i, j], tau)
            }
          }
        }
      }
      
      if (is.R()){
        filename <- file.path(tempdir(), "MOVE.bug")}
      write.model(MOVE, filename)
      inits <- list(list(gamma=c(0.5, NA), x=x, D=matrix(c(0,0,0,0), 2, 2)))
      data <- list("Xidx", "Yidx", "y", "idx","N")
      parameters <- c("gamma", "Sigma", "D", "b", "x")
      MOVE <- jags(data=data, inits=inits, parameters.to.save=parameters, filename,
                   n.chains=1, n.burnin=n.burnin, n.iter=n.iter+n.burnin, n.thin=1, DIC=TRUE)
      MOVE.mcmc <- as.mcmc(MOVE)
    }
    
  
    ## use posterior estimates from JAGS model to obtain estimates for user-defined spatial units ####################################
    # convert output from JAGS to create arrays (one for lat; one for lon) that have the posterior samples (rows) 
    # for each time step of the movement model (columns), and for each individual (depth)
    posts_lat <- array(0, c(n.iter, num_day, N))
    posts_lon <- array(0, c(n.iter, num_day, N))
    for(e in 1:N){
      for(i in 1:num_day){
        posts_lat[,i,e] <- MOVE$BUGSoutput$sims.array[,,paste("x[", ((e-1)*num_day + i), "," , 2, "]", sep="")]
        posts_lon[,i,e] <- MOVE$BUGSoutput$sims.array[,,paste("x[", ((e-1)*num_day + i), "," , 1, "]", sep="")]
      }
    }
    posts_lat <- posts_lat[1:n.iter,,]
    posts_lon <- posts_lon[1:n.iter,,]
    
    # this script makes it possible to loop through every nth value
    # useful for running through a spatial input file that has too many rows to run wihtin a reasonable time 
    nth <- function(x,n){
      x[x%%n==0]
    }
    xs = 1:length(BOEM[,1])
    index <- nth(xs,1)
    
    # option 2b: extraplolation only includes uncertainty from estimating the proportion of 
    # individuals who likely crossed into the specified extent; gives estimates by year
    posts_latb <- posts_lat[1:n.iter, , ]
    posts_lonb <- posts_lon[1:n.iter, , ]
    setwd("~/Documents/folders/GitProjects/URI_CRM/")
    source('Pr_occupancy_season.R')
    occ_post <- mat.or.vec(n.iter, length(BOEM[,1]))
    occ_post_pred <- mat.or.vec(n.iter, length(BOEM[,1]))
    for(z in 1:length(index)){
      lat_input_min <- BOEM$ymin[z]
      lat_input_max <- BOEM$ymax[z]
      lon_input_min <- BOEM$xmin[z]
      lon_input_max <- BOEM$xmax[z]
      occ_post[, z] <- Pr_occupancy_season(posts_latb, posts_lonb, lat_input_min, lat_input_max, lon_input_min, lon_input_max, side, n.iter, N, num_day, "time_step")
      occ_post_pred[, z] <- rbinom(n.iter, num_ind, occ_post[, z])/num_ind
    }
    
    # option 2b: extraplolation only includes uncertainty from estimating the proportion of 
    # individuals who likely crossed into the specified extent; gives estimates by year
    posts_latb <- array(0, c(2, num_day, num_ind))
    posts_latb[1,,] <- XD[, 2, ]
    posts_latb[2,,] <- XD[, 2, ]
    posts_lonb <- array(0, c(2, num_day, num_ind))
    posts_lonb[1,,] <- XD[, 1, ]
    posts_lonb[2,,] <- XD[, 1, ]
    setwd("~/Documents/folders/GitProjects/URI_CRM/")
    source('Pr_occupancy_season.R')
    occ_post_obs <- mat.or.vec(2, length(BOEM[,1]))
    for(z in 1:length(index)){
      lat_input_min <- BOEM$ymin[z]
      lat_input_max <- BOEM$ymax[z]
      lon_input_min <- BOEM$xmin[z]
      lon_input_max <- BOEM$xmax[z]
      occ_post_obs[, z] <- Pr_occupancy_season(posts_latb, posts_lonb, lat_input_min, lat_input_max, lon_input_min, lon_input_max, side, n.iter, N, num_day, "time_step")
    }
    
    occ_post_plot <- occ_post[, order(min_motus)]
    occ_post_pred_plot <- occ_post_pred[, order(min_motus)]
    occ_post_obs_plot <- occ_post_obs[, order(min_motus)]
    
    for(i in 1:length(min_motus)){
      a1 <- occ_post_obs_plot[1, i] - quantile(occ_post_pred_plot[,i], c(0.025))
      b1 <- quantile(occ_post_pred_plot[,i], c(0.975)) - occ_post_obs_plot[1,i] 
      out_bounds[m, i] <- length(which(c(a1, b1) < 0))
    }
    
    if(mute3[p] != 9){
    params[1, 1, m] <- mean(MOVE$BUGSoutput$sims.array[, , "gamma"])
    params[1, 2:3, m] <- quantile(MOVE$BUGSoutput$sims.array[, , "gamma"], c(.025, 0.975))
    params[2, 1, m] <- mean(MOVE$BUGSoutput$sims.array[, , "D[1]"])
    params[2, 2:3, m] <- quantile(MOVE$BUGSoutput$sims.array[, , "D[1]"], c(.025, 0.975))
    params[3, 1, m] <- mean(MOVE$BUGSoutput$sims.array[, , "D[2]"])
    params[3, 2:3, m] <- quantile(MOVE$BUGSoutput$sims.array[, , "D[2]"], c(.025, 0.975))
    params[4, 1, m] <- mean(MOVE$BUGSoutput$sims.array[, , "Sigma[1,1]"])
    params[4, 2:3, m] <- quantile(MOVE$BUGSoutput$sims.array[, , "Sigma[1,1]"], c(.025, 0.975))
    params[5, 1, m] <- mean(MOVE$BUGSoutput$sims.array[, , "Sigma[2,1]"])
    params[5, 2:3, m] <- quantile(MOVE$BUGSoutput$sims.array[, , "Sigma[2,1]"], c(.025, 0.975))
    params[6, 1, m] <- mean(MOVE$BUGSoutput$sims.array[, , "Sigma[1,2]"])
    params[6, 2:3, m] <- quantile(MOVE$BUGSoutput$sims.array[, , "Sigma[1,2]"], c(.025, 0.975))
    params[7, 1, m] <- mean(MOVE$BUGSoutput$sims.array[, , "Sigma[2,2]"])
    params[7, 2:3, m] <- quantile(MOVE$BUGSoutput$sims.array[, , "Sigma[2,2]"], c(.025, 0.975))
    }else{
      params[1, 1, m] <- mean(MOVE$BUGSoutput$sims.array[, , "gamma[1]"])
      params[1, 2:3, m] <- quantile(MOVE$BUGSoutput$sims.array[, , "gamma[1]"], c(.025, 0.975))
      params[2, 1, m] <- mean(MOVE$BUGSoutput$sims.array[, , "D[1,1]"])
      params[2, 2:3, m] <- quantile(MOVE$BUGSoutput$sims.array[, , "D[1,1]"], c(.025, 0.975))
      params[3, 1, m] <- mean(MOVE$BUGSoutput$sims.array[, , "D[2,1]"])
      params[3, 2:3, m] <- quantile(MOVE$BUGSoutput$sims.array[, , "D[2,1]"], c(.025, 0.975))
      params[4, 1, m] <- mean(MOVE$BUGSoutput$sims.array[, , "Sigma[1,1,1]"])
      params[4, 2:3, m] <- quantile(MOVE$BUGSoutput$sims.array[, , "Sigma[1,1,1]"], c(.025, 0.975))
      params[5, 1, m] <- mean(MOVE$BUGSoutput$sims.array[, , "Sigma[2,1,1]"])
      params[5, 2:3, m] <- quantile(MOVE$BUGSoutput$sims.array[, , "Sigma[2,1,1]"], c(.025, 0.975))
      params[6, 1, m] <- mean(MOVE$BUGSoutput$sims.array[, , "Sigma[1,2,1]"])
      params[6, 2:3, m] <- quantile(MOVE$BUGSoutput$sims.array[, , "Sigma[1,2,1]"], c(.025, 0.975))
      params[7, 1, m] <- mean(MOVE$BUGSoutput$sims.array[, , "Sigma[2,2,1]"])
      params[7, 2:3, m] <- quantile(MOVE$BUGSoutput$sims.array[, , "Sigma[2,2,1]"], c(.025, 0.975))
    }
    #j
    out_bounds_array[m, , p] <- out_bounds[m, i]
    post_pred_obs[, m, p] <- occ_post_obs_plot[1,]
    post_pred_mod[, m, p] <- colMeans(occ_post_plot)
    post_pred_mod_check[, m, p] <- colMeans(occ_post_pred_plot)
    post_pred_ecdf[, m, p] <- unlist(lapply(1:length(occ_post_obs_plot[1,]), function(x) ecdf(occ_post_pred_plot[,x])(occ_post_obs_plot[1,x])))
    #rm(MOVE)
  }
  #p
  params_out[[p]] <- params
}


scenario <- c("No missing", "Missing at random", "Motus", "Motus 2-state")

layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), 4, 3, byrow = FALSE))
#layout(matrix(c(1, 2, 3), 3, 1, byrow = TRUE))
par(mar=c(3, 4, 2, 1))
for(p in 4:4){
plot(1, 1, col=rgb(0, 0, 0, 0), xlim=c(0.5, 7.5), ylim=c(-0.5, 0.5), xlab=" ", ylab=" ")
for(i in 1:7){
  segments(i, quantile(params_out[[p]][i, 1, ], c(0.025)), i, quantile(params_out[[p]][i, 1, ], c(0.975)), 
           col=rgb(0, 0, 0, 0.5), lend="butt", lwd=8)
  segments(i, quantile(params_out[[p]][i, 2, ], c(.025)), i, quantile(params_out[[p]][i, 3, ], c(.975)), 
           col=rgb(1, 0, 0, 0.2), lend="butt", lwd=8)
}
if(mute3 == 0){
points(1:7, c(0.2, -0.11, -0.11, 0.01, 0, 0, 0.01), pch=16)
}else{
  points(1:7, c(0.2, -0.11, -0.11, 0.01, 0, 0, 0.01), pch=16) 
}
mtext(side=3, line=0.2, scenario[p], cex=0.6, adj=0, font=2)
mtext(side=1, line=2, "Movement parameters", cex=0.6)
mtext(side=2, line=2, "Parameter value", cex=0.6)
}

# column 2
par(mar=c(3, 4, 2, 1))
for(p in 4:4){
  plot(0, 0, col=rgb(0, 0, 0, 0), ylim=c(0, 1), xlim=c(min(min_motus)*0.95, max(min_motus)*1.05), xlab=" ", ylab=" ")
  ecfd_pts <- apply(post_pred_ecdf, 1, colMeans)[p, ]
  ecfd_pts <- post_pred_ecdf[,1,p]
points(min_motus, ecfd_pts, pch=16, cex=0.7)
mtext(side=2, line=2, "ECDF", cex=0.6)
mtext(side=1, line=2, "Distance from closest receiver", cex=0.6)
}

#column 3
par(mar=c(3, 4, 2, 1))
for(p in 4:4){
plot(apply(post_pred_obs, 1, colMeans)[p,], apply(post_pred_mod_check[,,p], 1, median), xlab=" ", ylab=" ", pch=16, cex=0.6)
mtext(side=2, line=2, "Occupancy using observed", cex=0.6)
mtext(side=1, line=2, "Occupancy using modeled", cex=0.6)
abline(0, 1)
}

for(p in 1:3){
  plot(min_motus, apply(out_bounds_array[,,p], 2, mean), pch=16)
  mtext(side=2, line=2, "Prop. of iterations outside bounds", cex=0.6)
  mtext(side=1, line=2, "Distance from closest reciever", cex=0.6)
}
