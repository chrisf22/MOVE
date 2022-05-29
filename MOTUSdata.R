##################################################################################################################################
### Script to 1) prepare Motus data files for analyses, 2) run movement model in JAGS, 
### 3) use posterior estimates to obtain occupancy estimates for user-specified spatial units, and 
### 4) create simple plots showing mean and variance of results over space. 

### Depending on which components are run, requires: 'false_pos_filter.R', 'remove_false_pos.R', 
### 'remove_dup_bursts.R', 'serial_time.R', 'MOTUS_JAGS.R', 'Pr_occupancy_season.R'
##################################################################################################################################


## load required packages
library(R2jags)
library(R2WinBUGS)
library(sp)
library(raster)
library(rgeos)
library(sf)


## read GPS data for Red Knot for estimating flight heights ######################################################################
rekn_GPS <- read.csv("~/Documents/folders/URI/MOTUSdata/GPS/GPSTrackingData.csv", header=TRUE, stringsAsFactors = FALSE)
rekn_id <- unique(rekn_GPS$Red.Knot.individual.Identity)
for(i in 1:length(rekn_id)){
  rekn_GPS[rekn_GPS$Red.Knot.individual.Identity == rekn_id[i], 6]
}


## read and combine Motus data for Red Knot, Roseate Tern, Common Tern, Piping Plover for movement analyses ######################
# rekn
proj_14_2016_rekn_final <- readRDS("~/Documents/folders/URI/MOTUSData/detection_data/proj_14_2016_rekn_final.rds")
proj_15_2016_rekn_final <- readRDS("~/Documents/folders/URI/MOTUSData/detection_data/proj_15_2016_rekn_final.rds")
proj_38_2016_rekn_final <- readRDS("~/Documents/folders/URI/MOTUSData/detection_data/proj_38_2016_rekn_final.rds")
#proj_47_2016_rekn_final <- readRDS("~/Documents/folders/URI/MOTUSData/detection_data/proj_47_2016_rekn_spring_final.rds")
proj_88_2016_rekn_final <- readRDS("~/Documents/folders/URI/MOTUSData/detection_data/proj_88_2016_rekn_final.rds")
allindvs <- rbind(proj_14_2016_rekn_final, proj_15_2016_rekn_final, proj_38_2016_rekn_final, proj_88_2016_rekn_final)
# remove the one detection from 2017 so that the filters do not think the earliest detection was in January
allindvs <- allindvs[-which(as.numeric(substr(allindvs$ts_gmt, 1, 4))==2017), ]
#allindvs <- proj_14_2016_rekn_final
allindvs <- allindvs[!is.na(allindvs$lat)&!is.na(allindvs$lon), ]

# rote
proj_14_2015_rote_final <- read.csv("~/Documents/folders/URI/MOTUSData/Appendix_H-Motus_Detection_Data/Appendix_H-Motus_Detection_Data_ROST_2015_v20190708.csv", header=TRUE, stringsAsFactors = FALSE)
proj_14_2016_rote_final <- read.csv("~/Documents/folders/URI/MOTUSData/Appendix_H-Motus_Detection_Data/Appendix_H-Motus_Detection_Data_ROST_2016_v20190708.csv", header=TRUE, stringsAsFactors = FALSE)
proj_14_2017_rote_final <- read.csv("~/Documents/folders/URI/MOTUSData/Appendix_H-Motus_Detection_Data/Appendix_H-Motus_Detection_Data_ROST_2017_v20190708.csv", header=TRUE, stringsAsFactors = FALSE)
allindvs <- rbind(proj_14_2015_rote_final, proj_14_2016_rote_final, proj_14_2017_rote_final)
allindvs <- allindvs[!is.na(allindvs$lat)&!is.na(allindvs$lon), ]

# pipl
proj_14_2015_pipl_final <- read.csv("~/Documents/folders/URI/MOTUSData/Appendix_H-Motus_Detection_Data/Appendix_H-Motus_Detection_Data_PIPL_2015.csv", header=TRUE, stringsAsFactors = FALSE)
proj_14_2016_pipl_final <- read.csv("~/Documents/folders/URI/MOTUSData/Appendix_H-Motus_Detection_Data/Appendix_H-Motus_Detection_Data_PIPL_2016.csv", header=TRUE, stringsAsFactors = FALSE)
proj_14_2017_pipl_final <- read.csv("~/Documents/folders/URI/MOTUSData/Appendix_H-Motus_Detection_Data/Appendix_H-Motus_Detection_Data_PIPL_2017.csv", header=TRUE, stringsAsFactors = FALSE)
allindvs <- rbind(proj_14_2015_pipl_final, proj_14_2016_pipl_final, proj_14_2017_pipl_final)
allindvs <- allindvs[!is.na(allindvs$lat)&!is.na(allindvs$lon), ]

# cote
proj_14_2015_cote_final <- read.csv("~/Documents/folders/URI/MOTUSData/Appendix_H-Motus_Detection_Data/Appendix_H-Motus_Detection_Data_COTE_2015_v20190708.csv", header=TRUE, stringsAsFactors = FALSE)
proj_14_2016_cote_final <- read.csv("~/Documents/folders/URI/MOTUSData/Appendix_H-Motus_Detection_Data/Appendix_H-Motus_Detection_Data_COTE_2016_v20190708.csv", header=TRUE, stringsAsFactors = FALSE)
proj_14_2017_cote_final <- read.csv("~/Documents/folders/URI/MOTUSData/Appendix_H-Motus_Detection_Data/Appendix_H-Motus_Detection_Data_COTE_2017_v20190708.csv", header=TRUE, stringsAsFactors = FALSE)
allindvs <- rbind(proj_14_2015_cote_final, proj_14_2016_cote_final, proj_14_2017_cote_final)
allindvs <- allindvs[!is.na(allindvs$lat)&!is.na(allindvs$lon), ]



## create columns for time variables for movement model ##########################################################################
# create a new column for hour of detection with minutes represented by a decimal
hours_dec <- round(as.numeric(substr(allindvs$ts_gmt, 12, 13)) + 
                     (as.numeric(substr(allindvs$ts_gmt, 15, 16))/60), 2)
allindvs <- cbind(allindvs, hours_dec)

# create a new column for month
month <- as.numeric(substr(allindvs$ts_gmt, 6, 7))
allindvs <- cbind(allindvs, month)

# create a new column for day
day <- as.numeric(substr(allindvs$ts_gmt, 9, 10))
allindvs <- cbind(allindvs, day)

# create a new column for year
year <- as.numeric(substr(allindvs$ts_gmt, 1, 4))
allindvs<- cbind(allindvs, year)

# order the detections by indiviudal, year, day, hour, and minute
allindvs_ordered <- allindvs[order(allindvs$id, allindvs$year, substr(allindvs$ts_gmt, 6, 11)), ]



## pass data through filters for false positives, duplicate bursts, and stationary individuals ###################################
# one option for filtering motus data is 'filterByActivity.R' from the 'motus' package
# filterByActivity() from 'motus' requires the data in an SQL table
# using filtering scripts that run on flat data frames instead until Motus releases pre-filtered downloadable data

# add column for burst length
setwd("~/Documents/folders/GitProjects/URI_CRM/")
source('false_pos_filter.R')
#burst_length <- unlist(lapply(1:length(allindvs_ordered[,1]), FUN=false_pos_filter, input_data = allindvs_ordered))
burst_length <- false_pos_filter(allindvs_ordered)
allindvs_ordered <- cbind(allindvs_ordered, burst_length)

# remove false positives 
source('remove_false_pos.R') 
allindvs_ordered_filtered <- remove_false_pos(allindvs_ordered, 3)

# remove duplicate detections by defining bursts as all detections within a 24-hour period
source('remove_dup_bursts.R')
allindvs_ordered_filtered_nodups <- remove_dup_bursts(allindvs_ordered_filtered, "day")

#source('serial_time.R')
#allindvs_ordered_filtered_nodups <- serial_time(allindvs_ordered_filtered_nodups, 'day')

# remove any individuals who do not leave the location of their initial detection
keep <- NA
mig <- mat.or.vec(length(unique(allindvs_ordered_filtered_nodups$id)), 1)
allindvs_ordered_filtered_nodups_mig <- allindvs_ordered_filtered_nodups
out <- mat.or.vec(length(unique(allindvs_ordered_filtered_nodups$id)), 1)
for(i in 1:length(unique(allindvs_ordered_filtered_nodups$id))){
  temp <- allindvs_ordered_filtered_nodups[allindvs_ordered_filtered_nodups$id == unique(allindvs_ordered_filtered_nodups$id)[i], ]
  temp_min <- 1
  temp_max <- length(temp$lat)
  out[i] <- length(temp$lat)
  temp_dist <- sqrt((round(temp$lat[temp_max], 4) - round(temp$lat[temp_min], 4))^2 + (round(temp$lon[temp_max], 4) - round(temp$lon[temp_min], 4))^2)
  if(temp_dist > 0){
    mig[i] <- 1
  }
  if(mig[i] == 1){
    keep <- c(keep, which(allindvs_ordered_filtered_nodups_mig$id == unique(allindvs_ordered_filtered_nodups$id)[i]))
  }
}
keep <- keep[2:length(keep)]
allindvs_ordered_filtered_nodups <- allindvs_ordered_filtered_nodups_mig[keep, ]



## create variables for indexing the JAGS movement model (sensu Baldwin et al. 2018) #############################################
# create a vector (Sind_obs) that references, for each individual, the position of the vector (for the observed data) that denotes its first detection
id_index <- mat.or.vec(length(allindvs_ordered_filtered_nodups[,1]), 1)
allindvs_ordered_filtered_nodups <- cbind(allindvs_ordered_filtered_nodups, id_index)
# create an index for individual
uni_inds <- unique(allindvs_ordered_filtered_nodups$id)
Sind_obs <- unique(allindvs_ordered_filtered_nodups$id)
for(i in 1:length(uni_inds)){
  allindvs_ordered_filtered_nodups[allindvs_ordered_filtered_nodups$id==uni_inds[i], 'id_index'] <- i
  Sind_obs[i] <- min(which(allindvs_ordered_filtered_nodups$id==uni_inds[i]))
}

source('serial_time.R')
allindvs_ordered_filtered_nodups <- serial_time(allindvs_ordered_filtered_nodups, 'day')
season_length <- max(allindvs_ordered_filtered_nodups$days_since)

# number of individuals
N <- length(uni_inds)
#N <- 30

# create a vector (Xidx) that references, for each individual, the position of the vector (for the regular time steps) that denotes its first detection
Xidx <- mat.or.vec(N, 1)
# Xidx2 indexes the last day of the season, with respect to x, for each individual
Xidx2 <- mat.or.vec(N, 1)
Xidx3 <- mat.or.vec(N, 1)
Zidx <- mat.or.vec(N, 1)
for(i in 1:N){
  Xidx[i] <- min(allindvs_ordered_filtered_nodups[allindvs_ordered_filtered_nodups$id==uni_inds[i], 'days_since']) + season_length*(i-1)
  Xidx2[i] <- season_length*i
  Xidx3[i] <- max(allindvs_ordered_filtered_nodups[allindvs_ordered_filtered_nodups$id==uni_inds[i], 'days_since'])
  Zidx[i] <- min(allindvs_ordered_filtered_nodups[allindvs_ordered_filtered_nodups$id==uni_inds[i], 'days_since'])
}

# create the idx file, which indexes for each observation, which "regular" day it is related to, 
# taking into account the fact that the data for all individuals is specified as one vector
a <- 1
for(i in 1:N){
  # season length * number of individuals already recorded is added to the days since the start of the season
  a <- c(a, (allindvs_ordered_filtered_nodups[allindvs_ordered_filtered_nodups$id==uni_inds[i], 'days_since']+((season_length)*(i-1))))
}
idx <- a[-1]

# locations
y <- allindvs_ordered_filtered_nodups[, c("lat", "lon")]
# there is an additional position added to the end of the vector to ensure that the index loops stay in bounds (as in Baldwin et al. 2018)
Xidx <- c(Xidx, (season_length*N+1))
Xidx2 <- c(Xidx2, (season_length*N+1))
# there is an additional position added to the end of the vector to ensure that the index loops stay in bounds (as in Baldwin et al. 2018)
Yidx <- c(Sind_obs, (length(y[,1])+1))[1:(N+1)]
Y <- season_length

# priors for x and b, specifying NAs for nodes that are never observed
x1 <- rep(NA, season_length*N)
x2 <- rep(NA, season_length*N)
x <- cbind(x1, x2)
b <- rep(NA, season_length*N)
for(k in 1:N){
x[(Xidx[k]+1):Xidx2[k], 1] <- mean(y[(Yidx[k]+1):(Yidx[k+1]-1),1])
x[(Xidx[k]+1):Xidx2[k], 2] <- mean(y[(Yidx[k]+1):(Yidx[k+1]-1),2])
b[(Xidx[k]+1):(Xidx2[k]-1)] <- 1
}

lat_min <- min(y[,1]) - sd(y[,1])/3
lat_max <- max(y[,1]) + sd(y[,1])/3
lon_min <- min(y[,2]) - sd(y[,2])/3
lon_max <- max(y[,2]) + sd(y[,2])/3

# number of iterations for the JAGS model
n.iter <- 1000
n.burnin <- 5000



## run JAGS model ################################################################################################################
setwd("~/Documents/folders/GitProjects/URI_CRM/")
# 2 state
source('MOTUS_JAGS_2states.R')
# 1 state 
source('MOTUS_JAGS_1state.R')

# create a vector of the number of individuals tagged by each month, 
# to ensure number of individuals observed is divided by the correct denominator to get occupancy
first_month <- mat.or.vec(12, length(uni_inds))
for(i in 1:length(uni_inds)){
  # season length * number of individuals already recorded is added to the days since the start of the season
  first_month[min(allindvs_ordered_filtered_nodups[allindvs_ordered_filtered_nodups$id==uni_inds[i], 'month']):12, i] <- 1
}
N_bymonth <- rowSums(first_month)

## use posterior estimates from JAGS model to obtain estimates for user-defined spatial units ####################################
# convert output from JAGS to create arrays (one for lat; one for lon) that have the posterior samples (rows) 
# for each time step of the movement model (columns), and for each individual (depth)
posts_lat <- array(NA, c(n.iter, season_length, N))
posts_lon <- array(NA, c(n.iter, season_length, N))
for(e in 1:N){
  for(i in Zidx[e]:season_length){
    posts_lat[,i,e] <- MOVE$BUGSoutput$sims.array[,,paste("x[", ((e-1)*season_length + i), "," , 1, "]", sep="")]
    posts_lon[,i,e] <- MOVE$BUGSoutput$sims.array[,,paste("x[", ((e-1)*season_length + i), "," , 2, "]", sep="")]
  }
}

# read a file that has the extent of user-defined spatial units as rows, labeled "xmin", "ymax", "xmax", "ymin"
# can use the field calculator in QGIS to get a data file with the extent of each cell and centroid coordinates
# e.g. x_min($geometry) on polygons and $x on a .shp for centroids
setwd("~/Documents/folders/GitProjects/URI_CRM/motus_data/")
#BOEM <- read.csv('BOEM_XY.csv', header=TRUE)
BOEM <- read.csv('BOEM_halfdeg_grid_latlon_att.csv', header=TRUE)
# label columns
colnames(BOEM)[1] <- "xmin"
colnames(BOEM)[2] <- "ymax"
colnames(BOEM)[3] <- "xmax"
colnames(BOEM)[4] <- "ymin"

# this script makes it possible to loop through every nth value
# useful for running through a spatial input file that has too many rows to run wihtin a reasonable time 
nth <- function(x,n){
  x[x%%n==0]
}
x = 1:length(BOEM[,1])
index <- nth(x,1)

## run one of 4 options for propagating stochasticity and uncertainty of posterior estimates ###########################
# options vary by run time, with propagation of all sources taking the longest
# option 1: full propagation of uncertainty, using JAGS
# this version of extrapolation takes a very long time to run
setwd("~/Documents/folders/GitProjects/URI_CRM/")
source('Pr_occupancy_season.R')
occ_post <- mat.or.vec(n.iter, length(BOEM[,1]))
for(z in 1:10){
  lat_input_min <- BOEM$ymin[z]
  lat_input_max <- BOEM$ymax[z]
  lon_input_min <- BOEM$xmin[z]
  lon_input_max <- BOEM$xmax[z]
  occ_post[,z] <- Pr_occupancy_season(posts_lat, posts_lon, lat_input_min, lat_input_max, lon_input_min, lon_input_max, side, n.iter, N, season_length, "JAGS")
}

# option 2: extrapolation only includes uncertainty from estimating the proportion of 
# individuals who likely crossed into the specified extent; gives estimates by month instead of year
posts_latb <- posts_lat[1:1000, , ]
posts_lonb <- posts_lon[1:1000, , ]
setwd("~/Documents/folders/GitProjects/URI_CRM/")
source('Pr_occupancy_season.R')
occ_post <- array(0, c(1000, 12, length(BOEM[,1])))
for(z in 1:length(index)){
#for(z in 1:1){
  lat_input_min <- BOEM$ymin[z]
  lat_input_max <- BOEM$ymax[z]
  lon_input_min <- BOEM$xmin[z]
  lon_input_max <- BOEM$xmax[z]
  occ_post[, ,z] <- Pr_occupancy_season(posts_latb, posts_lonb, lat_input_min, lat_input_max, lon_input_min, lon_input_max, side, n.iter, N_bymonth, season_length, "time_step_bymo")
}


# option 2b: extrapolation only includes uncertainty from estimating the proportion of 
# individuals who likely crossed into the specified extent; gives estimates by year
posts_latb <- posts_lat[1:1000, , ]
posts_lonb <- posts_lon[1:1000, , ]
setwd("~/Documents/folders/GitProjects/URI_CRM/")
source('Pr_occupancy_season.R')
occ_post <- mat.or.vec(1000, length(BOEM[,1]))
for(z in 1:length(index)){
  #for(z in 1:1){
  lat_input_min <- BOEM$ymin[z]
  lat_input_max <- BOEM$ymax[z]
  lon_input_min <- BOEM$xmin[z]
  lon_input_max <- BOEM$xmax[z]
  occ_post[,z] <- Pr_occupancy_season(posts_latb, posts_lonb, lat_input_min, lat_input_max, lon_input_min, lon_input_max, side, n.iter, N, season_length, "time_step")
}


# option 3: extrapolation only includes uncertainty from estimating the proportion of 
# individuals who likely crossed into the specified extent; gives estimates by year
posts_latb <- posts_lat[1:100, , ]
posts_lonb <- posts_lon[1:100, , ]
setwd("~/Documents/folders/GitProjects/URI_CRM/")
source('Pr_occupancy_season.R')
occ_post <- mat.or.vec(100, length(BOEM[,1]))
for(z in 1:length(index)){
  #for(z in 1:10){
  lat_input_min <- BOEM$ymin[z]
  lat_input_max <- BOEM$ymax[z]
  lon_input_min <- BOEM$xmin[z]
  lon_input_max <- BOEM$xmax[z]
  occ_post[,z] <- Pr_occupancy_season(posts_latb, posts_lonb, lat_input_min, lat_input_max, lon_input_min, lon_input_max, side, n.iter, N, season_length, "time_step")
}


# option 4: extrapolation starts with posterior means to get rapid point estimates
posts_latb <- apply(posts_lat, 3, colMeans)
posts_lonb <- apply(posts_lon, 3, colMeans)
setwd("~/Documents/folders/GitProjects/URI_CRM/")
source('Pr_occupancy_season.R')
occ_post <- mat.or.vec(length(BOEM[,1]), 1)
for(z in 1:length(index)){
  lat_input_min <- BOEM$ymin[z]
  lat_input_max <- BOEM$ymax[z]
  lon_input_min <- BOEM$xmin[z]
  lon_input_max <- BOEM$xmax[z]
  occ_post[z] <- Pr_occupancy_season(posts_latb, posts_lonb, lat_input_min, lat_input_max, lon_input_min, lon_input_max, side, n.iter, N, season_length, "fast")
}

saveRDS(occ_post, file = "~/Documents/folders/GitProjects/URI_CRM/motus_data/rote_occ_post.rds")


## create a 2-panel plot (mean and variance) showing movement results for the specified region and spatial units #######
plot(density(occ_post, adjust=2), main=" ", xlab=" ", lwd=2, xlim=c(0.03, 1))
mtext(side=1, line=2, "Pr(occupancy)")

# get the mean and variance of posteriors of each spatial unit for plotting
# if running the full or timestep version of Pr_occupancy_season (options 1, 2, or 3)
occ_alphas <- colMeans(occ_post)
occ_alphas_sd <- apply(occ_post, 2, sd)
#occ_alphas[occ_alphas==0] <- 0.001
#alphas <- log(occ_alphas)/min(log(occ_alphas))
#alphas <- 1 - alphas
#occ_alphas[occ_alphas==0] <- 0.004
#alphas <- 1 - occ_alphas
alphas <- occ_alphas
alphas_sd <- occ_alphas_sd

# if running the rapid version (option 4)
occ_alphas <- occ_post
occ_alphas_sd <- rep(0, length(occ_post))
alphas <- occ_alphas
alphas_sd <- occ_alphas_sd

# load shapefile of spatial units for plot
#setwd("~/Documents/folders/URI/GIS/ATL_BLKCLIP/")
#BOEM.sf <- st_read("ATL_BLKCLP.shp")
setwd("~/Documents/folders/URI/GIS/")
#BOEM.sf <- st_read("BOEM_1deg_grid_latlon.shp")
BOEM.sf <- st_read("BOEM_halfdeg_grid_latlon_att.shp")
#BOEM.sf.alphas <- cbind(BOEM.sf, alphas, alphas_sd, BOEM$lat, BOEM$lon)
BOEM.sf.alphas <- cbind(BOEM.sf, alphas, alphas_sd, BOEM$y_UTM, BOEM$x_UTM)

setwd("~/Documents/folders/URI/GIS/")
motus_rec <- st_read("receiver-deployments_BOEM_XY.shp")

# 2-panel plot
quartz.options(width=7.48, height=5)
layout(matrix(c(1, 2), 1, 2, byrow = TRUE))
#quartz.options(width=5, height=7)
bound <- getData('GADM', country='USA', level=1)
e <- c(min(BOEM$xmin)-1, max(BOEM$xmax)+1, min(BOEM$ymin)-1, max(BOEM$ymax)+1)
boundcrop <- crop(bound, e)
bound.simplify <- gSimplify(boundcrop, tol=0.01, topologyPreserve=TRUE)

par(mar=c(4, 4, 2, 1))
#(level=0 is country boundary, level=1 is state boundaries, level=2 is county boundaries)
plot(0, 0, xlim=c(e[1], e[2]), ylim=c(e[3], e[4]), xaxt="n", yaxt="n", yaxs="i", xaxs="i", xlab=" ", ylab=" ", xpd=FALSE, main="Roseate Tern                 Pr(Passage)", cex.main=0.75)
plot(BOEM.sf.alphas["alphas"], add=TRUE, breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), border=NA)
plot(bound.simplify, add=TRUE, col="white", border="dark gray", lwd=.9)
#plot(st_geometry(BOEM.sf), add=TRUE)
#plot(BOEM.sf.alphas["alphas"], add=TRUE, breaks = c(0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.025, 0.05, 0.075, .1,.15,.2,.25, 0.5, 1), border=NA)
#segments(-68, 37.1, -68-(5*0.8983112), 37.1, lend="butt", lwd=1.2)
#text(-70, 37.25, "500 km")
#text(-74, 37.2, "N", font=2, cex=1.5)
#points(-74, 37.425, pch=17, cex=1.5)
axis(side=1, line=-3, cex.axis=0.75)
axis(side=2, line=-3, cex.axis=0.75)
mtext(side=1, line=-1, "Longitude", cex=0.75)
mtext(side=2, line=-0.5, "Latitude", cex=0.75)
rect(e[1], e[3], e[2], e[4])

par(mar=c(4, 4, 2, 1))
plot(0, 0, xlim=c(e[1], e[2]), ylim=c(e[3], e[4]), xaxt="n", yaxt="n", yaxs="i", xaxs="i", xlab=" ", ylab=" ", xpd=FALSE, main="Roseate Tern                 SD of Pr(Passage)", cex.main=0.75)
plot(BOEM.sf.alphas["alphas_sd"], add=TRUE, breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, .1), border=NA)
plot(bound.simplify, add=TRUE, col="white", border="dark gray", lwd=.9)
axis(side=1, line=-3, cex.axis=0.75)
axis(side=2, line=-3, cex.axis=0.75)
mtext(side=1, line=-1, "Longitude", cex=0.75)
mtext(side=2, line=-0.5, "Latitude", cex=0.75)
rect(e[1], e[3], e[2], e[4])


# plot that shows thresholded confidence in occupancy estimates
quartz.options(width=7.48/2, height=5)
par(mar=c(4, 4, 2, 1))
plot(0, 0, xlim=c(e[1], e[2]), ylim=c(e[3], e[4]), xaxt="n", yaxt="n", yaxs="i", xaxs="i", xlab=" ", ylab=" ", xpd=FALSE, main="Red Knot                 Confidence areas", cex.main=0.75)
plot(BOEM.sf.alphas["alphas_sd"], add=TRUE, breaks = c(0, 0.15, 1), border="black", lwd=0.5, pal=c("light gray", "white"))
plot(BOEM.sf.alphas["alphas"], add=TRUE, breaks = c(0, .5, 1), border="black", lwd=0.5, pal=c(rgb(0, 0, 0, 0), rgb(1, 0, 0, 0.4)))
plot(bound.simplify, add=TRUE, col="white", border="dark gray", lwd=.9)
axis(side=1, line=-3, cex.axis=0.75)
axis(side=2, line=-3, cex.axis=0.75)
mtext(side=1, line=-1, "Longitude", cex=0.75)
mtext(side=2, line=-0.5, "Latitude", cex=0.75)
