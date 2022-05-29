This repository contains scripts for implementing a general modeling framework for inferring movements automated telemetry data, like Motus or GPS tracking. The approach is a state-space implementation of a correlated random walk model, based on the modeling framework from Baldwin et al. (2018) and Jonsen (2016) – i.e. the model is built around a multivariate normal distribution that describes movements in x and y dimensions. The main difference between this implementation of the modeling framework and previous ones is that I use a more stable implementation of the prior for the variance-covariance matrix. Like other versions, my implementation allow several extensions that make it possible to relax assumptions that might be too restrictive, including multiple states (e.g. staging vs. migration), correlation in the distance travelled between successive movements, tendency to drift in a particular direction over time (e.g. an individual during migratory staging might make random movements that can be in any direction but have a greater tendency to be in a south-western direction), physiological or behavioral constraints that specifying the maximum distance an individual can travel per day, and allowing missing data. 
The scripts here prepare Motus data files for analyses, run a movement model in JAGS, use posterior estimates to obtain occupancy estimates for user-specified spatial units, conducted posterior predictive checks using several test statistics, create simple plots showing mean and variance of results over space, and simulate the potential bias from estimating movements in offshore areas from land-based receivers. These models were developed specifically to model movements of species listed on the Endangered Species Act across the Atlantic Outer Continental Shelf. 

'false_pos_filter.R' - add a column for burst length to a dataframe 

'remove_false_pos.R' - filter false positives using the same rules as the Motus project

'remove_dup_bursts.R' - remove duplicate detections by defining bursts as all detections within a 24-hour period

'serial_time.R' – add serial time to a dataframe

'MOTUS_JAGS.R' – JAGS pseudo-code for general movement modeling framework for 1, 2, and 3 state specifications. This version uses a more stable implementation of the prior for the variance-covariance matrix.

'Pr_occupancy_season.R' – run one of 4 options for propagating stochasticity and uncertainty of posterior estimates for user-defined spatial units. Options vary by run time, with propagation of all sources taking the longest. Option 1: full propagation of uncertainty, using JAGS (this version of extrapolation takes a very long time to run). Option 2: extrapolation only includes uncertainty from estimating the proportion of individuals who likely crossed into the specified extent (gives estimates by month instead of year). Option 2b: extrapolation only includes uncertainty from estimating the proportion of individuals who likely crossed into the specified extent (gives estimates by year). Option 3: extrapolation only includes uncertainty from estimating the proportion of individuals who likely crossed into the specified extent (gives estimates by year). Option 4: extrapolation starts with posterior means to get rapid point estimates

‘MOTUSmove_walkthrough.Rmd’ – R Markdown file that demonstrates basic principles of correlated random walk models and how model complexity is added

‘motus_array_sim.R’ – simulation for quantifying the bias from estimating offshore movements from land-based receivers 
