library(scoringRules)
library(doParallel)
library(foreach)
library(abind)


# ------ Data paths ------

# Path of the dataset containing the ensemble data and observations
data_path <- "path" # specify the correct path

# Path of results
results_path <- "path" # specify the correct path

# ------ Functions ------

# Load function from DualModelEmos_function.R file
source("path") # specify the correct path



# ------ Define variables ------

# Load the forecast-observation data
windData <- readRDS(data_path)

# Dates
allDates <- dimnames(windData)[[1]]
nDates <- length(allDates)

# Stations
allLocs <- dimnames(windData)[[2]]
nLocs <- length(allLocs)

# Lead times
LeadTimes <- dimnames(windData)[[4]]
nLeadTimes <- length(LeadTimes)

# Parameters of the EMOS models
# -----------------------------
# a and b: constant term and coefficient of the ensemble mean in the location of the fitted TN distribution
# c and d: constant term and coefficient of the ensemble variance in the scale of the fitted TN distribution
# conv: variable to record convergence status of the models
# reg: variable to record if regularization is applied
Parameters <- c("a", "b", "c", "d", "conv", "reg")
nParameters <- length(Parameters)

# Length of the rolling training period
tr_len <- 60

# Forecast dates
fDates <- allDates[-(1:tr_len)]
nfDates <- length(fDates)

# Array for results
windEmos <- array(NA, dim = c(nfDates, nLocs, nParameters, nLeadTimes), 
                  dimnames = list(fDates, allLocs, Parameters, LeadTimes))

# Regularization parameter
lambda <- 0.001


# ------ Fit local EMOS models ------

# Set up clusters for parallel modeling
ncores <- detectCores() - 4 # not to overload the computer
c1 <- makeCluster(ncores)
registerDoParallel(c1)

# Modeling for each lead time
for(lt in 1:nLeadTimes){
  
  # Starting time for current lead time
  start_time <- Sys.time()
  
  # Print which lead time is being processed
  print(paste0("Lead time: ", LeadTimes[lt], "h"))
  
  # Print starting time
  print(start_time)
  
  # Data for current lead time
  windDatalt <- windData[,,,lt]
  
  # Parallel modeling for each forecast date taking the lead time into account
  # Output: Three dimensional array with estimated EMOS parameters for each forecast date and station at the current lead time
  EmosPars <- foreach(d = (tr_len + lt):nDates, .combine = 'acomb3', .multicombine = TRUE, .packages = c("scoringRules")) %dopar% {
    
    # Training cases for current forecast date
    train_idx <- (d - tr_len - lt + 1):(d - lt)
    
    # Matrix to store results
    res <- matrix(NA, nrow = nLocs, ncol = 6)
    rownames(res) <- allLocs
    colnames(res) <- c("a", "b", "c", "d", "conv", "reg") 
    
    # Estimate EMOS parameters for each station
    for(i in 1:nLocs){
      
      # Training data for current station
      train_data <- windDatalt[train_idx, i,]
      train_data <- train_data[complete.cases(train_data),]
      train_em <- train_data[,"EnsMean"]
      train_ev <- train_data[,"EnsVar"]
      train_obs <- train_data[,"obs"]
      rm(train_data)
      
      # Store results
      res[i,] <- EMOSfit(train_obs = train_obs, train_em = train_em, train_ev = train_ev, lambda = lambda)
    }
    
    res
  }
  
  # Permute the dimensions of the array of results
  EmosPars <- aperm(EmosPars, perm = c(3, 1, 2))
  
  # Store results of current lead time
  windEmos[lt:nfDates,,,lt] <- EmosPars
  
  # Garbage collection
  gc()
  
  # Print elapsed time
  end_time <- Sys.time()
  print(end_time - start_time)
}

# Stopping the cluster for parallel computation
stopCluster(c1)

# Garbage collection
gc()

# Save results
saveRDS(windEmos, file = paste0(results_path, "windEmosParsLOC.rds"))