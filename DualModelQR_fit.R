library(quantregGrowth)
library(doParallel)
library(foreach)
library(progress)


# --- Functions ---

# Load functions from DualModelQR_functions.R code
source("path") # specify the correct path



# --- Model fitting ---

# Path of the dataset containing the ensemble data and observations
data_path <- "path" # specify the correct path

# Path of results
results_path <- "path" # specify the correct path

# Load the data
windDataSum <- readRDS(data_path)

allDates <- dimnames(windDataSum)[[1]]
nDates <- length(allDates)

allLocs <- dimnames(windDataSum)[[2]]
nLocs <- length(allLocs)

LeadTimes <- dimnames(windDataSum)[[4]]
nLeadTimes <- length(LeadTimes)


### Quantile regression ###

# Training length
tr_len <- 60

# Quantiles to estimate
taus <- seq(from = 1/51, to = 50/51, by = 1/51)
n_taus <- length(taus)

# Forecast dates
fDates <- allDates[-(1:tr_len)]
nfDates <- length(fDates)

# Quantile names
QuantileNames <- mapply(paste0, as.character(1:n_taus), rep("/51Q", n_taus), USE.NAMES = FALSE)

# Number of clusters
n_clust <- 100

# Number of features for clustering
n_feat <- 24

# Set up parallel computation
cores <- detectCores()
c1 <- makeCluster(cores - 4) # not to overload the computer
registerDoParallel(c1)

# Modeling for each lead time
for(lt in 1:nLeadTimes){
  start_time <- Sys.time()
  print(paste0("Lead time: ", LeadTimes[lt], "h"))
  print(start_time)
  
  # Data of current lead time
  windDataSumlt <- windDataSum[,,,lt]
  
  # Parallel modeling for each forecast date taking the lead time into account
  # Output: Three dimensional array with quantile regression predictions for each forecast date and station at the current lead time
  QRS <- foreach(d = (tr_len + lt):nDates, .combine = 'acomb3', .multicombine = TRUE, .packages = c("quantregGrowth")) %dopar% {
    
    # Set seed for reproducibility
    set.seed(d * 104729)
    
    # Training dates
    train_idx <- (d - tr_len - lt + 1):(d - lt)
    
    # Matrix for storing results
    res <- matrix(NA, nrow = nLocs, ncol = n_taus + 3)
    rownames(res) <- allLocs
    colnames(res) <- c(QuantileNames, "obs", "cluster", "model")
    
    # Cluster analysis of the stations
    clusters <- cluster.fit(obs.train = t(windDataSumlt[train_idx,,"obs"]), ens.train.mean = t(windDataSumlt[train_idx,,"EnsMean"]),
                            nStat = nLocs, no.feat = n_feat, no.clus = n_clust, method = "both")
    
    # Training observations for all stations
    obs_all <- windDataSumlt[train_idx,, "obs"]
    
    # 0.1-quantiles of wind speed ensemble forecasts for all stations
    EnsQ10_all <- windDataSumlt[train_idx,, "EnsQ10"]
    
    # 0.5-quantiles of wind speed ensemble forecasts for all stations
    EnsQ50_all <- windDataSumlt[train_idx,, "EnsQ50"]
    
    # 0.9-quantiles of wind speed ensemble forecasts for all stations
    EnsQ90_all <- windDataSumlt[train_idx,, "EnsQ90"]
    
    # Data for prediction for all stations
    nd_all <- as.data.frame(windDataSumlt[d,, c("EnsQ10", "EnsQ50", "EnsQ90")])
    
    # Estimating regression models for each cluster of stations
    for(cl in 1:n_clust){
      
      # Indicies of stations in the current cluster
      current_sts <- which(clusters == cl)
      
      # Observation at stations in the current cluster
      obs <- c(obs_all[,current_sts])
      
      # 0.1-quantile of ensembles at stations in the current cluster
      EnsQ10 <- c(EnsQ10_all[,current_sts])
      
      # 0.5-quantile of ensembles at stations in the current cluster
      EnsQ50 <- c(EnsQ50_all[,current_sts])
      
      # 0.9-quantile of ensembles at stations in the current cluster
      EnsQ90 <- c(EnsQ90_all[,current_sts])
      
      # Results of model fitting
      qr_results <- QRfit(EnsQ10 = EnsQ10, EnsQ50 = EnsQ50, EnsQ90 = EnsQ90, obs = obs, taus = taus)
      
      model <- qr_results$model
      qmod <- qr_results$qmod
      
      # Matrix to store predictions
      pred_matrix <- matrix(NA, nrow = length(current_sts), ncol = n_taus)
      rownames(pred_matrix) <- allLocs[current_sts]
      colnames(pred_matrix) <- QuantileNames
      
      # Prediction for the different model settings
      if(model %in% c(1, 2, 4)){
        
        # New data for prediction for stations in the current cluster
        nd <- nd_all[current_sts,]
        
        # Complete cases
        cc <- complete.cases(nd)
        
        # Handle the case when the current cluster has only one element
        if(sum(cc) > 1){
          pred_matrix[cc,] <- predict(qmod, newdata = nd[cc,])
          
        } else if(sum(cc) == 1){
          nd <- nd[cc,]
          nd <- rbind(nd, c(0, 0, 0))
          prediction <- predict(qmod, newdata = nd)
          pred_matrix[which(cc),] <- matrix(prediction[1,], nrow = 1, byrow = TRUE)
        }
        
      } else if(model == 3){
        
        # New data for prediction for stations in the current cluster
        nd <- nd_all[current_sts, "EnsQ50"]
        
        # Complete cases
        cc <- complete.cases(nd)
        
        # The case when the current cluster has only one element needs to be handled separately
        if(sum(cc) > 1){
          nd <- nd[cc]
          nd <- data.frame(EnsQ50 = nd)
          pred_matrix[cc,] <- predict(qmod, newdata = nd)
          
        } else if(sum(cc) == 1){
          nd <- nd[cc]
          nd <- data.frame(EnsQ50 = nd)
          nd <- rbind(nd, 0)
          prediction <- predict(qmod, newdata = nd)
          pred_matrix[which(cc),] <- matrix(prediction[1,], nrow = 1, byrow = TRUE)
        }
      }
      
      # Inverse transform of the logarithmic transformation applied in linear QR model
      # For the other cases this is handled automatically during prediction
      if(model == 4 && sum(cc) > 0){
        pred_matrix <- exp(pred_matrix) - 0.01
      }
      
      # Replace negative predictions
      if(model != 5){
        
        # Cases having negative predictions
        bad_cases <- which(rowSums(pred_matrix < 0, na.rm = T) > 0)
        
        # The case when there is exactly one such case needs to be handled separately
        if(length(bad_cases) > 1){
          pred_matrix[bad_cases,] <- t(apply(pred_matrix[bad_cases,], 1, replace.negative.quantiles)) 
        } else if(length(bad_cases) == 1){
          pred_matrix[bad_cases,] <- replace.negative.quantiles(pred_matrix[bad_cases,])
        }
      }
      
      # Store results
      res[current_sts, 1:50] <- pred_matrix
      res[current_sts, "cluster"] <- cl
      res[current_sts, "model"] <- model
      
      # Remove prediction matrix
      rm(pred_matrix)
      
      # Garbage collection
      gc()
    }
    
    res
  }
  
  # Permute the dimensions of the array of results
  QRS <- aperm(QRS, perm = c(3, 1, 2))
  
  # Set the names of the forecast dates
  dimnames(QRS)[[1]] <- fDates[lt:93]
  
  # Load observations into the array of results
  QRS[,,"obs"] <- windDataSumlt[fDates[lt:93],,"obs"]
  
  # Save results
  saveRDS(QRS, file = paste0(results_path, "windQR", LeadTimes[lt], ".rds"))
  
  # Garbage collection
  gc()
  
  # Elapsed time for current lead time
  end_time <- Sys.time()
  print(end_time - start_time)
}

# Stopping the cluster for parallel computation
stopCluster(c1)

# Garbage collection
gc() 
  