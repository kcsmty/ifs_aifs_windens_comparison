# Function to combine arrays by their 3rd dimensions
acomb3 <- function(...) abind::abind(..., along = 3)

# Function performing cluster analysis of stations
cluster.fit <- function(obs.train, ens.train.mean, nStat, no.feat, no.clus, method){
  
  # Features to use for clustering
  # ------------------------------
  # Options:
  #
  # climat: Features are specific quantiles of observations
  # error: Features are specific quantiles of forecast errors
  # forc: Features are specific quantiles of median forecasts
  # both: Features are specific quantiles of observations and forecast errors
  switch (method,
          climat = {
            cluster.data <- apply(obs.train, 1, quantile, probs=c(1:no.feat)/(no.feat+1), na.rm = T)
          },
          error = {
            cluster.data <- apply((ens.train.mean - obs.train), 1, quantile, probs = c(1:no.feat)/(no.feat+1), na.rm = T)
          },
          forc = {
            cluster.data <- apply(ens.train.mean, 1, quantile, probs = c(1:no.feat)/(no.feat+1), na.rm = T)
          },
          both = {
            n1 <- floor(no.feat/2)
            n2 <- no.feat - n1
            cluster.data <- 
              rbind(apply(obs.train, 1, quantile, probs = c(1:n1)/(n1+1), na.rm = T),
                    apply((ens.train.mean - obs.train), 1, quantile, probs = c(1:n2)/(n2+1), na.rm = T))
          }
  )
  
  # Complete cases
  good.stat <- !apply(cluster.data, 2, anyNA)
  
  # Vector for results
  cluster.id <- rep(NA, nStat)
  
  # Cluster analysis of stations using k-means algorithm
  cluster.id[good.stat] <- kmeans(t(cluster.data[,good.stat]), no.clus, iter.max = 20)$cluster
  
  return(cluster.id)
}

# Function to handle negative predicted wind speed quantiles
# in the input vector quants
replace.negative.quantiles <- function(quants){
  
  # Number of negative predictions
  num_neg <- length(which(quants < 0))
  
  if(num_neg > 0){
    
    # First nonnegative prediction 
    pos <- quants[which(quants >= 0)][1]
    
    # Replace negative predictions with equidistant values between 0 and pos
    quants[1:num_neg] <- (1:num_neg)/(num_neg + 1)*pos
    
    return(quants)
    
  } else {
    
    return(quants)
    
  }
}

# Function for fitting QR models
QRfit <- function(EnsQ10, EnsQ50, EnsQ90, obs, taus){
  
  # Model codes:
  #
  # 1: primary nonlinear model
  # 2: nonlinear model with increased smoothing parameter and decreased number of nodes
  # 3: nonlinear univariate model (independent variable: ensemble median)
  # 4: linear quantile regression model
  # 5: unsuccessful model fitting
  model <- 1
  qmod <- try(gcrq(obs ~ ps(EnsQ10, lambda = 1) +
                     ps(EnsQ50, lambda = 1) + 
                     ps(EnsQ90, lambda = 1),
                   tau = taus, transf = "log(y + 0.01)"),
              silent = TRUE)
  
  if(inherits(qmod, "try-error")){
    model <- 2
    qmod <- try(gcrq(obs ~ ps(EnsQ10, lambda = 1.2, ndx = 5) +
                       ps(EnsQ50, lambda = 1.2, ndx = 5) +
                       ps(EnsQ90, lambda = 1.2, ndx = 5),
                     tau = taus, transf = "log(y + 0.01)"),
                silent = TRUE)
    
    if(inherits(qmod, "try-error")){
      model <- 3
      qmod <- try(gcrq(obs ~ ps(EnsQ50, lambda =  1),
                       tau = taus, transf = "log(y + 0.01)"),
                  silent = TRUE)
      
      if(inherits(qmod, "try-error")){
        model <- 4
        qmod <- try(rq(log(obs + 0.01) ~ EnsQ10 + EnsQ50 + EnsQ90,
                       tau = taus, method = 'fn'))
        
        if(inherits(qmod, "try-error")){
          model <- 5
        }
      }
    }
  }
  
  return(list(model = model,
              qmod = qmod))
}
