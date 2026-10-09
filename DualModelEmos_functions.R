# Function to combine arrays by their 3rd dimensions
acomb3 <- function(...) abind(..., along=3)

# CRPS loss function for EMOS model fitting
CRPS_TN <- function(p, train_obs, train_em, train_ev){
  crps <- sum(crps_tnorm(y = train_obs, location = p[1] + p[2]^2*train_em, scale = sqrt(p[3]^2 + p[4]^2*train_ev), lower = 0), na.rm = TRUE)
  return(crps)
}

# Regularized CRPS loss function for EMOS model fitting
CRPS_TN_REG <- function(p, train_obs, train_em, train_ev, lambda){
  crps <- sum(crps_tnorm(y = train_obs, location = p[1] + p[2]^2*train_em, scale = sqrt(p[3]^2 + p[4]^2*train_ev), lower = 0), na.rm = TRUE)
  reg <- sum(p^2)
  val <- crps + lambda*reg
  return(val)
}

# EMOS model fitting function
EMOSfit <- function(train_obs, train_em, train_ev, lambda){
  
  # Variable to record regularized cases
  reg <- 0
  
  # Estimate initial parameters for optimization using linear regression
  in_lm <- lm(train_obs ~ train_em)
  inpars <- c(in_lm$coef[1],
              sqrt(max(in_lm$coef[2], 0)), c = 1, d = 1)
  
  # Parameter estimation by minimizing mean TN CRPS over the training data using BFGS algorithm
  opt <- optim(inpars, fn = CRPS_TN, method = "BFGS", control = list(maxit = 300),
               train_obs = train_obs, train_em = train_em, train_ev = train_ev)
  
  # Fall-back options in case of no convergence
  if(opt$convergence > 0){
    
    # Parameter estimation by minimizing mean TN CRPS over the training data using Nelder-Mead algorithm
    opt <- optim(inpars, fn = CRPS_TN, method = "Nelder-Mead",
                 train_obs = train_obs, train_em = train_em, train_ev = train_ev)
  }
  
  # Fall-back options in case of large coefficients - regularization
  if(sum(abs(opt$par)) > 50){
    
    # Record regularized case
    reg <- 1
    
    # Parameter estimation by minimizing regularized mean TN CRPS over the training data using BFGS algorithm
    opt <- optim(inpars, fn = CRPS_TN_REG, method = "BFGS", control = list(maxit = 300),
                 train_obs = train_obs, train_em = train_em, train_ev = train_ev, lambda = lambda)
    
    # Fall-back option in case of no convergence
    if(opt$convergence > 0){
      
      # Parameter estimation by minimizing regularized mean TN CRPS over the training data using Nelder-Mead algorithm
      opt <- optim(inpars, fn = CRPS_TN_REG, method = "Nelder-Mead",
                   train_obs = train_obs, train_em = train_em, train_ev = train_ev, lambda = lambda)
    }
  }
  
  return(c(opt$par, opt$convergence, reg))
}