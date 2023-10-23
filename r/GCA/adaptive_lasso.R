library(CVXR)


####

adaptive_lasso <- function(X, Z, beta0, adaptive=TRUE,
                         lambda, max.iter=5000,
                         max_k = 10, verbose = FALSE,
                         THRESHOLD=1e-5){
  
  p <- nrow(beta0)
  r <- ncol(beta0)
  ## Define_Parameters

  ## Build_Penalty_Terms
  
  if (adaptive){
    W0 = sqrt(apply(beta0^2, 1,sum))
    W0[which(W0< THRESHOLD)] = THRESHOLD
    W  <- diag(W0)
    #invW <- diag(sapply(W0, function(x){ifelse(x<ZERO_THRESHOLD, 1/ZERO_THRESHOLD, 1/x)}))
  }else{
    W <- diag(rep(1,p))
    #invW <- diag(rep(1,p))
  }
  
  Uhat <- Variable(p, r)
  penalty_term <- sum(cvxr_norm(Uhat, 2, axis = 2))
  
  #lambda=100
  objective <- 1 / 2 * sum_squares(X %*% W %*% Uhat  - Z) + lambda  * penalty_term
  
  ## Define_and_Solve_Problem
  prob <- Problem(Minimize(objective))
  thresh = 1e-8
  result <- solve(prob, FEASTOL = thresh, RELTOL = thresh, ABSTOL = thresh, verbose = verbose)
  #result <- psolve(prob, verbose = TRUE, num_iter =max.iter)
  ## Return_Values
  Uhat <- W%*% result$getValue(Uhat)
  #### This is just a quirk of CVXR
  #plot(beta0[,1],Uhat[,1])
  ## Zero out stuff before returning
  Uhat[abs(Uhat) < THRESHOLD] <- 0.0
  #### Normalize appropriately
  Uhat <- Uhat %*% sqrtm(t(Uhat) %*% t(X) %*% X %*% Uhat)$Binv
  list(
    Uhat = Uhat,
    criterion = result$value)
}
