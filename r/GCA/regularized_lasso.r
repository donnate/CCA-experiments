library(CVXR)


####

regularized_lasso <- function(X, Z,  Gamma,
                         lambda, lambda2, max.iter=5000,
                         max_k = 10, verbose = FALSE,
                         ZERO_THRESHOLD=1e-5){
  
  p <- nrow(beta0)
  r <- ncol(beta0)
  ## Define_Parameters

  ## Build_Penalty_Terms

  Uhat <- Variable(p, r)
  penalty_term <- sum(cvxr_norm(Gamma %*% Uhat, 2, axis = 1))
  
  #lambda=100
  objective <- 1 / 2 * sum_squares(X %*% Uhat  - Z) + lambda  * penalty_term + lambda1
  ## Define_and_Solve_Problem
  prob <- Problem(Minimize(objective))
  thresh = 1e-8
  result <- solve(prob, FEASTOL = thresh, RELTOL = thresh, ABSTOL = thresh, verbose = TRUE)
  #result <- psolve(prob, verbose = TRUE, num_iter =max.iter)
  ## Return_Values
  Uhat <- result$getValue(Uhat)
  #### This is just a quirk of CVXR
  #plot(beta0[,1],Uhat[,1])
  ## Zero out stuff before returning
  Uhat[abs(Uhat) < ZERO_THRESHOLD] <- 0.0
  list(
    Uhat = Uhat,
    criterion = result$value)
}
