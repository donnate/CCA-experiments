# Function for Thredholded Gradient Descent
# Inputs:
# =======
# A, B:         Pair of matrix to calculate generalized eigenspace (sigma, sigma0)
# r:            Latent dimension
# init:         Initialize estimator of generlized eigenspace, obtained from sgca_int. Need to be k-sparse
# lambda:       penalty parameter of Lagrangian function f(L), default to be 0.01
# convergence:  tolerance level for convergence in TGD
# maxiter:      maximum number of iterations in TGD
# plot:         if set to True will plot intermediate iterations, need to specify scale variable V (as shown in paper)
#               default set to False

# Outputs:
# ========
# final_estimate:   final estimation of leading r sparse generalized eigenspace
library(pracma)
sgca_tgd <-
  function(A, B, r, init, k, lambda = 0.01, eta=0.01, convergence=1e-3, maxiter=10000, plot = FALSE, scale=NULL){
    #perform hard thresholding
    init <- hard(init, k, r)
    u <- init
    criteria <- 1e10
    iter <- 0
    error <- rep(0, maxiter)
    # renormalization 
    ut <- init %*% sqrtm(diag(r)+t(u) %*% A %*% u/lambda)$B;
    
    while(criteria > convergence && iter <= maxiter){
      #perform gradient descent
      grad <- -A %*% ut + lambda * B %*% ut %*% (t(ut) %*% B %*% ut- diag(r));
      vt <- ut - eta * grad
      
      
      # Perform truncation
      vt <- hard(vt, k, r)
      
      criteria <- sqrt(sum((ut-vt)^2))
      if (is.na(criteria)){
         criteria =0
      } 
      ut <- vt
      iter <- iter+1
      if (plot & (is.null(scale)==FALSE)){
        error[iter] <- subdistance(vt, scale)
      }
    }
    if (plot){
      plot(error[1:iter], type='l',  main="Matrix distance and iterations", 
           xlab="Number of iterations", ylab="Matrix distance",)
    }
    final_estimate <- ut %*% sqrtm(t(ut) %*% B %*% ut)$Binv
    return(final_estimate)
  }
