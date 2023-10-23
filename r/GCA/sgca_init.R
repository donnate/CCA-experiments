# Function for Intialization via generalized fantope projection
# Inputs:
# =======
# A, B:       Pair of matrix to calculate generalized eigenspace (sigma, sigma0)
# nu:         Parameter of ADMM, default set to 1
# K:          nuclear norm constraint, equal to r
# rho:     penalty parameter on the l_1 norm of the solution, scaled by
#             sqrt(log(max(p1,p2))/n)
# epsilon:    tolerance level for convergence in ADMM
# maxiter:    maximum number of iterations in ADMM
# trace:      if set to True will print all iterations 

# Outputs:
# ========
# $Pi:     optimum of the convex program

sgca_init <-
  function(A,B,rho,K,nu=1,epsilon=5e-3,maxiter=1000,trace=FALSE){
    p <- nrow(B)
    eigenB <- eigen(B)
    sqB <- eigenB$vectors%*%sqrt(diag(pmax(eigenB$values,0)))%*%t(eigenB$vectors)	
    tau <- 4*nu*eigenB$values[1]^2	
    criteria <- 1e10
    i <- 1
    # Initialize parameters
    H <- Pi <- oldPi <-  diag(1,p,p)
    Gamma <- matrix(0,p,p)
    # While loop for the iterations
    while(criteria > epsilon && i <= maxiter){
      for (j in 1:20){
        Pi <- updatePi(B,sqB,A,H,Gamma,nu,rho,Pi,tau)
      }
      #Pi <- updatePi(B,sqB,A,H,Gamma,nu,lambda,Pi,tau)
      
      H <- updateH(sqB,Gamma,nu,Pi,K)
      Gamma <- Gamma + (sqB%*%Pi%*%sqB-H) * nu	
      criteria <- sqrt(sum((Pi-oldPi)^2))
      oldPi <- Pi
      i <- i+1
      if(trace==TRUE)
      {
        print(i)
        print(criteria)
      }
    }
    return(list(Pi=Pi,H=H,Gamma=Gamma,iteration=i,convergence=criteria))
    
  }



soft_threshold <- function(A, lambda){
  return( A * ((abs(A) >lambda) * (1-lambda* sign(A))))
}

sgca_init_proximal <-
  function(Sigma,Sigma0, lambda, rho, eta=0.001, nu=1,epsilon=5e-3,maxiter=1000,trace=FALSE){
    
    p <- nrow(Sigma)
    criteria <- 1e10
    i <- 1
    # Initialize parameters
    oldA<-A <- diag(nrow=p, ncol=p)
    # While loop for the iterations
    while(criteria > epsilon && i <= maxiter){
      print(c(i,criteria))
      gradient <- -Sigma + rho * ( Sigma0 %*% A *  Sigma0  - Sigma0 )
      A <- soft_threshold(A - eta * gradient, lambda)
      criteria <- sqrt(sum((A-oldA)^2))
      oldA <- A
      i<- i+1
    }
    return(list(A=A,iteration=i,convergence=criteria))
    
  }


