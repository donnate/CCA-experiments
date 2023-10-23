#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Huber : obtain matrix of Huber functions used in initialization              #
#   algorithm of FIT-SSVD, Algorithm 2 arXiv:1112.2433                         #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Inputs                                                                       #
#   x    : matrix of observed data                                             #
#   beta : degree of "Huberization"                                            #
# Outputs                                                                      #
#   Matrix of Huber functions                                                  #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
Huber <- function(x,
                  beta){
  
  X.abs <- abs(x)
  
  #--------------------------------------------------------------------------#
  # delta : the beta-quantile of the absolute values of all the entries in X #
  #--------------------------------------------------------------------------#
  delta <- quantile(x=X.abs, probs=beta)
  
  #--------------------------------------------------------------------------#
  # Create Huber rho function                                                #
  #--------------------------------------------------------------------------#
  Y <- matrix(0,nrow=nrow(x),ncol=ncol(x))
  
  matrix.le <- X.abs <= delta
  
  Y[matrix.le] <- (X.abs[matrix.le])^2
  Y[!matrix.le] <- (2*X.abs[!matrix.le]*delta - delta^2)
  
  return(Y)
}


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# select.indices                                                               #
#  Uses the Holm-Bonferroni method to identify rejected hypotheses             #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Inputs                                                                       #
#  t.vector  : hypotheses to be tested                                         #
#  alpha     : significance level (Default=0.05)                               #
# Outputs                                                                      #
#   vector, indices of t.vector that are rejected                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
select.indices <- function(t.vector, 
                           alpha=0.05){
  
  #--------------------------------------------------------------------------#
  # Population standard deviation of t-vector                                #
  #--------------------------------------------------------------------------#
  s_hat <- mad(t.vector)
  
  #--------------------------------------------------------------------------#
  # Median of t-vector                                                       #
  #--------------------------------------------------------------------------#
  mu_hat <- median(t.vector)
  
  #--------------------------------------------------------------------------#
  # Calculate p-value                                                        #
  #--------------------------------------------------------------------------#
  p.value <- 1 - pnorm(q=t.vector, mean=mu_hat, sd=s_hat)
  
  #--------------------------------------------------------------------------#
  # Adjusted p-value using Holm (1979)                                       #
  #--------------------------------------------------------------------------#
  p.holm  <- p.adjust(p=p.value, method="holm")
  
  #--------------------------------------------------------------------------#
  # Order adjusted p-values                                                  #
  #--------------------------------------------------------------------------#
  p.holm.ordered <- sort(p.holm)
  
  #--------------------------------------------------------------------------#
  # Identify rejected hypotheses                                             #
  #--------------------------------------------------------------------------#
  reject <- sum(p.holm.ordered <= alpha)
  
  if(reject == 0){
    reject <- length(p.holm)
  }
  
  index <- which( rank(p.holm) <= reject)
  
  return(index)
}