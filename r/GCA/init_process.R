# function to convert Pi from Fantope to input of TGD
# Inputs:
# =======
# Pi:         Output of sgca_init$Pi
# r:          Latent dimension

# Outputs:
# ========
# ainit:   Initialization for the generalized eigenspace
init_process <-
  function(Pi, r){
    ainit = svd(Pi)
    uest <- ainit$u
    dest <- diag(ainit$d)
    if (r == 1){
      ainit <- uest[,1] * sqrt(dest[1:r,1:r])
    } else
      ainit <- uest[,1:r] %*% pracma::sqrtm(dest[1:r,1:r])$B
    return(ainit)
  }