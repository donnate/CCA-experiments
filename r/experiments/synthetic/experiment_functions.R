library(MASS)
library(stats)
library(CVXR)
library(geigen)
library(pracma)
library(tidyverse)
library(CCA)
library(VGAM)
library(matlib)
library(PMA)
library(mvtnorm)
library(glmnet)
library(caret)
#Example of TGD on Sparse CCA
#n = 500, p1 = p2 = 100, s_u = s_v = 5
#k = 20, eta = 0.0025, lambda =0.01, T = 12000
wd = getwd()

source("r/GCA/utils.R")
source("r/GCA/gca_to_cca.R")
source("r/GCA/init_process.R")
source("r/GCA/sgca_init.R")
source("r/GCA/sgca_tgd.R")
source("r/GCA/subdistance.R")


cv_function <- function(X, Y, 
                        kfolds=10, initu, initv,
                        lambdax, adaptive=TRUE, normalize=FALSE,
                        criterion="prediction") {
  # define empty vector to store results
  folds <- createFolds(1:nrow(Y), k = kfolds, list = TRUE, returnTrain = FALSE)
  rmse <- numeric(length = kfolds)
  nnz  <- numeric(length = kfolds)
  # loop over folds
  for (i in seq_along(folds)) {
    # split data into training and validation sets
    X_train <- X[-folds[[i]], ]
    Y_train <- Y[-folds[[i]], ]
    X_val <- X[folds[[i]], ]
    Y_val <- Y[folds[[i]], ]
    
    # fit model on training data with hyperparameters
    tryCatch(
    {
    model <- adaptive_lasso(X_train, Y_train %*% initv, initu, adaptive=adaptive, 
                         lambdax, 
                         max.iter=5000, 
                         max_k = 10, verbose = FALSE, THRESHOLD=1e-5)
    
    # make predictions on validation data
    # compute RMSE on validation data
    ##### Normalize Uhat
    if (normalize == FALSE){
      if (criterion=="prediction"){
        rmse[i] <- sum((X_val %*% model$Uhat - Y_val%*% initv)^2)
        if (norm(model$Uhat) == 0){ ####prevents selecting values that would make everything 0
          rmse[i] <- 1e8
        }
      }else{
        rmse[i] <- sum(abs(cor(X_val %*% model$Uhat, Y_val%*% initv)))
      }
      nnz[i] <- sum(apply(model$Uhat^2, 1, sum) >1e-4)
    }else{
      sol <- gca_to_cca(rbind(model$Uhat, initv), 
                        cov(rbind(X_val, Y_val)), pp)
      if (criterion=="prediction"){
      rmse[i] <- sum((X_val %*% sol$u - Y_val%*% initv)^2)
      }else{
        rmse[i] <- sum(abs(cor(X_val %*%sol$u, Y_val%*% initv)))
      }
      nnz[i] <- sum(apply(sol$u^2, 1, sum) >1e-4)
     }
    },
    error = function(e) {
      # If an error occurs, assign NA to the result
      rmse[i] <- NA
    })
  }
  
  # return mean RMSE across folds
  if (mean(is.na(rmse)) == 1){
      return(1e8)
   }else{
  return(median(rmse, na.rm=TRUE))
   }
}


cv_function_tgd <- function(X, Y, Mask, kfolds=5, ainit,
                        lambda, r=2, k=20,  maxiter=1000, eta=0.001,
                        convergence=1e-3, normalize=FALSE, criterion="prediction") {
  # define empty vector to store results
  folds <- createFolds(1:nrow(Y), k = kfolds, list = TRUE, returnTrain = FALSE)
  rmse <- numeric(length = kfolds)
  p1 <- dim(X)[2]
  p2 <- dim(Y)[2]
  p <- p1 + p2;
  n <- nrow(X)
  pp <- c(p1,p2);
  S0 = cov(cbind(X, Y))
  
  #init <- gca_to_cca(ainit, S0, pp)
  # loop over folds
  for (i in seq_along(folds)) {
    # split data into training and validation sets
    X_train <- X[-folds[[i]], ]
    Y_train <- Y[-folds[[i]], ]
    X_val <- X[folds[[i]], ]
    Y_val <- Y[folds[[i]], ]
    
    S = cov(cbind(X_train, Y_train))
    sigma0hat = S * Mask
    
    # fit model on training data with hyperparameters
    tryCatch(
    {
    final = sgca_tgd(A=S, B=sigma0hat,
             r=r,ainit, k, lambda = lambda, eta=eta,
             convergence=convergence,
             maxiter=maxiter, plot = FALSE, 
             scale=NULL)
    final <- gca_to_cca(final, S, pp)
    
    # make predictions on validation data
    # compute RMSE on validation data
    if (normalize == FALSE){
      if (criterion=="prediction"){
        rmse[i] <- sum((X_val %*% final$u - Y_val%*% initv)^2)
      }else{
        rmse[i] <- sum(abs(cor(X_val %*% final$u, Y_val%*% initv)))
      }
    }else{
      sol <- gca_to_cca(rbind(final$u, initv), 
                        cov(rbind(X_val, Y_val)), pp)
      if (criterion=="prediction"){
        rmse[i] <- sum((X_val %*% sol$u - Y_val%*% initv)^2)
      }else{
        rmse[i] <- sum(abs(cor(X_val %*%sol$u, Y_val%*% initv)))
      }
    }
    },
    error = function(e) {
      # If an error occurs, assign NA to the result
      rmse[i] <- NA
    })
  }
  
  # return mean RMSE across folds
  if (mean(is.na(rmse)) == 1){
      return(1e8)
   }else{
  return(median(rmse, na.rm=TRUE))
   }
}

preselection <-function(Data, CorrelationMat, p1, r, alpha){
  p = ncol(Data)
  n=  nrow(Data)
  p2 =  p-p1
  t = apply(CorrelationMat -diag(diag(CorrelationMat)), 1, 
            function(x){max(x^2)})
  J = order(-t)[1: ceiling(alpha  * n/(r))]
  set_u = J[which(J <= p1)]
  set_v = J[which(J > p1)]
  t=CCA::cc(as.matrix(Data[,set_u]), as.matrix(Data[, set_v]))
  Uhat = matrix(0, p, r)
  Uhat[set_u, ] =  t$xcoef[,1:r]
  Uhat[set_v, ] =  t$ycoef[,1:r]
  return(Uhat)
}


## Running initialization using convex relaxation

#(Sigma,Sigma0, lambda, rho, eta=0.001, nu=1,epsilon=5e-3,maxiter=1000,trace=FALSE)


