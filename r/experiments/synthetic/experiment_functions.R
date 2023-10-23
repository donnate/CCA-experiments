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
setwd("GCA")
source("utils.R")
source("gca_to_cca.R")
source("init_process.R")
source("sgca_init.R")
source("sgca_tgd.R")
source("subdistance.R")
source("adaptive_lasso.R")

setwd(wd)
source('experiments/alternative_methods/SAR.R')
source('experiments/alternative_methods/Parkhomenko.R')
source('experiments/alternative_methods/Witten_CrossValidation.R')
source('experiments/alternative_methods/Waaijenborg.R')

setwd(wd)

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


evaluate_results <- function(Uhat, Vhat, example, name_method, overlapping_amount, 
                             thres = 0.0001, lambdax= NULL,lambday = NULL, it=1,
                             normalize_diagonal=TRUE,
                             criterion="prediction", r_pca = 0, nnz=0){
  Uhat_tot = rbind(Uhat, Vhat)
  U_tot = rbind(example$u, example$v)
  p1 = ncol(example$X)
  p2 = ncol(example$Y)
  n = nrow(example$X)
  r = ncol(example$u)
  sparsity = length(which(apply(example$u^2, 1, sum)>0))
  silly_benchmark = subdistance(matrix(0, p1, 2), example$u)
  data.frame("method" = name_method,
             "exp" = it,
             "n" = n,
             "nnz" = nnz,
             "p1" = p1,
             "p2" = p2,
             "sparsity" = sparsity,
             "r" = r,
             "r_pca" = r_pca,
             "criterion" = criterion,
             "overlapping_amount" = overlapping_amount,
             "zero_benchmark" = silly_benchmark,
             "nb_discoveries" = sum(apply(Uhat_tot^2, 1, sum)>0),
             "nb_real_discoveries" = sum(apply(Uhat_tot^2, 1, sum)>thres),
             "param1" = lambdax,
             "param2" = lambday,
             "normalize_diagonal" = normalize_diagonal,
             "distance_tot" = subdistance(Uhat_tot, U_tot),
             "distance_U" = subdistance(Uhat, example$u),
             "distance_V" = subdistance(Vhat, example$v),
             "sinTheta_tot" = sinTheta(Uhat_tot, U_tot),
             "sinTheta_U" = sinTheta(Uhat, example$u),
             "sinTheta_V" = sinTheta(Vhat, example$v),
             "sinTheta_U_0" = sinTheta(matrix(0, p1, r), example$u),
             "sinTheta_V_0" = sinTheta(matrix(0, p2, r), example$v),
             "prediction_tot" = mean((example$X %*% Uhat - example$Y %*% Vhat)^2),
             "prediction_U" = mean((example$X %*% Uhat - example$X %*% example$u)^2),
             "prediction_V" = mean((example$Y %*% Vhat - example$Y %*% example$v)^2),
             "TPR" =TPR(apply(Uhat_tot^2, 1, sum), apply(U_tot^2, 1, sum), tol=thres),
             "TNR" = TNR(apply(Uhat_tot^2, 1, sum), apply(U_tot^2, 1, sum), tol=thres),
             "FPR" = FPR(apply(Uhat_tot^2, 1, sum), apply(U_tot^2, 1, sum), tol=thres),
             "FNR" = FPR(apply(U_tot^2, 1, sum),apply(Uhat_tot^2, 1, sum), tol=thres)
  )
}
