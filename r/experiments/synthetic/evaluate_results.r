source("r/metrics.R")

evaluate_results <- function(Uhat, Vhat, example, 
                             name_method, overlapping_amount, 
                             thres = 0.0001, lambdax= NULL,
                             lambday = NULL, it=1,
                             normalize_diagonal=TRUE,
                             criterion="prediction", r_pca = 0, 
                             nnz=0){
  Uhat_tot = rbind(Uhat, Vhat)
  U_tot = rbind(example$u, example$v)
  p1 = ncol(example$X)
  p2 = ncol(example$Y)
  n = nrow(example$X)
  r = ncol(example$u)
  sparsity = length(which(apply(example$u^2, 1, sum)>0))
  silly_benchmark = subdistance(matrix(0, p1, r), example$u)
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