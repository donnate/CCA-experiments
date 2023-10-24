#########################################
# AUXILIARLY FUNCTIONS SIMULATION STUDY #
#########################################


principal_angles <- function(a, b){
  ### Calculate principal angles between subspace spanned by the columns of a and the subspace spanned by the columns of b
  
  angles=matrix(0, ncol=ncol(a), nrow=1)
  qa= qr.Q(qr(a))
  qb= qr.Q(qr(b))
  C=svd(t(qa)%*%qb)$d
  rkA=qr(a)$rank;rkB=qr(b)$rank
  if (rkA<=rkB){
    B = qb - qa%*%(t(qa)%*%qb);
  } else {B = qa - qb%*%(t(qb)%*%qa)}
  S=svd(B)$d
  S=sort(S)
  
  for (i in 1:min(rkA, rkB)){
    if (C[i]^2 < 0.5) {angles[1, i]=acos(C[i])}
    else if (S[i]^2 <=0.5) {angles[1, i]=asin(S[i])}
  }
  angles=t(angles)
  
  ##OUTPPUT
  out=list(angles=angles)
}

TPR  <-  function(A, B, tol=1e-4){
  # This is a function that compares the structure of two matrices A and B
  # It outputs the number of entries that A and B have in common that are different from zero
  # A and B need to have the same number of rows and columns
  A[which(abs(A) <tol)] =0
  B[which(abs(B) <tol)] =0
  out  <-  sum((A!=0)*(B!=0))/max(1,sum(A!=0))
  return(out)
}

FPR  <-  function(A, B, tol=1e-4){
  # This is a function that compares the structure of two matrices A and B
  # It outputs the number of entries that A and B have in common that are different from zero
  # A and B need to have the same number of rows and columns
  A[which(abs(A) <tol)] =0
  B[which(abs(B) <tol)] =0
  out  <-  sum((A!=0)*(B==0))/max(1,sum(A!=0))
  return(out)
}

TNR  <-  function(A, B, tol=1e-4){
  # This is a function that compares the structure of two matrices A and B
  # It outputs the number of entries that A and B have in common that are zero #
  # A and B need to have the same number of rows and columns
  A[which(abs(A) <tol)] =0
  B[which(abs(B) <tol)] =0
  out  <-  sum((A==0)*(B==0))/max(1,sum(A==0))
  return(out)
}


subdistance <- function(A, B){
  svdresult = svd(t(A) %*% B);
  U = svdresult$u
  V = svdresult$v
  O = U %*% t(V);
  l = norm(A %*% O-B, type = c('F'));
  return(l)
}

sinTheta<- function(U, V){
  l = 1/sqrt(2) * norm(U %*% t(U)-V %*% t(V), type = c('F'));
  return(l)
}



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

