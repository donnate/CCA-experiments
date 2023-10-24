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
