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