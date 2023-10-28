library(tidyverse)
example <- generate_example_none_trivial_pca(n, p1, p2, 
                                             r_pca = r_pca, 
                                             nnzeros = nnz,
                                             theta = diag(seq(from = 0.95, to=0.85, 
                                                              length.out = r)),
                                             lambda_pca = 1,
                                             r=r, 
                                             overlapping_amount = overlapping_amount,
                                             normalize_diagonal=normalize_diagonal)


center_D = scale(example$Data, scale=TRUE)
test = (center_D[1, ] %o% center_D[1, ])^2
for (i in 2: n){
  test = test + (center_D[i, ] %o% center_D[i, ])^2
}
test = test/(n-1) 
  
ggplot(data.frame("x" = 1:nrow(example$u),
                  "true" = c(rep(1, nnz), rep(0, p1-nnz)),
                  #"norm_u"= apply(example$u^2, 1, sum),
                  #"norm_Sigma"= apply(example$Sigma[1:p1, (p1+1):(p1+p2)]^2, 1, sum),
                  "norm_sxy"= apply(cov(example$Data)[1:p1, (p1+1):(p1+p2)]^2, 1, sum),
                  "norm_sxy_tr"= apply(atanh(cov(example$Data)[1:p1, (p1+1):(p1+p2)])^2, 1, sum),
                  "norm_test"= apply(test[1:p1, (p1+1):(p1+p2)], 1, sum),
                  "norm_sxy4"= apply(cov(example$Data)[1:p1, (p1+1):(p1+p2)]^4, 1, sum),
                  "norm_sxy2_2"= apply(cov(example$Data)[1:p1, (p1+1):(p1+p2)]^2, 1, sum)^2,
       "norm_sxy6"= apply(cov(example$Data)[1:p1, (p1+1):(p1+p2)]^2, 1, sum)),
       aes(x=x, y=norm_sxy4, colour=as.factor(true)))+
  geom_point() 
  #geom_point(aes(x=x, y=norm_sxy, colour=as.factor(true))) 
  #geom_point(aes(x=x, y=norm_sxy4- norm_sxy2_2) 
  #geom_point(aes(x=x, y=norm_sxy2_2, colour=as.factor(true)))  
  #geom_point(aes(x=x, y=norm_Sigma, colour="Sigma_xy")) 



