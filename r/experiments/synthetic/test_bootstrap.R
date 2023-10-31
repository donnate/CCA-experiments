x.huber_X <- cor(example$X, example$X)^2
rownorm2_X <- apply(x.huber_X, 1, sum)
mean(rownorm2_X[1:r])
mean(rownorm2_X[1:(r_pca * nnz)])
mean(rownorm2_X[(r_pca * nnz + 1):p1])
plot(rownorm2_X)

n_bootstrap <- 1000
bootstrap_norms <- matrix(0, n_bootstrap, p1)


x.huber <- atanh(cor(example$Data) -diag(p))^2
x.huber <- atanh(cor(example$X,example$Y)) ^2
rownorm2 <- apply(x.huber, 1, sum)

sigma.hat <- mad(as.vector(atanh(cor(example$Data) -diag(p))))
I.row <- get.subset(rownorm2, method = method, alpha.method = alpha.method,
                    alpha.theory = alpha.theory, sigma = sigma.hat, df = pv)

plot(rownorm2)

mean(rownorm2[1:12])
mean(rownorm2[1:(r_pca * nnz)])
mean(rownorm2[(r_pca * nnz + 1):p1])
n_bootstrap <- 1000
bootstrap_norms <- matrix(0, n_bootstrap, p1)

for(i in 1:n_bootstrap){
  print(i)
  permuted_cols<- sample(1:(p1+p2), (p1+p2), replace = TRUE)
  permuted_cov <- cor(example$Data[, permuted_cols])
  x.huber.boot <- atanh(permuted_cov)^2
  bootstrap_norms[i, ] <- apply(x.huber.boot, 1, sum) # L2 norm
}

permuted_colsX

p.boot <- sapply(rownorm2, function(x){mean(rownorm2 > x)}) 
quantile(bootstrap_norms[, 1], 0.95)
rownorm2

which(rownorm2 > quantile(rownorm2, 0.95))

U = matrix(0, p1, r)
U[set_sel[which(set_sel<(p1+1))],] = test$xcoef[,1:r]
V = matrix(0, p2, r)
V[set_sel[which(set_sel>(p1))] - p1,] = test$ycoef[,1:r]
evaluate_results(Uhat= U, 
                 Vhat = V, 
                 example = example, 
                 name_method="test", 
                 overlapping_amount=overlapping_amount,
                 lambdax= NA,
                 lambday = NA, 
                 thres = THRES,
                 it=it,
                 normalize_diagonal=normalize_diagonal,
                 criterion=criterion,
                 r_pca = r_pca, nnz=nnz,
                 signal_strength=signal_strength)

U = matrix(0, p2, r)