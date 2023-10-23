library(MASS)
library(stats)
library(geigen)
library(pracma)
# Example on Sparse Generalized Correlation Analysis
print('  ');
print('--------------------------------------');
print('SGCA Example: 3 high dimensional datasets with Toeplitz covariance matrix n = 300, p = 200'); 
pp <- c(100, 50, 50);
p <- sum(pp)
s  <- c(1, 6, 11, 24, 31);
r <- 3;
k <- 20;
eta <- 0.001;
n = 300;
lambda <- 0.01;
max_iter <- 15000;
# Generate sparse generalized canonical vectors
u1 <- matrix(0, pp[1],r);
u2 <- matrix(0, pp[2],r);
u3 <- matrix(0, pp[3],r);

# Generate Sigma and Sigma0
Sigma = diag(sum(pp));
T1 = toeplitz(0.5^(0:(pp[1]-1)));
T2 = toeplitz(0.7^(0:(pp[2]-1)));
T3 = toeplitz(0.9^(0:(pp[3]-1)));
u1[s,1:r] <- matrix( rnorm(length(s)*r,mean=0,sd=1), length(s), r)
u1 <- u1 %*%(sqrtm(t(u1[s,1:r]) %*% T1[s,s] %*% u1[s,1:r])$Binv);
u2[s,1:r] = matrix(rnorm(length(s)*r,mean=0,sd=1), length(s), r) 
u2 <- u2 %*% (sqrtm(t(u2[s,1:r]) %*% T2[s,s] %*% u2[s,1:r])$Binv);
u3[s,1:r] = matrix( rnorm(length(s)*r,mean=0,sd=1), length(s), r) 
u3 <- u3 %*%(sqrtm(t(u3[s,1:r]) %*% T3[s,s] %*% u3[s,1:r])$Binv);

idx1 = 1:pp[1];
idx2 = (pp[1]+1):(pp[1]+pp[2]);
idx3 = (pp[1]+pp[2]+1):(pp[1]+pp[2]+pp[3]);
Sigma[idx1, idx1] = T1;
Sigma[idx2, idx2] = T2;
Sigma[idx3, idx3] = T3;
SigmaD = Sigma;

Sigma[idx1, idx2] <- T1 %*% u1 %*% t(u2) %*% T2;
Sigma[idx1, idx3] <- T1 %*% u1 %*% t(u3) %*% T3;
Sigma[idx2, idx3] <- T2 %*% u2 %*% t(u3) %*% T3;
Sigma = Sigma + t(Sigma) - SigmaD;

Mask = matrix(0, sum(pp),sum(pp));
Mask[idx1,idx1] <- matrix(1,pp[1],pp[1]);
Mask[idx2,idx2] <- matrix(1,pp[2],pp[2]);
Mask[idx3,idx3] <- matrix(1,pp[3],pp[3]);
Sigma0 = Sigma * Mask;

# Generate data for SGCA
X <-mvrnorm(n, rep(0, p) , Sigma)
S <- cov(X)
sigma0hat <- S * Mask

# Estimate the subspace spanned by the largest eigenvector using convex relaxation
# First calculate ground truth
result = geigen::geigen(Sigma,Sigma0)
evalues <- result$values
evectors <-result$vectors
evectors <- evectors[,p:1]
a <- evectors[,1:r]
scale <- a %*% sqrtm(diag(r)+t(a) %*% Sigma %*% a/lambda)$B;



# Running initialization using convex relaxation
ag <- sgca_init(A=S, B=sigma0hat, rho=0.5*sqrt(log(p)/n),K=r ,nu=1,trace=FALSE, plot=TRUE,
                scale=scale)
ainit <- init_process(ag$Pi, r)



print('The initial error is')
print(subdistance(ainit, a))



## Running Thresholded Gradient Descent
final <- sgca_tgd(A=S, B=sigma0hat,r,ainit,k,lambda = 0.01, eta=0.001,convergence=1e-6,maxiter=15000, plot = TRUE)
print('The final error is')
print(subdistance(final, a))








