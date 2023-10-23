library(magic)
library(MASS)
library(stats)
library(geigen)
library(pracma)
# Example pf Sparse PCA using TGD
# r=3, s=20
# we first generate an orthogonal matrix of size 20*3 whose
# row l2 norm are same for all rows, let it be A1
# to generate A1, we notice A0 = [1,1,-1;1,-1,-1;1,1,1;1,-1,1] has same row norm
# We use A0 to get one block for A1. 
# Then we generate 5 blocks in this way and concantenate them to get A1

disp('  ');
disp('--------------------------------------');
disp('SPCA Example: Toeplitz covariance matrix'); 
disp('n = 500, p = 200, s=20, leading eigenvalues are 7, 5, 3.');

n <- 500;
p <- 200;
k <- 40;
s <- 20;
r <- 3

A0=t(matrix(c(1,1,-1,1,-1,-1,1,1,1,1,-1,1), nrow = 3, ncol = 4));
A1 <- rbind(A0, A0,A0,A0,A0);

d <- diag(c(1, 1, 1)) #eigenvalues are 5, 5, 5
A1 = sqrt(4/17) * A1;
upper_block <- A1 %*% d %*% t(A1) + (5/17) * diag(s)
R <- adiag(upper_block, diag(p-s))

#generate sigma_0 with random numbers in [0.1,1] on the diagonal 
Sigma0 = diag(as.vector(rand(p,1)*0.9+0.1));
sq = sqrtm(Sigma0)$B
Sigma = sq %*% R %*% sq

# Generate data for SPCA
X <-mvrnorm(n, rep(0, p) , Sigma)
S <- cov(X)
sigma0hat <- S * Mask

# Estimate the subspace spanned by the largest eigenvector using convex relaxation
# First calculate ground truth
result = geigen(Sigma,Sigma0)
evalues <- result$ values
evectors <-result$vectors
evectors <- evectors[,p:1]
a <- evectors[,1:r]
scale <- a %*% sqrtm(diag(r)+t(a) %*% Sigma %*% a/lambda)$B


# Running initialization using convex relaxation
ag <- sgca_init(A=S, B=sigma0hat, lambda=0.5*sqrt(log(p)/n),K=r ,nu=1,trace=FALSE)
ainit = svd(ag$Pi)
uest <- ainit$u
dest <- diag(ainit$d)
if (r == 1){
  ainit <- uest[,1] * sqrt(dest[1:r,1:r])
} else
  ainit <- uest[,1:r] %*% sqrtm(dest[1:r,1:r])$B


print('The initial error is')
print(subdistance(ainit, a))

#perform hard thresholding
ainit <- hard(ainit, k)

## Running Thresholded Gradient Descent
final <- sgca_tgd(A=S, B=sigma0hat,r,ainit,k,lambda = 0.01, eta=0.001,convergence=1e-6,maxiter=15000, plot = TRUE)
print('The final error is')
print(subdistance(final, a))



