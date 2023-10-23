library(MASS)
library(stats)
library(geigen)
library(pracma)
#Example of TGD on Sparse CCA
#n = 500, p1 = p2 = 100, s_u = s_v = 5
#k = 20, eta = 0.0025, lambda =0.01, T = 12000
setwd("~/Documents/Sparse-Generalized-Correlation-Analysis/src/R Code/")
source("utils.R")
source("gca_to_cca.R")
source("init_process.R")
source("sgca_init.R")
source("sgca_tgd.R")
source("subdistance.R")
print('--------------------------------------');
print('SCCA Example: Toeplitz covariance matrix'); 
print('n = 500, p1 = p2 = 100, s_u = s_v = 5');
print('lambda1 = 0.9, lambda2 = 0.8');

n  <- 500;
p1 <- 100;
p2 <- 100;
p <- p1 + p2;
pp <- c(p1,p2);
s  <- sample(1:min(p1,p2),5);
theta <- diag( c(0.9,  0.8) );
r <- 2
print('--------------------------------------');
print('Generating data ...');
a <- 0.3;
Sigma <- diag(p1+p2)

# generate covariance matrix for X and Y
T1 = toeplitz(a^(0:(pp[1]-1)));
T1[which(T1<1e-6)] = 0
Sigma[1:p1, 1:p1] = T1;
Tss = T1[s, s];
u = matrix(0, pp[1], r)
u[s,1:r] <- randi(c(5), n = length(s), m = r) - 3
u <- u %*%(sqrtm(t(u[s,1:r]) %*% Tss %*% u[s,1:r])$Binv)

T2 = toeplitz(a^(0:(pp[2]-1)));
Sigma[(p1+1):(p1+p2), (p1+1):(p1+p2)] = T2;
Tss = T2[s, s];
v = matrix(0, pp[2], r)
v[s,1:r] <- randi(c(5), n = length(s), m = r) - 3
v <- v %*%(sqrtm(t(v[s,1:r]) %*% Tss %*% v[s,1:r])$Binv)

Sigma[(p1+1):(p1+p2), 1:p1] = T2 %*%  v  %*% theta %*% t(u) %*% T1;
Sigma[1:p1, (p1+1):(p1+p2)] = t(Sigma[(p1+1):(p1+p2), 1:p1])

Sigmax = Sigma[1:p1,1:p1];
Sigmay = Sigma[(p1+1):p,(p1+1):p];

#Generate Multivariate Normal Data According to Sigma
Data = mvrnorm(n, rep(0, p), Sigma);

X = Data[,1:p1];
Y = Data[,(p1+1):(p1+p2)];

print('Data generated.');
print('--------------------------------------');

Mask = matrix(0, p, p);
idx1 = 1:pp[1];
idx2 = (pp[1]+1):(pp[1]+pp[2]);
Mask[idx1,idx1] <- matrix(1,pp[1],pp[1]);
Mask[idx2,idx2] <- matrix(1,pp[2],pp[2]);
Sigma0 = Sigma * Mask;

S <- cov(Data)
sigma0hat <- S * Mask

# Estimate the subspace spanned by the largest eigenvector using convex relaxation and TGD
# First calculate ground truth
result = geigen(Sigma,Sigma0)
evalues <- result$ values
evectors <-result$vectors
evectors <- evectors[,p:1]
a <- evectors[,1:r]

## Running initialization using convex relaxation
ag <- sgca_init(A=S, B=sigma0hat, rho=sqrt(log(p)/n),K=r ,nu=1,trace=FALSE, maxiter = 30)
ainit <- init_process(ag$Pi, r)


## Perform TGD

lambda <- 0.01
k <- 20
scale <- a %*% sqrtm(diag(r)+t(a) %*% Sigma %*% a/lambda)$B;
final <- sgca_tgd(A=S, B=sigma0hat,r,ainit,k,lambda = 0.01, eta=0.00025,convergence=1e-6,maxiter=12000, plot = TRUE)

init <- gca_to_cca(ainit, S, pp)
final <- gca_to_cca(final, S, pp)
initu<- init$u
initv <- init$v
finalu <- final$u
finalv <- final$v
sqx <- sqrtm(Sigmax)$B
sqy <- sqrtm(Sigmay)$B

print('Initial prediction error of U is')
print( subdistance(sqx %*% initu, sqx %*% u)^2)
print('Initial prediction error of V is')
print( subdistance(sqy %*% initv, sqy %*% v)^2)

print('The final prediction error on U is')
print( subdistance(sqx %*% finalu , sqx %*% u)^2)

print('The final prediction error on V is')
print( subdistance(sqy %*% finalv, sqy %*% v)^2)


#### Try Chao's approach instead
#### 

res = group_lasso(S[1:100, 1:100], S[1:100, 101:200], ainit[1:100,],
                 lambda= 0.01, max.iter=1000,
                  verbose = TRUE)

# res_y = group_lasso(X,  Y%*%ainit[1:100,],
#                     0.1,
#                     max_k = max_k, max.iter=max.iter,
#                     verbose = verbose)


print('The final prediction error on U is')
print( subdistance(sqx %*% Uhat2 , sqx %*% u)^2)

group_lasso <- function(X, YV,
                        lambda, max.iter=5000, verbose = FALSE,
                        ZERO_THRESHOLD=1e-5){

  p <- ncol(X)
  k <- ncol(YV)
  ## Define_Parameters
  Uhat <- Variable(p, k)
  ## Build_Penalty_Terms
  penalty_term2 <- sum(cvxr_norm(Uhat, 2, axis = 1))
  
  ## Compute_Fitted_Value
  ## Build_Objective
  objective <- 1 / 2 * sum_squares(X %*% Uhat - YV) + lambda  * penalty_term2
  ## Define_and_Solve_Problem
  prob <- Problem(Minimize(objective))
  result <- psolve(prob, verbose = TRUE, num_iter =max.iter)
  ## Return_Values
  Uhat <- result$getValue(Uhat)
  
  ## Zero out stuff before returning
  Uhat[abs(Uhat) < ZERO_THRESHOLD] <- 0.0
  list(
    Uhat = Uhat,
    criterion = result$value)
}


Uhat2 =  Uhat %*% (sqrtm(t(Uhat) %*% S[1:100, 1:100] %*% (Uhat)))$Binv
for (lambda in c(0.0001, 0.001, 0.01, 0.1, 1, 10, 15, 20,25, 30, 40, 50)){
  res = group_lasso(X, YV, ainit[1:100,],
                    lambda= lambda, max.iter=100,
                    verbose = TRUE)
  Uhat = res$Uhat
  Uhat2 =  Uhat %*% (sqrtm(t(Uhat) %*% S[1:100, 1:100] %*% (Uhat)))$Binv
  print('The final prediction error on U is')
  print( subdistance(sqx %*% Uhat2 , sqx %*% u)^2)
  
}

print('Initial prediction error of U is')
print( subdistance(sqx %*% initu, sqx %*% u)^2)
print('Initial prediction error of V is')
print( subdistance(sqy %*% initv, sqy %*% v)^2)

print('The final prediction error on U is')
print( subdistance(sqx %*% Uhat2 , sqx %*% u)^2)

print('The final prediction error on V is')
print( subdistance(sqy %*% Uhat, sqy %*% v)^2)
FIT_SAR_BIC<-SparseCCA(X=X,Y=Y,rank=k,
                       lambdaAseq=seq(from=20,to=0.02,length=20),
                       lambdaBseq=seq(from=20,to=0.02,length=20),
                       max.iter=100,conv=10^-2,
                       selection.criterion=1,n.cv=5)
print('The final prediction error on U is')
print( subdistance(sqx %*% FIT_SAR_BIC$ALPHA , sqx %*% u)^2)
#t(sqx %*% FIT_SAR_BIC$ALPHA)%*% sqx %*% FIT_SAR_BIC$ALPHA


Witten_Perm <- CCA.permute(x=X,z=Y,typex="standard",typez="standard", nperms=50)
WittenSCCA_Perm<-CCA(x=X_train, z=Y_train, typex="standard",typez="standard",K=rank,penaltyx=Witten_Perm$bestpenaltyx,penaltyz=Witten_Perm$bestpenaltyz,trace=F)

