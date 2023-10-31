seed =1
r=3
name_exp <- "test.csv"
n <- 300
r_pca = 3
criterion = "prediction"
normalize_diagonal <- T
ratio <- 1.5
set.seed(seed)
psize = ratio * n
nnz = 12
it = 1
overlapping_amount = 0
THRES = 1e-3
signal_strength  = "medium"


alpha.method = 0.05;
alpha.theory = 1.5; huber.beta = 0.95; sigma = NA;
gamma.u = sqrt(2); gamma.v = sqrt(2);  dothres = "hard"; tol = 1e-08;
n.iter = 100; n.boot = 100; non.orth = FALSE; reps = 1; init_norm = FALSE

