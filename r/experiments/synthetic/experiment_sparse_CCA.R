wd = getwd()
#wd = "~/Documents/CCA-experiments/"
print(wd)
setwd(wd)
source("r/experiments/synthetic/experiment_functions.R")
source("r/experiments/synthetic/generate_examples.R")
source("r/experiments/synthetic/evaluate_results.r")
source("r/experiments/synthetic/generate_examples.R")
source("r/ssvd/my_ssvd.R")
results <- c()

args <- commandArgs(trailingOnly=TRUE)
seed <- as.integer(as.numeric(args[1]))
#seed=1
print(seed)
name_exp <- args[2]
n <- as.integer(as.numeric(args[3]))
r <- ceiling(as.numeric(args[4]))
r_pca <- ceiling(as.numeric(args[5]))
criterion <- args[6]
normalize_diagonal <- as.logical(as.numeric(args[7]))
ratio <- as.numeric(args[8])
set.seed(seed)
it = seed
THRES = 1e-3
print(args)

#sparsity = 0.05


for (psize in c(ratio * n)) {
  for (sparsity in c(0.01, 0.05, 0.1, 0.2)){
      nnz = ceil(sparsity * psize)
      print(c(max(r_pca, r) * nnz, psize, max(r_pca, r), nnz))
      if ((max(r_pca, r) * nnz < psize) && (max(r_pca, r) < nnz)) { 
        if (nnz > 2){
          for(signal_strength in c("high", "medium", "low")){
            if ( signal_strength == "high"){
                thetas <- diag(seq(0.95, 0.8, length.out = r))
              }else{
                if ( signal_strength == "medium"){
                  thetas <- diag(seq(0.65, 0.55, length.out = r))
                }
                else{
                  thetas <- diag(seq(0.5, 0.35, length.out = r))
                }
              }
          for (overlapping_amount in c(0, 1)){
            p1 <- as.integer(psize); 
            p2 <- as.integer(psize)
            print(c(n, p1, p2, nnz))
            p <-  p1 + p2
            example <- generate_example_sparse_product(n, p1, p2, 
                                                r_pca = r_pca, 
                                                nnzeros = nnz,
                                                theta = thetas,
                                                lambda_pca = 1,
                                                r=r, 
                                                overlapping_amount = overlapping_amount,
                                                normalize_diagonal=normalize_diagonal)
            print("here")
            max1 = 500 * sqrt(log(p1)/n)
            min1 = 0.001 * sqrt(log(p1)/n) 
            param1 = exp(seq(log(min1), log(max1), length.out=20))
            maxk = 0.25 * p
            mink = 0.01 * p 
            param2 = ceiling(seq(max(ceiling(mink),5), ceiling(maxk), length.out = 10))
            fantope_solution = NULL
            
            ssvd_results <- tryCatch({
              test1<-my.ssvd((example$S[1:p1, (p1+1):p]),
                            Sigma_u = example$S[1:p1, 1:p1],
                            Sigma_v=example$S[(p1+1):p, (p1+1):p],
                            r=r, method = "theory", reps=1)
              Uhat = rbind(test1$u, test1$v)
              temp <- evaluate_results(Uhat= test1$u, 
                                      Vhat = test1$v, 
                                      example = example, 
                                      name_method="SSVD-theory", 
                                      overlapping_amount=overlapping_amount,
                                      lambdax= NA,
                                      lambday = NA, 
                                      thres = THRES,
                                      it=it,
                                      normalize_diagonal=normalize_diagonal,
                                      criterion=criterion,
                                      r_pca = r_pca, nnz=nnz,
                                      signal_strength=signal_strength)
              if (length(results)==0){
                results=temp
              }else{
                results <- rbind(results, temp )
               
              }
              #### Evaluate goodness of initial fit
              temp2 <- evaluate_results(Uhat= test1$u_init, 
                                        Vhat = test1$v_init, 
                                        example = example, 
                                        name_method="SSVD-theory-just-init", 
                                        overlapping_amount=overlapping_amount,
                                        lambdax= NA,
                                        lambday = NA, 
                                        thres = THRES,
                                        it=it,
                                        normalize_diagonal=normalize_diagonal,
                                        criterion=criterion,
                                        r_pca = r_pca, nnz=nnz,
                                        signal_strength=signal_strength)
              results <- rbind(results, temp2)

            }, error = function(e) {
              # Print the error message
              cat("Error occurred in method sparse SSVD (theory)", ":", conditionMessage(e), "\n")
              # Skip to the next iteration
            })
            
            #
            ssvd_results <- tryCatch({
              test1<-my.ssvd((example$S[1:p1, (p1+1):p]),
                             Sigma_u = example$S[1:p1, 1:p1],
                             Sigma_v=example$S[(p1+1):p, (p1+1):p],
                             x_full = example$S,
                             r=r, method = "theory", reps=1,
                             init_type = "full")
              Uhat = rbind(test1$u, test1$v)
              temp <- evaluate_results(Uhat= test1$u, 
                                       Vhat = test1$v, 
                                       example = example, 
                                       name_method="SSVD-row_selection_full_x", 
                                       overlapping_amount=overlapping_amount,
                                       lambdax= NA,
                                       lambday = NA, 
                                       thres = THRES,
                                       it=it,
                                       normalize_diagonal=normalize_diagonal,
                                       criterion=criterion,
                                       r_pca = r_pca, nnz=nnz,
                                       signal_strength=signal_strength)
              if (length(results)==0){
                results=temp
              }else{
                results <- rbind(results, temp )
                
              }
              #### Evaluate goodness of initial fit
              temp2 <- evaluate_results(Uhat= test1$u_init, 
                                        Vhat = test1$v_init, 
                                        example = example, 
                                        name_method="SSVD-row_selection_full_x-just-init", 
                                        overlapping_amount=overlapping_amount,
                                        lambdax= NA,
                                        lambday = NA, 
                                        thres = THRES,
                                        it=it,
                                        normalize_diagonal=normalize_diagonal,
                                        criterion=criterion,
                                        r_pca = r_pca, nnz=nnz,
                                        signal_strength=signal_strength)
              results <- rbind(results, temp2)
              
            }, error = function(e) {
              # Print the error message
              cat("Error occurred in method sparse SSVD (theory)", ":", conditionMessage(e), "\n")
              # Skip to the next iteration
            })
            
            ssvd_results <- tryCatch({
              test1<-my.ssvd((example$S[1:p1, (p1+1):p]),
                             Sigma_u = example$S[1:p1, 1:p1],
                             Sigma_v=example$S[(p1+1):p, (p1+1):p],
                             x_full = NULL,
                             r=r, method = "theory", reps=1,
                             init_type = "full")
              Uhat = rbind(test1$u, test1$v)
              temp <- evaluate_results(Uhat= test1$u, 
                                       Vhat = test1$v, 
                                       example = example, 
                                       name_method="SSVD-top-rows-original", 
                                       overlapping_amount=overlapping_amount,
                                       lambdax= NA,
                                       lambday = NA, 
                                       thres = THRES,
                                       it=it,
                                       normalize_diagonal=normalize_diagonal,
                                       criterion=criterion,
                                       r_pca = r_pca, nnz=nnz,
                                       signal_strength=signal_strength)
              if (length(results)==0){
                results=temp
              }else{
                results <- rbind(results, temp )
                
              }
              #### Evaluate goodness of initial fit
              temp2 <- evaluate_results(Uhat= test1$u_init, 
                                        Vhat = test1$v_init, 
                                        example = example, 
                                        name_method="top-rows-original-just-init", 
                                        overlapping_amount=overlapping_amount,
                                        lambdax= NA,
                                        lambday = NA, 
                                        thres = THRES,
                                        it=it,
                                        normalize_diagonal=normalize_diagonal,
                                        criterion=criterion,
                                        r_pca = r_pca, nnz=nnz,
                                        signal_strength=signal_strength)
              results <- rbind(results, temp2)
              
            }, error = function(e) {
              # Print the error message
              cat("Error occurred in method sparse SSVD (theory)", ":", conditionMessage(e), "\n")
              # Skip to the next iteration
            })
          
            
      
          
            result <- tryCatch({
              normalize=FALSE
              name_method = paste0("TG", ifelse(normalize, "_normalized", ""))
              res_tg <- pipeline_thresholded_gradient(example$Data, example$Mask, 
                                                      example$sigma0hat, 
                                                      r=r, nu=1,
                                                      Sigmax=example$Sigmax, 
                                                      Sigmay=example$Sigmay, 
                                                      maxiter.init=100, 
                                                      lambda=NULL,k=NULL,
                                                      kfolds=5, maxiter=2000, 
                                                      convergence=1e-3, eta=1e-3,
                                                      param1=param1,
                                                      param2=param2, normalize=normalize,
                                                      criterion=criterion,
                                                      fantope_solution=fantope_solution)
              Uhat = rbind(res_tg$ufinal, res_tg$vfinal)
              temp <- evaluate_results(Uhat= res_tg$ufinal, 
                                       Vhat = res_tg$vfinal, 
                                       example = example, 
                                       name_method=name_method, 
                                       overlapping_amount=overlapping_amount,
                                       lambdax=  res_tg$lambda,
                                       lambday =   res_tg$k, 
                                       thres = THRES,
                                       it=it,
                                       normalize_diagonal=normalize_diagonal,
                                       criterion=criterion,
                                       r_pca = r_pca, nnz=nnz,
                                       signal_strength=signal_strength)
              
              results <- rbind(results, temp)
              
              #Uhat = rbind(res_tg$initu, res_tg$initv)
              temp <- evaluate_results(Uhat= res_tg$initu, 
                                       Vhat = res_tg$initv, 
                                       example = example, 
                                       name_method="Fantope", 
                                       overlapping_amount=overlapping_amount,
                                       lambdax=  res_tg$lambda,
                                       lambday =   res_tg$k, 
                                       thres = THRES,
                                       it=it,
                                       normalize_diagonal=normalize_diagonal,
                                       criterion=criterion,
                                       r_pca = r_pca, nnz=nnz,
                                       signal_strength=signal_strength)
              results <- rbind(results, temp)
              selected_rows = which(apply(res_tg$initu^2, 1, sum)>0)
              selected_rows.v = which(apply(res_tg$initv^2, 1, sum)>0)
              print("Selected rows.v")
              print(selected_rows.v)
              
              print("Selected rows")
              print(selected_rows)
              #### Try one SSVD with Fantope init
              test <- my.ssvd.iter.thresh(example$S[1:p1, (p1+1):p], 
                                          Sigma_u = example$S[1:p1, 1:p1],
                                          Sigma_v=example$S[(p1+1):p, (p1+1):p],
                                          method = "theory", u.old = res_tg$initu,
                                          v.old = res_tg$initv, r = r,
                                          gamma.u = sqrt(2), gamma.v = sqrt(2), dothres = "hard", tol = 1e-08,
                                          n.iter = 100, n.boot = 100, non.orth = FALSE, reps = 1)
              temp <- evaluate_results(Uhat= test$u, 
                                       Vhat = test$v, 
                                       example = example, 
                                       name_method="thresholding-Fantope-init", 
                                       overlapping_amount=overlapping_amount,
                                       lambdax=  res_tg$lambda,
                                       lambday =   res_tg$k, 
                                       thres = THRES,
                                       it=it,
                                       normalize_diagonal=normalize_diagonal,
                                       criterion=criterion,
                                       r_pca = r_pca, nnz=nnz,
                                       signal_strength=signal_strength)
              results <- rbind(results, temp)
              
              
              
              if (min(length(selected_rows), length(selected_rows.v))<n & min(length(selected_rows), length(selected_rows.v))>0){
                print("here CCA")
                t = cancor(as.matrix(example$Data[,selected_rows]), as.matrix(example$Data[, (selected_rows.v + p1)]))
                Uhat = matrix(0, p, r)
                Uhat[selected_rows,]  = t$xcoef[,1:r]
                Uhat[selected_rows.v + p1,]  = t$ycoef[,1:r]
                print(paste0("Done with Uhat: dimension:"))
                
                temp <- evaluate_results(Uhat= Uhat[1:p1,],
                                         Vhat = Uhat[(p1+1):p,], 
                                         example = example, 
                                         name_method="cancor", 
                                         overlapping_amount=overlapping_amount,
                                         lambdax= res_tg$lambda,
                                         lambday =  res_tg$k, 
                                         thres = THRES,
                                         it=it,
                                         normalize_diagonal=normalize_diagonal,
                                         criterion=criterion,
                                         r_pca = r_pca, nnz=nnz,
                                         signal_strength=signal_strength)
                
                print(dim(temp))
                print(dim(results))
                print("done")
              }else{
                print("Couldnt threshold")
                
              }
              results <- rbind(results, temp )
              
              
            }, error = function(e) {
              # Error handling code goes here
              
              # Print the error message
              cat("Error occurred in method", name_method, ":", conditionMessage(e), "\n")
              
              # Skip to the next iteration
              
            })
            
            
         
            
          #### Oracle
          set_u =  which(apply(example$u,1, norm)>0)
          set_v = nrow(example$u) + which(apply(example$v,1, norm)>0)
          t=CCA::cc(as.matrix(example$Data[,set_u]), as.matrix(example$Data[, set_v]))
          Uhat = matrix(0, p, r)
          Uhat[set_u, ] =  t$xcoef[,1:r]
          Uhat[set_v, ] =  t$ycoef[,1:r]
          temp <- evaluate_results(Uhat= Uhat[1:p1,],
                                  Vhat = Uhat[(p1+1):p,], 
                                  example = example, 
                                  name_method = "Oracle", 
                                  overlapping_amount=overlapping_amount,
                                  lambdax= NA,
                                  lambday =  NA, 
                                  thres = THRES,
                                  it=it,
                                  normalize_diagonal=normalize_diagonal,
                                  criterion=criterion,
                                  r_pca = r_pca, nnz=nnz,
                                  signal_strength=signal_strength)
          results <- rbind(results, temp )
        write_excel_csv(results, paste0("r/experiments/synthetic/results/", name_exp, "_", criterion, ".csv"))
        }
      }
        }
      }
  }
}




