library(tidyverse)

file_list <- list.files(path = "~/Documents/CCA-experiments/r/experiments/synthetic/results/", 
                        pattern = "new_exp_9*", full.names = TRUE)
# Read and combine CSV files into a single data frame
results <- file_list %>%
  map_dfr(~ read_csv(.x) %>% mutate(filename = .x)) %>%
  mutate("selection" = ifelse(str_detect(filename, "correlation"), "correlation", "prediction"))

res = results %>% 
  mutate(FDR=1 - nb_real_discoveries/max(1,nb_discoveries)) %>%
  group_by(method, selection, criterion, n, normalize_diagonal,
           nnz, p1, p2, r, r_pca,
           overlapping_amount) %>% 
  summarise_if(is.numeric, median) %>%
  arrange(n, nnz, p1, p2, distance_tot)
t = as.numeric(sapply(res$method, function(x){ifelse(str_count(x, "folds" )>0, 1, 0)}))
t = t  + as.numeric(sapply(res$method, function(x){ifelse(str_count(x, "CV" )>0, 2, 0)}))
res["shape"] = as.factor(t)


legend_order <- c( "SSVD-theory" ,  "SSVD-method","SSVD-theory-init-norm",
                   "SSVD-method-init-norm",
                   "TG" , "Fantope" , "Oracle",  "cancor")

my_colors <- c( "gray", "black",
                "dodgerblue", "cyan",  "red","orange","orange4","navy", "grey", "yellow",  
                "green", "green3", "navajowhite3", "plum2", "plum3",
                "cadetblue1", "lightskyblue", "brown", "beige","whitesmoke", "white")

