library(tidyverse)

file_list <- list.files(path = "~/Documents/CCA-experiments/r/experiments/synthetic/results/", 
                        pattern = "new_signal_strength*", full.names = TRUE)
# Read and combine CSV files into a single data frame
results <- file_list %>%
  map_dfr(~ read_csv(.x) %>% mutate(filename = .x)) %>%
  mutate("selection" = ifelse(str_detect(filename, "correlation"), "correlation", "prediction"))

res = results %>% 
  mutate(FDR=1 - nb_real_discoveries/max(1,nb_discoveries)) %>%
  group_by(method, selection, criterion, n, normalize_diagonal,
           signal_strength,
           nnz, p1, p2, r, r_pca,
           sparsity,
           overlapping_amount) %>% 
  summarise_if(is.numeric, mean) %>%
  arrange(n, nnz, p1, p2, distance_tot)
t = as.numeric(sapply(res$method, function(x){ifelse(str_count(x, "folds" )>0, 1, 0)}))
t = t  + as.numeric(sapply(res$method, function(x){ifelse(str_count(x, "CV" )>0, 2, 0)}))
res["shape"] = as.factor(t)


unique(res$method)
legend_order <- c( "SSVD-theory" ,  "SSVD-method", "SSVD-theory-init-norm",
                   "SSVD-method-init-norm", "SSVD-theory-just-init", "SSVD-method-just-init",
                   "SSVD-top-rows-original", "top-rows-original-just-init"  ,
                   "SSVD-row_selection_full_x", "SSVD-row_selection_full_x-just-init",
                   "thresholding-Fantope-init" ,
                   "TG" , "Fantope" , "cancor",  "Oracle")

my_colors <- c( 
                "dodgerblue", "navy", "cyan", "lightblue3",  
                 "orange", 
                "yellow",
                "navajowhite1",  "navajowhite3",
                "purple1", "orchid4",
                "seagreen4",
                "red", "pink1", #"green1", "green4", 
                "beige", "black")

unique(res$r)
unique(res$r_pca)
unique(res$n)
unique(res$p1)
unique(res$overlapping_amount)
ggplot(res %>% filter(r == 3 , r_pca ==0,  overlapping_amount==0,
                      n ==700),
       aes(x=nnz/p1, y=distance_tot, colour=method))+
  #geom_line(aes(x=nnz/p1, y= distance_tot, colour=method), linewidth=1.)+
  #geom_point()+
  geom_jitter(width = 0.015, height = 0., size=2)+
  #geom_jitter(aes(x=p1/n, y=distance_tot, colour=method, shape=shape), size=2.2)+
  #geom_line(data=res %>% group_by(n, nnz, p1, p2) %>% summarise(b=mean(zero_benchmark)),
  #          aes(y=b, x=p1/n, colour="Zero Benchmark"), colour="black", linewidth=0.5)+
  scale_color_manual(values = my_colors, breaks = legend_order) +
  facet_grid(p1/n ~ signal_strength, scales="free") +
  theme_bw()

unique(res$r)
ggplot(res %>% filter(r == 5 , r_pca ==5,  overlapping_amount==0,
                      n ==300),
       aes(x=nnz/p1, y=distance_tot, colour=method))+
  #geom_line(aes(x=nnz/p1, y= distance_tot, colour=method), linewidth=1.)+
  #geom_point()+
  geom_jitter(width = 0.02, height = 0., size=2)+
  #geom_jitter(aes(x=p1/n, y=distance_tot, colour=method, shape=shape), size=2.2)+
  #geom_line(data=res %>% group_by(n, nnz, p1, p2) %>% summarise(b=mean(zero_benchmark)),
  #          aes(y=b, x=p1/n, colour="Zero Benchmark"), colour="black", linewidth=0.5)+
  scale_color_manual(values = my_colors, breaks = legend_order) +
  facet_grid(p1/n ~ signal_strength, scales="free") +
  theme_bw()
 ggplot(res %>% filter(r == 5, r_pca ==5, 
                      n == 500) )+
  geom_line(aes(x=nnz/p1, y=1/sqrt(p1) * distance_tot, colour=method), linewidth=1.)+
  geom_point(aes(x=nnz/p1, y=1/sqrt(p1) * distance_tot, colour=method), size=2.2)+
  #geom_jitter(aes(x=p1/n, y=distance_tot, colour=method, shape=shape), size=2.2)+
  #geom_line(data=res %>% group_by(n, nnz, p1, p2) %>% summarise(b=mean(zero_benchmark)),
  #          aes(y=b, x=p1/n, colour="Zero Benchmark"), colour="black", linewidth=0.5)+
  scale_color_manual(values = my_colors, breaks = legend_order) +
  facet_grid(n/p1~overlapping_amount, scales="free") +
  theme_bw()

unique(res$sparsity)
ggplot(res %>% filter(r == 5, r_pca ==0, 
                      n == 300,
                      overlapping_amount ==0) )+
  geom_line(aes(x=nnz, y=distance_tot, colour=method), linewidth=1.)+
  geom_point(aes(x=nnz, y=distance_tot, colour=method), size=2.2)+
  #geom_jitter(aes(x=p1/n, y=distance_tot, colour=method, shape=shape), size=2.2)+
  #geom_line(data=res %>% group_by(n, nnz, p1, p2) %>% summarise(b=mean(zero_benchmark)),
  #          aes(y=b, x=p1/n, colour="Zero Benchmark"), colour="black", linewidth=0.5)+
  scale_color_manual(values = my_colors, breaks = legend_order) +
  facet_grid(n~p1, scales="free") +
  theme_bw()

ggplot(res %>% filter(r == 5, r_pca ==5, 
                      n == 300,
                      overlapping_amount ==0) )+
  geom_line(aes(x=nnz, y=distance_tot, colour=method), linewidth=1.)+
  geom_point(aes(x=nnz, y=distance_tot, colour=method), size=2.2)+
  #geom_jitter(aes(x=p1/n, y=distance_tot, colour=method, shape=shape), size=2.2)+
  #geom_line(data=res %>% group_by(n, nnz, p1, p2) %>% summarise(b=mean(zero_benchmark)),
  #          aes(y=b, x=p1/n, colour="Zero Benchmark"), colour="black", linewidth=0.5)+
  scale_color_manual(values = my_colors, breaks = legend_order) +
  facet_grid(n~p1, scales="free") +
  theme_bw()

unique(res$n)
ggplot(res %>% filter(r == 2, r_pca ==0, overlapping_amount == 1) )+
  geom_line(aes(x=nnz/p1, y=distance_tot, colour=method), linewidth=1.)+
  geom_point(aes(x=nnz/p1, y=distance_tot, colour=method), size=2.2)+
  #geom_jitter(aes(x=p1/n, y=distance_tot, colour=method, shape=shape), size=2.2)+
  #geom_line(data=res %>% group_by(n, nnz, p1, p2) %>% summarise(b=mean(zero_benchmark)),
  #          aes(y=b, x=p1/n, colour="Zero Benchmark"), colour="black", linewidth=0.5)+
  scale_color_manual(values = my_colors, breaks = legend_order) +
  facet_grid(n~p1, scales="free") +
  theme_bw()
