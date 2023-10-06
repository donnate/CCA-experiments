library(MASS)
library(geigen)
library(pracma)
library(VGAM)
library(mvtnorm)
n= 500
p1=210
p2=230
source("experiments/synthetic/generate_synthetic_examples.R")

example <- generate_example_none_trivial_pca( n, p1, p2, nnzeros = 10,
                                              theta = diag( c(0.8,  0.9)), 
                                              overlapping_amount = 0,
                                              r=2, r_pca = 3)

example <- generate_simple_example( n, p1, p2, nnzeros = 10,
                                              theta = diag( c(0.8,  0.9)), 
                                              r=2)
