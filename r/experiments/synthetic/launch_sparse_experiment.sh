#!/bin/bash
  
# Define the values for the variables
n_values="300 500 1000 10000"
r_values="2 5"
rpca_values="0 5 10"
ratio_values="0.3 0.5 1 1.5 2"

for n in $n_values; do
  for r in $r_values; do
    for r_pca in $rpca_values; do
       for ratio in $ratio_values; do
          sbatch experiment/synthetic/sparse_experiment.sh "$n" "$r" "$r_pca" prediction 1 "$ratio"
        done
    done
  done
done
# $1 : N
# $2 : r
# $3 : r_pcas
# $4 : criterion (prediction/ correlation) for CV
# $5 : normalized diagonal (0/1)
# $6 : ratio p/n