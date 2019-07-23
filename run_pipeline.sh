#!/bin/bash



snakemake \
--config bID="./config.txt" \
--snakefile Snakefile_star.smk \
--use-conda \
-j 32 \
--cluster "sbatch -t 24:00:00 --nodes 1 --ntasks 1 --cpus-per-task 10 --mem 100G"
