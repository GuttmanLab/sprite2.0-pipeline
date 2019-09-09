#!/bin/bash

mkdir -p slurms

snakemake \
--snakefile Snakefile \
--use-conda \
-j 32 \
--cluster "sbatch -t 24:00:00 --nodes 1 \
--ntasks 1 --cpus-per-task 10 --mem 100G \
--output slurms/job-%A_%a.out"

#--config bID="./config.txt" assembly="mm10" type="RNA-DNA" num_tags="5" \