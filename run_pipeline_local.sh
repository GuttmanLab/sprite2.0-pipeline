#!/bin/bash

snakemake \
--snakefile Snakefile \
--use-conda \
--cores 10 \
--configfile config_local.yaml
