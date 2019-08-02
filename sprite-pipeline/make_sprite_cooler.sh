#!/bin/bash

# $1 : input file
# $2 : chrom.sizes
# $3 : cluster size lower bound
# $4 : cluster size upper bound

set -euxo pipefail

INPUT_CLUSTERS=$1
CHROMSIZES=$2
OPREFIX=${1%.clusters.gz}

zcat $INPUT_CLUSTERS \
    | cut -f 2- \
    | awk "NF >= ${3} && NF <= ${4}" \
    | python <(cat <<EOF
from itertools import combinations
import fileinput
for line in fileinput.input():

    coords = [coord.replace(':', '\t') for coord in line.strip().split('\t')]
    pairs = list(combinations(coords, 2))
    score = 2.0 / len(coords)
    for a, b in pairs:
        print(a, b, score, sep='\t')
EOF
)  \
   | cooler cload pairs \
        -c1 1 -p1 2 -c2 3 -p2 4 \
        --zero-based --chunksize 10000000 --field count=5:dtype=float32 \
        $CHROMSIZES:1000 - "${OPREFIX}_${3}-${4}.norm_n2.1000.cool"
