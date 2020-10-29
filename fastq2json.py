#!/usr/bin/env python3

import json
import os
import re
from os.path import join
import argparse
from collections import defaultdict

'''
Modified from https://github.com/crazyhottommy/pyflow-RNAseq/blob/master/fastq2json.py

'''

parser = argparse.ArgumentParser()
parser.add_argument("--fastq_dir", nargs='+',
help="Required. the FULL path to the fastq folder(s)")
args = parser.parse_args()

assert args.fastq_dir is not None, "please provide the path to the fastq folder"


## default dictionary is quite useful!

FILES = defaultdict(lambda: defaultdict(list))

## build the dictionary with full path for each fastq.gz file
for folder in args.fastq_dir:
	for root, dirs, files in os.walk(folder):
		for f in files:
			if f.endswith("fastq.gz") or f.endswith("fq.gz") or f.endswith("fastq") or f.endswith("fq"):
				full_path = join(root, f)
				if '_00' in f:
					if 'L00' in f:
						m = re.search(r"(.+)_(L[0-9]{3})_(R[12])_00[0-9].(fastq.gz|fq.gz|fastq|fq)", f)
						reads_g = 3
					else:
						m = re.search(r"(.+)_(R[12])_00[0-9].(fastq.gz|fq.gz|fastq|fq)", f)
						reads_g = 2
				else:	
					if 'L00' in f:
						m = re.search(r"(.+)_(L[0-9]{3})_(R[12]).(fastq.gz|fq.gz|fastq|fq)", f)
						reads_g = 3
					else:
						m = re.search(r"(.+)_(R[12]).(fastq.gz|fq.gz|fastq|fq)", f)
						reads_g = 2
				#R1 will be forward reads, R2 will be reverse reads	
				if m:
					sample = m.group(1)
					reads = m.group(reads_g)  
					FILES[sample][reads].append(full_path)

#Make sure file from different lanes are in correct order
FILES_sorted = defaultdict(lambda: defaultdict(list))

for sample in FILES.keys():
		for read in FILES[sample]:
			FILES_sorted[sample][read] = sorted(FILES[sample][read])

				
print()
print ("total {} unique samples will be processed".format(len(FILES.keys())))
print ("------------------------------------------")
for sample in FILES.keys():
	for read in FILES[sample]:
		print ("{sample} {read} has {n} fastq".format(sample = sample, read = read, n = len(FILES[sample][read])))
print ("------------------------------------------")
print("check the samples.json file for fastqs belong to each sample")
print()

js = json.dumps(FILES_sorted, indent = 4, sort_keys=True)
open('samples.json', 'w').writelines(js)

