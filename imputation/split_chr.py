#!/usr/bin/python

import sys
import os

file_path = sys.argv[1]
l = file_path.split('/')
filename = l[len(l)-1]
len_filename = len(filename)
ext = filename[-3:]
file_dir = file_path[:-len_filename] + "unphased/"
# print file_dir

os.system("mkdir "+file_dir)

# copy paste meta-info from original VCF to sub-VCFs
os.system("for chr in $(seq 1 22); \
	do cat "+file_path+" \
	| awk '{ if (length($1) > 2) { print $0 } }' \
	> "+file_dir+"chr$chr\.unphased."+ext+"; \
	done")

# split original VCF based on chromosome, append to existing files
os.system("for chr in $(seq 1 22); \
	do cat "+file_path+" \
	| awk -v c=$chr '{ if ($1 == c) { print $0 } }' \
	>> "+file_dir+"chr$chr\.unphased."+ext+"; \
	done")