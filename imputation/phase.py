#!/usr/bin/python

import sys
import os

i = sys.argv[1]
unphased_dir = sys.argv[2]
map_dir = sys.argv[3]
ref_dir = sys.argv[4]

l = unphased_dir.split('/')
dir_name = l[len(l)-2]
phased_dir = unphased_dir[:-len(dir_name)-1] + "phased/"

os.system("mkdir "+phased_dir)

gwas= unphased_dir+"chr"+str(i)+".unphased.qc0.vcf"
genetic_map = map_dir+"genetic_map_chr"+str(i)+"_combined_b37.txt"
exclude = unphased_dir+"chr"+str(i)+".exclude.site"
ref_hap = ref_dir+"1000GP_Phase3_chr"+str(i)+".hap.gz"
ref_legend = ref_dir+"1000GP_Phase3_chr"+str(i)+".legend.gz"
ref_sample = ref_dir+"1000GP_Phase3.sample"

# modify output from quality control
output = phased_dir+"chr"+str(i)+".phased.with.ref"
log = phased_dir+"chr"+str(i)+".phased.with.ref"

allExist = True
		
if not os.path.isfile(gwas):
	allExist = False
	print gwas+" does not exist"
	
if not os.path.isfile(genetic_map):
	allExist = False
	print genetic_map+" does not exist"

if not os.path.isfile(exclude):
	allExist = False
	print exclude+" does not exist"

if not os.path.isfile(ref_hap):
	allExist = False
	print ref_hap+" does not exist"

if not os.path.isfile(ref_legend):
	allExist = False
	print ref_legend+" does not exist"

if not os.path.isfile(ref_sample):
	allExist = False
	print ref_sample+" does not exist"

if allExist:
	os.system("shapeit -V "+gwas+\
		" -M "+genetic_map+\
		" --exclude-snp "+exclude+\
		" --input-ref "+ref_hap+" "+ref_legend+" "+ref_sample+\
		" -O "+output+\
		" --output-log "+log)