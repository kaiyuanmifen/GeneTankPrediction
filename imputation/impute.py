#!/usr/bin/python

import sys
import os

i = sys.argv[1]
lower = sys.argv[2]
upper = sys.argv[3]
phased_dir = sys.argv[4]
map_dir = sys.argv[5]
ref_dir = sys.argv[6] # consider merge

l = phased_dir.split('/')
dir_name = l[len(l)-2]
imputed_dir = phased_dir[:-len(dir_name)-1] + "imputed/"

os.system("mkdir "+imputed_dir)

# input
phased = phased_dir+"chr"+str(i)+".phased.with.ref.haps"
genetic_map = map_dir+"genetic_map_chr"+str(i)+"_combined_b37.txt"
ref_hap = ref_dir+"1000GP_Phase3_chr"+str(i)+".hap.gz"
ref_legend = ref_dir+"1000GP_Phase3_chr"+str(i)+".legend.gz"

# output
output = imputed_dir+"chr"+str(i)+".pos"+lower+"_to_"+upper+".imputed"

allExist = True
		
if not os.path.isfile(phased):
	allExist = False
	print phased+" does not exist"
	
if not os.path.isfile(genetic_map):
	allExist = False
	print genetic_map+" does not exist"

if not os.path.isfile(ref_hap):
	allExist = False
	print ref_hap+" does not exist"

if not os.path.isfile(ref_legend):
	allExist = False
	print ref_legend+" does not exist"

if allExist:
	os.system("impute2 -use_prephased_g"+\
        " -known_haps_g "+phased+\
        " -m "+genetic_map+\
        " -h "+ref_hap+\
        " -l "+ref_legend+\
        " -int "+lower+" "+upper+\
        " -Ne 20000"+\
        " -o "+output)