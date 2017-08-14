import sys
import os

i = sys.argv[1]
unphased_dir = sys.argv[2]
map_dir = sys.argv[3]
ref_dir = sys.argv[4]

log_dir = unphased_dir+"log/"
os.system("mkdir "+log_dir)

gwas= unphased_dir+"chr"+str(i)+".unphased.qc0.vcf"
genetic_map = map_dir+"genetic_map_chr"+str(i)+"_combined_b37.txt"
duplicate = unphased_dir+"chr"+str(i)+".duplicate.site"
ref_hap = ref_dir+"1000GP_Phase3_chr"+str(i)+".hap.gz"
ref_legend = ref_dir+"1000GP_Phase3_chr"+str(i)+".legend.gz"
ref_sample = ref_dir+"1000GP_Phase3.sample"
alignments = log_dir+"chr"+str(i)+".unphased.alignments"

allExist = True
		
if not os.path.isfile(gwas):
	allExist = False
	print gwas+" does not exist"
	
if not os.path.isfile(genetic_map):
	allExist = False
	print genetic_map+" does not exist"

if not os.path.isfile(duplicate):
	allExist = False
	print duplicate+" does not exist"

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
	os.system("shapeit -check -V "+gwas+\
		" -M "+genetic_map+\
		" --exclude-snp "+duplicate+\
		" --input-ref "+ref_hap+" "+ref_legend+" "+ref_sample+\
		" --output-log "+alignments)