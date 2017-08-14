#!/usr/bin/python

import sys
import zipfile
from datetime import datetime
import os

months = {'Jan':1, 'Feb':2, 'Mar':3, 'Apr':4, 'May':5, 'Jun':6, 'Jul':7,
          'Aug':8, 'Sep':9, 'Oct':10, 'Nov':11, 'Dec':12}
v2_date = datetime(2008, 9, 1)
v3_date = datetime(2011, 1, 1)
v4_date = datetime(2013, 11, 1)
year = ''
month = ''
day = ''
file_date = ''
filename = ''
versioin = ''
build = ''
hg = ''
v_detected = False
b_detected = False

raw_path = sys.argv[1]
output_dir = sys.argv[2]

# detect verison and reference build
with zipfile.ZipFile(raw_path, 'r') as zf:
    filename = zf.namelist()[0]
    
    with zf.open(filename, 'r') as f:
        line = f.readline()
        words = line.split()
        
        while words[0] == '#'and (not v_detected or not b_detected): 
            # this line contains meta-info and some info not detected yet
            if "generated" in words:
                year = int(words[-1])
                month = months[words[-4]]
                day = int(words[-3])
                file_date = datetime(year, month, day)
                if file_date < v2_date:
                    version = '1'
                elif file_date < v3_date:
                    version = '2'
                elif file_date < v4_date:
                    version = '3'
                else:
                    version = '4'                
                v_detected = True
            
            if "build" in words:
                i = words.index("build")
                build = words[i+1][:2]
                if build == '36':
                	hg = '18'
                elif build == '37':
                	hg = '19'
                b_detected = True
                
            # update
            line = f.readline()
            words = line.split()
         
        if not v_detected or not b_detected:
            # not more meta-info entry, but some info not detected yet
            print "This file does not contain version or build info"
            sys.exit()
        else: 
            print "[version, reference]: [%s, hg%s]" % (version, hg)

# convert
dir_name = filename[:-4]
output_dir = output_dir + dir_name

os.system("mkdir " + output_dir)
print "output_dir:", output_dir
os.system("perl 23andme2vcf_edited_3.pl "+raw_path+" "+output_dir+" "+version+" "+hg)

# liftover hg18ToHg19
if hg == '18':
	new_filename = dir_name + ".hg18ToHg19.vcf"
	os.system("cat "+output_dir+"/"+dir_name+".vcf "+\
		"| python ~/host/lift-over-vcf-master/lift_over.py --chain hg18ToHg19 "+\
		"> "+output_dir+"/"+new_filename)