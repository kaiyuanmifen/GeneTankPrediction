#!/usr/bin/python

import sys
import os

file_path = sys.argv[1]
l = file_path.split('/')
filename = l[len(l)-1]
chrN = filename.split('.')[0]
len_filename = len(filename)
file_dir = file_path[:-len_filename]
site_path = file_dir+chrN+".duplicate.site"

os.system("touch "+site_path)

r = open(file_path, "r")
w = open(site_path, "w")

# skip mete-info
line = r.readline()
while line[0] == '#':
	line = r.readline()

# now at 1st entry
pre_pos = ''

# read line by line until EOF, write all duplicate positions
while line != "":
	cur_pos = line.split('\t')[1]
	
	if cur_pos == pre_pos:
		w.write(cur_pos+"\n")

	pre_pos = cur_pos
	line = r.readline()

r.close()
w.close()