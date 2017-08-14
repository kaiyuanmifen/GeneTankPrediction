#!/usr/bin/python

import sys
import zipfile
# import os

version = ''
filename = ''
v_detected = False

raw_path = sys.argv[1]
# output_dir = sys.argv[2]

with zipfile.ZipFile(raw_path, 'r') as zf:
	filename = zf.namelist()[0]

	with zf.open(filename, 'r') as f:
		line = f.readline()
		words = line.split()

		while '#' in words[0] and not v_detected:
			# print words
			if "collected" in words:
				if words[-1] == 'V1.0':
					version = '1'
				else:
					version = '2'
				v_detected = True

			line = f.readline()
			words = line.split()

		if not v_detected:
			print "This file does not contain version info"
			sys.exit()
		else:
			print "version: V%s.0" % version