#!/usr/bin/env python3

#########
# Author: Harvey Orel 
# Date: 28 Nov 2022
# Copies files to 
# Usage: 'alignment_locus_selector.py <./source_folder> </destination_folder> <gene_list.txt> <"keep" or "remove">'
#
# Note ::: Currently configured to work on files with the format 'LOCUS.any_suffix'
#########

import os
import re
import sys
import shutil

source = str(sys.argv[1])
destination = str(sys.argv[2])
gene_list = str(sys.argv[3])
method = str(sys.argv[4])

# Get suffix for files in source_directory
suffix = os.listdir(source)[0]
suffix = suffix[suffix.find('.'):]

# Run Main
if method == "keep":
	print("Keeping loci in the file '%s'" % gene_list)
	with open(gene_list, "r") as lines: # adapt the name of the file to open to your exact location.
	    filenames_to_copy = set(line.rstrip() for line in lines)

	for filename in filenames_to_copy:
	    source_path = os.path.join(source, filename)
	    source_path = source_path + suffix
	    if os.path.exists(source_path):
	        print("copying {} to {}".format(source_path, destination))
	        shutil.copy(source_path, destination)

if method == "remove":
	print("Removing loci in the file '%s'" % gene_list)
	with open(gene_list, "r") as lines: # adapt the name of the file to open to your exact location.
	    filenames_not_to_copy = set(line.rstrip() for line in lines)

	# Loop through all files in the source folder and copy them to the destination folder
	for filename in os.listdir(source):
	    src_file = os.path.join(source, filename)
	    dst_file = os.path.join(destination, filename)
	    shutil.copy2(src_file, dst_file)

	# Loop through the list of filenames to be removed and delete them from the destination folder
	for filename in filenames_not_to_copy:
	    file_path = os.path.join(destination, filename)
	    file_path = file_path + suffix
	    print("removing %s" % file_path)
	    os.remove(file_path)