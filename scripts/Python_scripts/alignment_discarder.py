#!/usr/bin/env python3

###############################
# Harvey Orel
# Date: 24 Oct 2022
# Modified: 4 Nov 2022
#			- Included an additional argument for minimum number of samples required to keep an alignment.
# Description: Filters a folder of alignments and discards those with less than n samples
# Note: Input is a directory of fasta alignments. Output is to "final_alignments/" in launch directory;
#
#       Arguments are: <input_alignment_directory> <min_number_samples>
###############################


import sys
import os
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO

directory = str(os.getcwd())

input_directory = str(sys.argv[1])
min_number_samples = int(sys.argv[2])

# Create output directory
if not os.path.exists(directory + "/final_alignments/"):
	os.makedirs(directory + "/final_alignments/")
else:
	print("Error: directory 'final_alignments' already exists... Cancelling")
	sys.exit(1)


removed_alignment_counter = 0

# Iterate through files and remove samples in samples_2_remove.txt, write to new fasta files
for filename in os.listdir(input_directory):
	if filename.endswith(".fasta") or filename.endswith(".fas"):
		new_filename = filename.partition(".fasta")[0]
		print("Testing alignment: " + filename)
		input_filename = open(input_directory + filename, 'r')
		output_directory = str(directory + "/final_alignments/")

		alignment = SeqIO.parse(input_filename, 'fasta')

		sample_names = []
		for seqs in alignment:
			sample_name = seqs.id
			sample_names.append(sample_names)
		
		if len(sample_names) < min_number_samples:
			removed = ('%s was removed.' % filename)
			indicator = "***"
			print(f"{removed:<70}{indicator:>20}")
			removed_alignment_counter += 1
		else:
			print(str(len(sample_names)) + " samples in alignment")
			shutil.copy(str(input_directory + filename), output_directory)

		input_filename.close()

print("\nDiscarded %s alignments" % removed_alignment_counter)

print("DONE!!!")

