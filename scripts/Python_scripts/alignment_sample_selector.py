#!/usr/bin/env python3

###############################
# Harvey Orel
# Date: 24 Oct 2022
# Description: Removes samples from alignments based on a list in a text file (can handle tab, new-line or comma delimited text)
# Note: Input is a directory of fasta alignments. Output is to "final_alignments/"" in launch directory;
#
#       Arguments are: <input_alignment_directory> <samples_2_remove.txt> <"keep" or "remove">
###############################


import sys
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO

directory = str(os.getcwd())

input_directory = str(sys.argv[1])
samples_2_remove = str(sys.argv[2])
method = str(sys.argv[3])


# Logic for arguments
if samples_2_remove is None:
	print("Error: no list of samples to remove... Exiting")
	sys.exit(1)
else:
	pass

# Logic for arguments
if method is None:
	print("Error: no preferred method ('keep' or 'remove'... Exiting")
	sys.exit(1)
else:
	pass

# Format sample list
with open(samples_2_remove) as f:
    samples_2_remove = f.read()

samples_2_remove = re.split(';|,|\n|\t', samples_2_remove)
print(samples_2_remove)

# Create output directory
if not os.path.exists(directory + "/final_alignments/"):
	os.makedirs(directory + "/final_alignments/")
else:
	print("Error: directory 'final_alignments' already exists... Cancelling")
	sys.exit(1)


# Iterate through files and remove samples in samples_2_remove.txt, write to new fasta files
for filename in os.listdir(input_directory):
	if filename.endswith(".fasta") or filename.endswith(".fas"):
		new_filename = filename.partition(".fasta")[0]
		print("Removing samples from: " + filename)
		input_filename = open(input_directory + filename, 'r')
		output_filename = open(directory + "/final_alignments/" + new_filename + ".fasta", 'w')

		if method == "remove":
			for seqs in SeqIO.parse(input_filename, 'fasta'):
				sample_name = seqs.id.split('.')[0]
				if re.search(rf"{sample_name}", '|'.join(samples_2_remove), re.IGNORECASE): # Update to accomodate .paralog numbers from HybPiper
					print('%s was removed.' % sample_name)
				else:
					SeqIO.write(seqs, output_filename, 'fasta')
					print('%s was kept.' % sample_name)
		elif method == "keep":
			for seqs in SeqIO.parse(input_filename, 'fasta'):
				sample_name = seqs.id.split('.')[0]
				if re.search(rf"{sample_name}", '|'.join(samples_2_remove), re.IGNORECASE): # Update to accomodate .paralog numbers from HybPiper
					print('%s was kept.' % sample_name)
					SeqIO.write(seqs, output_filename, 'fasta')
				else:
					print('%s was removed.' % sample_name)

		input_filename.close()
		output_filename.close()

print("DONE")

