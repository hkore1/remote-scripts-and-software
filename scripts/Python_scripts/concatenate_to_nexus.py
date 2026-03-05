#!/usr/bin/env python3

###############################
# Harvey Orel
# Date: 24 Sep 2024
# Description: Concatenates faasta alignments into a nexus file and records partitions as charsets.
# Note: Input is a directory of fasta alignments. Output is to "concatenated_alignments.nex" in the launch directory;
#
#       Arguments are: <input_alignment_directory>
###############################

from Bio import SeqIO
from Bio.Nexus import Nexus
from Bio.Align import MultipleSeqAlignment

import os
import re
import sys
import argparse

parser = argparse.ArgumentParser(description='Concatenate alignments and convert to nexus with partitions')

parser.add_argument('-a', metavar='--alignment_folder', type=str, nargs=1, help='Path to directory with .fasta alignments')

args = parser.parse_args()


#### Confirm script parameters ####
align_folder = str(args.a[0])

#### Define functions ####
def convert_to_nexus(directory,file):
	base = os.path.splitext(file)[0]
	return SeqIO.convert(in_file = directory + file, in_format = "fasta", out_file = directory + "/" + base + ".nex", out_format = "nexus", molecule_type = "DNA")



# Run the main loops

for filename in os.listdir(align_folder):
	if filename.endswith(".fasta"):
		convert_to_nexus(directory=str(os.getcwd()) + "/" + align_folder, file=filename)

file_list = list()
for filename in os.listdir(align_folder):
	if filename.endswith(".nex"):
		file_list.append(str(align_folder + filename))


concat = [(fname, Nexus.Nexus(fname)) for fname in file_list]
combined = Nexus.combine(concat)
combined.write_nexus_data(filename=open("concatenated_alignments.nex", "w"))

for filename in os.listdir(align_folder):
	if filename.endswith(".nex"):
		os.remove(os.path.join(align_folder, filename))