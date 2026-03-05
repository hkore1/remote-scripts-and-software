#!/usr/bin/env python3

###############################
# Harvey Orel
# Date: 28 Nov 2022
# Description: Runs a head-to-tail alignment comparison as per Simmons et al. (2022) Celastrales ms.
# Notes: Requires mafft and trimAl. Outputs to working directory.
#			-- On Spartan first 'module load mafft/7.453-with-extensions trimal/1.4.1 python/3.8.2'
#       Arguments are: <input_alignment_directory/>
###############################


import sys
import os
import shutil
import subprocess
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from io import StringIO


# Logic for arguments
if not len(sys.argv) > 1:
	sys.exit("Error: no linput alignments provided... Exiting")
else:
	pass

directory = str(os.getcwd())
input_directory = str(sys.argv[1])

# Create temp output directory (deleted in last few lines)
if not os.path.exists(directory + "/h2t_alignments/"):
	os.makedirs(directory + "/h2t_alignments/")
else:
	print("\nWarning: directory 'h2t_alignments' already exists... Overwriting")
	shutil.rmtree(directory + "/h2t_alignments")
	os.makedirs(directory + "/h2t_alignments/")


# Move seq files to h2t_alignments and create reverse orientation sequences
print("\nGenerating sequence files... \n")
for filename in os.listdir(input_directory):
	if filename.endswith(".fasta") or filename.endswith(".fas"):
		new_filename = filename.partition(".fasta")[0]
		input_filename = open(input_directory + filename, 'r')
		forward_output_filename = open(directory + "/h2t_alignments/" + new_filename + "_seqs.fasta", 'w')
		reverse_output_filename = open(directory + "/h2t_alignments/" + new_filename + "_seqs_reversed.fasta", 'w')

		for seqs in SeqIO.parse(input_filename, 'fasta'):
			sequence = seqs.seq
			SeqIO.write(seqs, forward_output_filename, 'fasta')
			reverse_sequence = sequence[::-1] # reverses the seq string 
			seqs.seq = reverse_sequence
			SeqIO.write(seqs, reverse_output_filename, 'fasta')

		input_filename.close()
		forward_output_filename.close()
		reverse_output_filename.close()


# Align each gene file with MAFFT - based on: https://github.com/mousepixels/sanbomics_scripts/blob/main/python_sequence_alignment.ipynb
print("Performing alignments... \n")
os.makedirs(directory + "/h2t_alignments/mafft/")

for filename in os.listdir("h2t_alignments/"):
	if filename.endswith(".fasta") or filename.endswith(".fas"):
		print("Aligning: " + filename)
		seqs = list(SeqIO.parse("h2t_alignments/" + filename, 'fasta'))

		seq_str = ''
		for seq in seqs:
			seq_str += '>' + seq.id + '\n'
			seq_str += str(seq.seq) + '\n'

		child = subprocess.Popen(['mafft', '--quiet', '-'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
		child.stdin.write(seq_str.encode())
		child_out = child.communicate()[0].decode('utf8')
		alignment = list(SeqIO.parse(StringIO(child_out), 'fasta'))
		child.stdin.close()

		SeqIO.write(alignment, "h2t_alignments/mafft/" + filename, "fasta")

	# Unreverse reversed alignments
	for filename in os.listdir("h2t_alignments/mafft"):
		if "reversed.fasta" in filename:
			print("Unreversing: " + filename)
			new_filename = filename.partition("_seqs_reversed.fasta")[0]
			input_filename = open("h2t_alignments/mafft/" + filename, 'r')
			output_filename = open("h2t_alignments/mafft/" + new_filename + "_seqs_unrev.fasta", 'w')

			for seqs in SeqIO.parse(input_filename, 'fasta'):
				sequence = seqs.seq
				reverse_sequence = sequence[::-1] # reverses the seq string 
				seqs.seq = reverse_sequence
				SeqIO.write(seqs, output_filename, 'fasta')

			input_filename.close()
			forward_output_filename.close()
			os.remove("h2t_alignments/mafft/" + filename)


# Read 'mafft/' filenames, create a list of paired forward & reverse alignments, run trimAl -compareset -ct
print("Running trimAl... \n")
os.makedirs(directory + "/h2t_output/")
os.makedirs(directory + "/h2t_output/alignments/")
os.makedirs(directory + "/h2t_output/stats/")

filelist = []
for filename in sorted(os.listdir("h2t_alignments/mafft/")):
	filelist.append(filename)

if len(filelist) % 2 == 0: # checks that filelist is an even number (i.e. all alignments are paired)
	ali_array = np.array_split(filelist, len(filelist)/2)
else:
	sys.exit("Error: generated an uneven number of alignments... Exiting")

for pair in ali_array:
	print("Comparing pair: " + str(pair))
	pair1 = "h2t_alignments/mafft/" + str(pair[0])
	pair2 = "h2t_alignments/mafft/" + str(pair[1])

	output_alignment_name = "h2t_output/alignments/" + str(pair[0].partition(".paralogs")[0]) + ".paralogs_H2T_alignment.fasta"
	output_stats_name = "h2t_output/stats/" + str(pair[0].partition(".paralogs")[0]) + ".paralogs_H2T_stats.html"
	output_log_name = "h2t_output/stats/" + str(pair[0].partition(".paralogs")[0]) + ".paralogs_H2T_log.log"

	with open('trimal_alignment_paths.txt', 'w') as f: # writes the pair to a temp txt file for trimal input
		f.write("%s\n%s" % (pair1,pair2))
	f.close()

	with open(output_log_name, 'w') as f:
		p2=subprocess.Popen("trimal -compareset trimal_alignment_paths.txt -ct 0.5 -out %s -htmlout %s >> %s" % (output_alignment_name,output_stats_name,output_log_name), shell=True).wait()
	f.close()

shutil.rmtree(directory + "/h2t_alignments")
os.remove("trimal_alignment_paths.txt")

print("\nDONE!!!")


