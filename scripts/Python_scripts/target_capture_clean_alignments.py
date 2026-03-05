#!/usr/bin/env python3

###############################
# Harvey Orel
# Date: Oct 2022
# Description: Cleans alignments output from Yang and Smith: 
#					* By sample (removes samples with some % of missing sites across alignment) &
#					* By site (removes sites with some % of gaps across each site)
# Note: Input is a directory of fasta alignments. Output is to "trimmed_alignments" folder in launch directory;
#		1_trimmed_sample_alignments = trimmed sample alignments only
#		2_trimmed_site_alignments = alignments trimmed by sample and site
#
#       Arguments are: <input_alignment_directory> <sample_gap_threshold> <site_gap_threshold>
#
# Threshold values retain samples/sites that are BELOW the threshold value, i.e.
# "python3 target_capture_filter_alignments.py g/data/nm31/.../18_alignments_stripped_names_MO/ 0.5 0.5"
# retains samples with less than 50% of sites missing (higher value = more lax), and sites with less than 50% gaps (higher values = more strict)
###############################


import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO

directory = str(os.getcwd())

input_directory = str(sys.argv[1])
sample_gap_threshold = float(sys.argv[2])
site_gap_threshold = float(sys.argv[3])

# Create output directory
if not os.path.exists(directory + "/trimmed_alignments/"):
	os.makedirs(directory + "/trimmed_alignments/")
	os.makedirs(directory + "/trimmed_alignments/1_trimmed_sample_alignments/")
	os.makedirs(directory + "/trimmed_alignments/2_trimmed_site_alignments/")
else:
	print("Error: directory 'trimmed_alignments' already exists")
	sys.exit(1)

# Check gap thresholds correctly specified
if (sample_gap_threshold > 1) or (sample_gap_threshold < 0):
	print('\n Sample gap threshold must be between 0 and 1!\n')
	sys.exit(1)

if (site_gap_threshold > 1) or (site_gap_threshold < 0):
	print('\n Site gap threshold must be between 0 and 1!\n')
	sys.exit(1)


# Iterate through files and remove samples below the sample gap threshold, write to new fasta files
for filename in os.listdir(input_directory):
	if filename.endswith(".fasta") or filename.endswith(".fas"):
		new_filename = filename.partition(".fasta")[0]
		print("Trimming samples from: " + filename)
		input_filename = open(input_directory + filename, 'r')
		output_filename = open(directory + "/trimmed_alignments/1_trimmed_sample_alignments/" + new_filename + "_samples_trimmed.fasta", 'w')

		for seqs in SeqIO.parse(input_filename, 'fasta'):
			sample_name = seqs.id
			sample_seq = seqs.seq
			alignment_length = len(seqs)
			gap_count = sample_seq.count("-")
			if (gap_count/float(alignment_length)) > sample_gap_threshold:
				print(' %s was removed.' % sample_name)
			else:
				SeqIO.write(seqs, output_filename, 'fasta')

		input_filename.close()
		output_filename.close()


# Iterate through stripped alignments and strip sites below the site gap threshold
directory = str(os.getcwd() + "/trimmed_alignments/")

# Iterate through sample-stripped alignments and remove sites below the site gap threshold
for filename in os.listdir(directory + "/1_trimmed_sample_alignments/"):
	if filename.endswith(".fasta") or filename.endswith(".fas"):
		new_filename = filename.partition("_samples_")[0]
		input_filename = directory + "/1_trimmed_sample_alignments/" + filename
		output_filename = open(directory + "/2_trimmed_site_alignments/" + new_filename + "_gaps_trimmed.fasta", 'w')

		try:
			alignment = AlignIO.read(open(input_filename, 'r'), "fasta")
		except ValueError:
			break
		n_samples = len(alignment)
		n = 0 # First column/site in alignment
		
		# Create index list of sites to keep in all samples
		sites2keep = []
		while n != alignment.get_alignment_length():
			site = alignment[:, n]
			n_gaps = site.count("-")
			if (n_gaps/n_samples) > 1 - site_gap_threshold: # Test for site gap threshold
				pass
			else:
				sites2keep.append(n)
			n += 1

		# Read in alignments with SeqIO (as sequences) and remove sites by index number in sites2keep
		alignment = SeqIO.parse(open(input_filename, 'r'), 'fasta')
		for sample in alignment:
			sample_name = sample.id
			new_sample_seq = list(sample.seq)
			new_sample_seq = [new_sample_seq[i] for i in sites2keep]

			sample.seq = Seq(''.join(new_sample_seq))
			SeqIO.write(sample, output_filename, 'fasta')
		
		output_filename.close()

