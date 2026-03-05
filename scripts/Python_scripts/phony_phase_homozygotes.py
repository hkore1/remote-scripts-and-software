#!/usr/bin/env python3

###############################
# Harvey Orel
# Date: 26 Apr 2023
#
# Description: Duplicates homozygote samples (from my Allele Phasing workflow) so that they have 2 identical copies, labelled with
#               '_h1' and '_h2' for better amalgamation of phylogenetic results
#
# Note: Input is a directory of fasta alignments. Output is to "PhonyPhasedAlignments/" in launch directory;
#
#       Arguments are: <input_alignment_directory>
###############################

import os
import sys
import re

# Get cwd
directory = str(os.getcwd())

# input directory containing alignment files
input_directory = str(sys.argv[1])

# Create output directory
if not os.path.exists(directory + "/PhonyPhasedAlignments/"):
    os.makedirs(directory + "/PhonyPhasedAlignments/")
else:
    print("Error: directory 'PhonyPhasedAlignments' already exists... Cancelling")
    sys.exit(1)

for filename in os.listdir(input_directory):
    if filename.endswith(".fasta") or filename.endswith(".fas") or filename.endswith(".FNA"):
        basename = os.path.basename(filename)
        with open(input_directory + filename, 'r') as infile, open(directory + "/PhonyPhasedAlignments/" + basename + '_no_homozygotes.fasta', 'w') as outfile:
            sequences = infile.read().split('>')[1:]
            for seq in sequences:
                seq_parts = seq.split('\n')
                label = seq_parts[0]
                sequence = ''.join(seq_parts[1:])
                if '_h1' not in label and '_h2' not in label:
                    print("Homozygote %s detected in file %s" % (label, basename))
                    # duplicate sequence
                    label_h1 = re.sub(r'\s', '_h1 ', label, count=1)
                    sequence_h1 = sequence
                    label_h2 = re.sub(r'\s', '_h2 ', label, count=1)
                    sequence_h2 = sequence
                    # update labels
                    seq_parts[0] = label_h1
                    seq_parts_h2 = seq_parts.copy()
                    seq_parts_h2[0] = label_h2
                    # write new sequences to output file
                    outfile.write('>' + '\n'.join(seq_parts) + '\n')
                    outfile.write('>' + '\n'.join(seq_parts_h2) + '\n')
                else:
                    # write original sequence to output file
                    outfile.write('>' + seq + '\n')