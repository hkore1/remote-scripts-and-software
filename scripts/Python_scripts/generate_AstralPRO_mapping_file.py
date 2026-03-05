#!/usr/bin/env python

import argparse

# create an ArgumentParser object
parser = argparse.ArgumentParser(description='Process input file')

# add argument for input file
parser.add_argument('input_file', type=str, help='path to input file')

# parse the command-line arguments
args = parser.parse_args()

outfile = args.input_file.split(".")[0] + "_apro_mapping.txt"

with open(args.input_file, "r") as input_file, open(outfile, "w") as output_file:
    for line in input_file:
        if line.startswith(">"):
            sample_name = line.strip().split()[0]
            if "_h1" in sample_name:
                sample_name_full = sample_name
                sample_name_full = sample_name_full.replace(">", "")
                sample_name = sample_name.replace("_h1", "")
                sample_name = sample_name.replace(">", "")
            elif "_h2" in sample_name:
                sample_name_full = sample_name
                sample_name_full = sample_name_full.replace(">", "")
                sample_name = sample_name.replace("_h2", "")
                sample_name = sample_name.replace(">", "")
            else:
                sample_name_full = sample_name
                sample_name_full = sample_name_full.replace(">", "")
                sample_name = sample_name.replace(">", "")
            #gene_name = args.input_file.split(".")[0]
            #sample_name_full = sample_name_full.replace(f"{sample_name}", "")
            output_file.write(f"{sample_name_full} {sample_name}\n")
