#!/usr/bin/env python3

################################
# Harvey Orel
# Date: 1 June 2023
# Description: Runs IQ-TREE2 on a folder of fasta alignments in parallel.
# Usage: Run from the command line, in the following format:
#
#				python3 run_iqtree_in_parallel.py -af <###> -s "<###>" -c <###> -t <###>
#
#		Arguments: 	-af = The path to directory containing .fasta alignments for analysis.
#                   -s = The settings for IQ-TREE.
#					-c = The number of CPUs to devote to each alignment.
#					-t = The number of samples to run at a time. The product of this and -c should equal
#						 the total number of available CPUs.
#
# Note: Assumes iqtree2 is installed in the users system path.
#
# Example commands to run (note quotations around assembly settings): 
#	python3 ...
################################

import os
import argparse
import subprocess
import concurrent.futures


parser = argparse.ArgumentParser(description='Batch run IQ-TREE2 in parallel')

parser.add_argument('-af', metavar='--alignment_folder', type=str, nargs=1, help='Path to directory with .fasta alignments')
parser.add_argument('-s', metavar='--iqtree_settings', type=str, nargs='+', help='Additional settings for IQ-TREE, as a string')
parser.add_argument('-c', metavar='--num_cpus', type=int, nargs=1, help='Number of CPUs to devote to each alignment')
parser.add_argument('-t', metavar='--num_threads', type=int, nargs=1, help='Number of threads (i.e. alignments to run at a time)')

args = parser.parse_args()


#### Confirm script parameters ####
align_folder = str(args.af[0])
iqtree_settings = ' '.join(args.s)
cpus_per_run = int(args.c[0]) # Number of CPUs to use per IQ-TREE run
num_parallel = int(args.t[0]) # Number of parallel IQ-TREE runs


#### Define Functions ####

def run_iqtree(align_file, cpus, additional_commands):
    # Define the IQ-TREE command with desired options
    iqtree_cmd = f"iqtree2 -s {align_file} -nt {cpus} {additional_commands}"
    
    # Run IQ-TREE command using subprocess
    subprocess.run(iqtree_cmd, shell=True)
    
def run_iqtree_parallel(align_folder, iq_tree_commands, cpus_per_run, num_parallel):
    # Get a list of alignment files in the folder
    align_files = [f for f in os.listdir(align_folder) if f.endswith(".fasta")]
    
    # Create a ThreadPoolExecutor for parallel execution
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_parallel) as executor:
        # Submit IQ-TREE runs in parallel
        futures = []
        for align_file in align_files:
            futures.append(executor.submit(run_iqtree, align_folder + "/" + align_file, cpus_per_run, iq_tree_commands))
        
        # Wait for all IQ-TREE runs to complete
        concurrent.futures.wait(futures)


######################
### Run the script ###
######################

if __name__ == '__main__':

	print("\n")
	print("**************************************")
	print("**** Running IQ-TREE2 in parallel ****")
	print("**************************************")
	print("\n")

	run_iqtree_parallel(align_folder, iqtree_settings, cpus_per_run, num_parallel)
