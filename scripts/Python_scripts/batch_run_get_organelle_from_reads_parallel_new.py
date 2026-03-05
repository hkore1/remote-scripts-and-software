#!/usr/bin/env python3

###########################################
# Harvey Orel
# Date: 15 Oct 2022
# 	Modified: 25 May 2023 - Major update to include running samples in parallel.
# Description: Batch runs GetOrganelle (Jin et al. 2020) on a directory of paired-end read files. 
# Usage: Run from the command line, in the following format:
#
#				python3 batch_run_get_organelle_from_reads.py -go <###> -i <###> -a "<###>" -o <###> -t <###>
#
#		Arguments: 	-go = Path to directory with get_organelle_from_reads.py
#					-i = Path to directory with input read files
#					-a = A string of settings for GetOrganelle assembly, as they would be entered for a single GetOrganelle run
#					-o = The name for the output parent directory (containing subfolders of assemblies)
#					-t = The number of samples to run at a time. The product of this and the -t value applied to GetOrganelle settings should equal
#						 the total number of available CPUs
#					-c = compress .fq and .sam files? "yes" or "no"
#					-l = use a list of readfile names to run only on certain samples in the input directory
#
# Note: If no assembly settings are provided, it will run with some default settings I've found work for me.
#		It will automatically exclude any non fastq or fastq.gz files, or any files that don't have corresponding
# 		R1 and R2 files, that are in the input directory. Outputs to the '-o' folder in the
#		working directory. The settings --spades-options "--careful --memory 1024" and --zip-files are employed by default, so do not 
#		specify these in -a
#
# Example commands to run (note quotations around assmebly settings): 
#	python3 batch_run_get_organelle_from_reads.py -go /Users/JohnSmith/opt/miniconda3/bin/ -i read_files/
#	python3 batch_run_get_organelle_from_reads.py -go /Users/JohnSmith/opt/miniconda3/bin/ -i read_files/ -a "-R 5 -k 21,45 -F embplant_pt -t 4 --out-per-round --memory-save" -o getorg_aseemblies -t 2
###########################################

import sys
import glob
import time
import os
import argparse
import re
import logging
import psutil
import concurrent.futures
import subprocess

parser = argparse.ArgumentParser(description='Batch run GetOrganelle')

parser.add_argument('-go', metavar='--GetOrganelle_directory', type=str, nargs=1, help='path to directory with get_organelle_from_reads.py')
parser.add_argument('-i', metavar='--input_directory', type=str, nargs=1, help='path to directory with input read files')
parser.add_argument('-a', metavar='--assembly_settings', type=str, nargs='+', help='settings for GetOrganelle assembly, as a string')
parser.add_argument('-o', metavar='--output_directory', type=str, nargs=1, help='name for output directory, as a string')
parser.add_argument('-t', metavar='--threads', type=int, nargs=1, help='number of threads (i.e. samples to run at a time)')
parser.add_argument('-c', metavar='--compress', type=str, nargs=1, help='compress .fq and .sam files? "yes or "no"')
parser.add_argument('-l',metavar='--list', type=str, nargs=1, help='optional .txt list to subset runs to certain samples in the input_directory, one line per read file to include')

args = parser.parse_args()

print("\nListing input arguments...\n")




###############################
### Logic for  arguments... ###
###############################

# Test path to getorganelle directory
if args.go is None:
	print("No direction to get_organelle_from_reads.py, EXITING...!!!")
	sys.exit()
else:
	print("Path to get_organelle_from_reads.py: '" + str(args.go[0]) + "'")
	pass

# Test input file directory
if args.i is None:
	print("No path to input files, EXITING...!!!")
	sys.exit()
else:
	print("Path to input files: '" + str(args.i[0]) + "'")
	pass

# Check assembly settings
if args.a is None:
	print("No assembly settings provided, using defaults: '-R 15 -k 21,45,65,85,105 -F embplant_pt -t 1 --out-per-round --memory-save -w 0.67'\n")
	args.a = ["-R","15","-k","21,45,65,85,105","-F","embplant_pt","-t","1","--out-per-round","--memory-save","-w","0.67"]
elif args.a == 'arg_was_not_given':
	print("Option given, but no command-line argument '-a', using defaults: '-R 15 -k 21,45,65,85,105 -F embplant_pt -t 1 --out-per-round --memory-save -w 0.67'\n")
	args.a = ["-R","15","-k","21,45,65,85,105","-F","embplant_pt","-t","1","--out-per-round","--memory-save","-w","0.67"]
else:
	print("Using assembly settings provided: '" + str(' '.join(args.a)) + "'\n")

# Test output directory
if args.o is None:
	print("No output directory name provided, using 'getorganelle_assemblies'")
	args.o = "getorganelle_assemblies"
else:
	print("Name of output directory: '" + str(args.o[0]) + "'")
	pass

# Test threads
if args.t is None:
	print("No thread number specified, using 1")
	args.t = 1
else:
	print("Number of threads: '" + str(args.t[0]) + "'")
	pass

# Test for list
if args.l is None:
	pass
else:
	print("Running on a subset of samples in input directory based on the list: '" + str(args.l[0]) + "'")
	subset_list = str(args.l[0])


### Confirm parameters for script ###
getorganelle_dir = str(args.go[0])
input_dir = str(args.i[0])
assembly_settings = ' '.join(args.a)
output_dir = str(os.getcwd() + '/' + args.o[0])
samplethreads = int(args.t[0])
compression = str(args.c[0])




########################################
### Load get_organelle_from_reads.py ###
########################################

sys.path.insert(1, getorganelle_dir)

try:
	import get_organelle_from_reads
	print("Successfully loaded 'get_organelle_from_reads.py'")
except ModuleNotFoundError:
	print("Error importing 'get_organelle_from_reads.py', check file path...")
	sys.exit()




##########################################################
###  Organise list of read files to run the program on ###
##########################################################

### Coerce read files in input directory to list of tuples ###
#create list of fastq files
filenames = []
for filename in os.listdir(input_dir):
	if filename.endswith(".fastq") or filename.endswith(".fastq.gz"):
		filenames.append(filename)

#test if list argument was provided and, if so, remove items from filenames not matching those in the list
if 'subset_list' in globals():
	with open(subset_list) as txtfile:
		to_run = [line.rstrip() for line in txtfile]
	filenames = set(filenames).intersection(to_run)
else:
	pass

#get base file names (filename preceding R1 or R2) and put matches in tuple
tuples = []
for file in filenames:
	pattern = re.findall(r"(.+)(_R\d{1})", file)[0][0]
	matches = tuple([x for x in filenames if x.startswith(pattern)])
	tuples.append(matches)

#remove any items from list that arent a tuple with 2 elements
removed_files = []
tuples_cleaned = []
for i in tuples:
	if len(i) != 2:
		removed_files.append(i)
	else:
		tuples_cleaned.append(i)

#print some feedback (and sort tuples so that R1 is always first)
print("\nRemoved " +  str(len(removed_files)) + " files:")
for i in removed_files:
	print(i)

#sorting
tuples_cleaned_sorted = []
for i in tuples_cleaned:
	i = sorted(i)
	tuples_cleaned_sorted.append(tuple(i))

#finalise tuple list
tuples_cleaned_sorted = list(set(tuples_cleaned_sorted))
print("\nLocated " + str(len(tuples_cleaned_sorted)) + " sets of paired read files for assembly:")
print(tuples_cleaned_sorted)




###############################
### Create output directory ###
###############################

if not os.path.exists(output_dir):
	os.makedirs(output_dir)
elif os.path.exists(output_dir):
	print("Warning: prexisting " + str(args.o[0]) + " folder in the working directory, new files will write to that directory")




########################
### Define functions ###
########################

def compress_with_tar(source_path,source):
	print(source_path + "/" + source)
	if os.path.isfile(source_path + "/" + source):
		print("Filename: " + source)
		run_command = 'tar -czvf ' + source_path + '/' + source + '.tar.gz ' + source_path + '/' + source + ' --remove-files'
		p3 = subprocess.run(run_command, shell=True)
	else:
		pass

def run_get_org(sample,compress):
	# First, get the information for an individual sample...
	pattern = str(re.findall(r".*(?=_R\d{1})", sample[0])[0])
	readfile1 = input_dir + str(sample[0])
	readfile2 = input_dir + str(sample[1])
	output = output_dir + "/" + pattern + "_output"
	sample_name=re.split(r"_H.{8}_", pattern)[0]

	# If output for a sample exists already (e.g. from a previous run) then exit the function and move on to next sample
	if os.path.isdir(output):
		print("The directory: " + output + " already exists, moving on to next sample")
		return

	# Define timeout length of 3 hrs for moving on from failed samples
	timeout = time.time() + 60*180 # 180 mins in future

	# Run get_organelle_from_reads.py
	p1 = subprocess.run('python3 %s/get_organelle_from_reads.py -1 %s -2 %s -o %s %s --spades-options "--careful --memory 1024" --prefix %s_' %(getorganelle_dir,readfile1,readfile2,output,assembly_settings,sample_name), shell=True)

	# Loop to stop function exiting before an assembly graph is generated, or timeout if assembly takes longer than 3 hrs
	while True:
		if glob.glob(str(output + '/*assembly_graph.fastg')):
			break
		elif time.time() > timeout:
			return

	# If assembly is successful...
	print("Running mapping evaluation on: " + sample_name)
	if time.time() < timeout:
		
		# Get output info from previous command, for input in next command
		getorg_fasta = str([s for s in os.listdir(output) if "path_sequence.fasta" in s][0])

		# Run evaluate_assembly_using_mapping.py
		p2=subprocess.run("python3 %s/evaluate_assembly_using_mapping.py -f %s/%s -1 %s/%s_extended_1_paired.fq -2 %s/%s_extended_2_paired.fq -o %s/%s_assembly_evaluation -c yes --draw" %(getorganelle_dir,output,getorg_fasta,output,sample_name,output,sample_name,output,sample_name), shell=True).wait()

		if compress == "yes":
			print("Compressing .fq and anti-seed .sam files (if they are present)")
			compress_with_tar(output, sample_name + "_extended_1_paired.fq")
			compress_with_tar(output, sample_name + "_extended_1_unpaired.fq")
			compress_with_tar(output, sample_name + "_extended_2_paired.fq")
			compress_with_tar(output, sample_name + "_extended_2_unpaired.fq")
			compress_with_tar(output, sample_name + "_anti_seed_bowtie.sam")
		else:
			pass




######################
### Run the script ###
######################

if __name__ == '__main__':

	print("\n")
	print("******************************")
	print("*** Running GetOrganelle *****")
	print("******************************")
	print("\n")

	logging.basicConfig(level=logging.DEBUG, filename="logfile.log", filemode="a+",
						format="%(asctime)-15s %(levelname)-8s %(message)s")
	logging.info("Begin execution of 'batch_run_get_organelle_from_reads.py'")

	executor = concurrent.futures.ProcessPoolExecutor(samplethreads)
	futures = [executor.submit(run_get_org,sample,compression) for sample in tuples_cleaned_sorted]
	concurrent.futures.as_completed(futures)


