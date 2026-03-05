#!/bin/bash

######################
# Author: Harvey Orel (modelled on Ben Anderson's hybphaser_evaluate.sh script on GADI)
# Date: 6 Nov 2022 
# Description: run HybPhaser evaluation of the results of a HybPiper run
# Note: Unlike Ben's script, it doesn't use his container (as I was having problems with it)
# Note: requires a file structure and consensus reads already present (from hybphaser_readmap.sh or hybphaser_generate.sh)
#		Below, enter args for the output directory (where the file structure is, e.g. has the folder 01_data present),
#		the targets file (used for the HP2 run), and the namelist file (one sample ID per line)
#		also, specify whether the output to evaluate is intronerate or not [default]
######################

# Partition for the job:
#SBATCH --partition=physical

# Multithreaded (SMP) job: must run on one node 
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name="HybPhaser_evaluate"

# The project ID which this job should run under:
#SBATCH --account="punim1517"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=5G

# Use this email address:
#SBATCH --mail-user=horel@student.unimelb.edu.au

# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL
# begins
#SBATCH --mail-type=BEGIN
# ends successfully
#SBATCH --mail-type=END

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-00:10:00

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi


# define the container being used
container="/data/gpfs/projects/punim1517/bin/hybphaser.sif"

# path to HybPhaser folder in bin
hybphaser_path="/data/gpfs/projects/punim1517/bin/HybPhaser-main"

# check args
out_dir="$(pwd)/HybPhaser"
target_file="/data/gpfs/projects/punim1517/target_files/lamiales_targetfile.fasta"
namelist="/data/gpfs/projects/punim1517/harvey_orel/projects/Myoporeae/hp2run/done.txt"
intronerate="y"


#############################
### Dont alter below this ###
#############################


# Set the original launch directory
launch_dir="$(pwd)"


# set up a master log file
logfile="$(pwd)/hybphaser_evaluate_${SLURM_JOB_ID}.log"


# load the modules
module load gcccore/8.3.0 >> $logfile 2>&1
module load r/4.0.0 >> $logfile 2>&1
module list >> $logfile 2>&1


# Record and echo start time and input
start="$(date +%s)"
echo -e "\n**********\nStarting HybPhaser evaluation at $(date)\n**********" >> $logfile 2>&1


# We need to create the config.txt file in the out_dir if it doesn't exist already
config_file="${out_dir}/config.txt"
if [ ! -f $config_file ]; then
	echo "# General settings" >> $config_file
	echo "path_to_output_folder = \"${out_dir}\"" >> $config_file
	echo "fasta_file_with_targets = \"${target_file}\"" >> $config_file
	echo "targets_file_format = \"DNA\"" >> $config_file
	echo "path_to_namelist = \"${namelist}\"" >> $config_file
	echo "intronerated_contig = \"$intronerate\"" >> $config_file
	echo "" >> $config_file

	echo "name_for_dataset_optimization_subset = \"\"" >> $config_file
	echo "" >> $config_file
	echo "# Missing data" >> $config_file
	echo "remove_samples_with_less_than_this_propotion_of_loci_recovered = 0" >> $config_file
	echo "remove_samples_with_less_than_this_propotion_of_target_sequence_length_recovered = 0" >> $config_file
	echo "remove_loci_with_less_than_this_propotion_of_samples_recovered = 0" >> $config_file
	echo "remove_loci_with_less_than_this_propotion_of_target_sequence_length_recovered = 0" >> $config_file
	echo "" >> $config_file

	echo "# Paralogs" >> $config_file
	echo "remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs = \"none\"" >> $config_file
	echo "file_with_putative_paralogs_to_remove_for_all_samples = \"\"" >> $config_file
	echo "remove_outlier_loci_for_each_sample = \"no\"" >> $config_file
	echo "" >> $config_file
else
	echo -e "\nDetected a config file, so will use it for the run\n"
fi


# Now, we change to the outdirectory and run the R script to count SNPs
cd $out_dir
Rscript ${hybphaser_path}/1a_count_snps.R >> $logfile 2>&1


# Finally, we run the R script to assess the dataset
Rscript ${hybphaser_path}/1b_assess_dataset.R >> $logfile 2>&1
cd $launch_dir


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished HybPhaser runs at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" >> $logfile 2>&1
