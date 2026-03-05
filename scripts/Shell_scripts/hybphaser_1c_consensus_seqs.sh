#!/bin/bash

######################
# Author: Harvey Orel (modelled on Ben Anderson's hybphaser_evaluate.sh script on GADI)
# Date: 7 Nov 2022 
# Description: run HybPhaser initial sequence list generation
# Note: This must be done from shell, as it requires also calling scripts 1a and 1b to generate the "failed_loci" object (see https://github.com/LarsNauheimer/HybPhaser/issues/6)
#		Initial results of hybphaser_1a1b_evaluate.sh are deleted and run again, this time with parameters for sequence generation employed.
######################

# Partition for the job:
#SBATCH --partition=physical

# Multithreaded (SMP) job: must run on one node 
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name="HybPhaser_1c_consensus"

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
#SBATCH --time=0-00:20:00

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
intronerate="n"

# Set the parameters for HybPhaser trimming
name_for_dataset_optimization_subset="" #string

remove_samples_with_less_than_this_propotion_of_loci_recovered=0 #numeric
remove_samples_with_less_than_this_propotion_of_target_sequence_length_recovered=0 #numeric
remove_loci_with_less_than_this_propotion_of_samples_recovered=0 #numeric
remove_loci_with_less_than_this_propotion_of_target_sequence_length_recovered=0 #numeric
remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs=0.02 #numeric

file_with_putative_paralogs_to_remove_for_all_samples="" #string
remove_outlier_loci_for_each_sample="no" #string



#############################
### Dont alter below this ###
#############################


# Set the original launch directory
launch_dir="$(pwd)"


# set up a master log file
logfile="$(pwd)/hybphaser_1c_consensus_${SLURM_JOB_ID}.log"


# load the modules
module load gcccore/8.3.0 >> $logfile 2>&1
module load r/4.0.0 >> $logfile 2>&1
module list >> $logfile 2>&1


# Record and echo start time and input
start="$(date +%s)"
echo -e "\n**********\nStarting HybPhaser evaluation at $(date)\n**********" >> $logfile 2>&1

rm -r "${out_dir}/00_R_objects" >> $logfile 2>&1
rm -r "${out_dir}/02_assessment" >> $logfile 2>&1


# We need to create the config.txt file in the out_dir if it doesn't exist already
config_file="${out_dir}/config.txt"
if [ ! -f $config_file ]; then
	echo -e "######################################################"  >> $config_file
	echo -e "### Configuration File for all HybPhaser R scripts ###"  >> $config_file
	echo -e "######################################################\n"  >> $config_file

	echo "# General settings" >> $config_file
	echo "path_to_output_folder = \"${out_dir}\"" >> $config_file
	echo "fasta_file_with_targets = \"${target_file}\"" >> $config_file
	echo "targets_file_format = \"DNA\"" >> $config_file
	echo "path_to_namelist = \"${namelist}\"" >> $config_file
	echo "intronerated_contig = \"$intronerate\"" >> $config_file
	echo "" >> $config_file

	echo -e "###############################"  >> $config_file
	echo -e "### Part 1:  SNP Assessment ###"  >> $config_file
	echo -e "###############################\n"  >> $config_file

	echo "name_for_dataset_optimization_subset = \"${name_for_dataset_optimization_subset}\"" >> $config_file
	echo "" >> $config_file
	echo "# Missing data" >> $config_file
	echo "remove_samples_with_less_than_this_propotion_of_loci_recovered = "${remove_samples_with_less_than_this_propotion_of_loci_recovered}"" >> $config_file
	echo "remove_samples_with_less_than_this_propotion_of_target_sequence_length_recovered = "${remove_samples_with_less_than_this_propotion_of_target_sequence_length_recovered}"" >> $config_file
	echo "remove_loci_with_less_than_this_propotion_of_samples_recovered = "${remove_loci_with_less_than_this_propotion_of_samples_recovered}"" >> $config_file
	echo "remove_loci_with_less_than_this_propotion_of_target_sequence_length_recovered = "${remove_loci_with_less_than_this_propotion_of_target_sequence_length_recovered}"" >> $config_file
	echo "" >> $config_file

	echo "# Paralogs" >> $config_file
	echo "remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs = "${remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs}"" >> $config_file
	echo "file_with_putative_paralogs_to_remove_for_all_samples = \"${file_with_putative_paralogs_to_remove_for_all_samples}\"" >> $config_file
	echo "remove_outlier_loci_for_each_sample = \"${remove_outlier_loci_for_each_sample}\"" >> $config_file
	echo "" >> $config_file
else
	rm -r $config_file
	echo -e "######################################################"  >> $config_file
	echo -e "### Configuration File for all HybPhaser R scripts ###"  >> $config_file
	echo -e "######################################################\n"  >> $config_file

	echo "# General settings" >> $config_file
	echo "path_to_output_folder = \"${out_dir}\"" >> $config_file
	echo "fasta_file_with_targets = \"${target_file}\"" >> $config_file
	echo "targets_file_format = \"DNA\"" >> $config_file
	echo "path_to_namelist = \"${namelist}\"" >> $config_file
	echo "intronerated_contig = \"$intronerate\"" >> $config_file
	echo "" >> $config_file
	
	echo -e "###############################"  >> $config_file
	echo -e "### Part 1:  SNP Assessment ###"  >> $config_file
	echo -e "###############################\n"  >> $config_file

	echo "name_for_dataset_optimization_subset = \"${name_for_dataset_optimization_subset}\"" >> $config_file
	echo "" >> $config_file
	echo "# Missing data" >> $config_file
	echo "remove_samples_with_less_than_this_propotion_of_loci_recovered = "${remove_samples_with_less_than_this_propotion_of_loci_recovered}"" >> $config_file
	echo "remove_samples_with_less_than_this_propotion_of_target_sequence_length_recovered = "${remove_samples_with_less_than_this_propotion_of_target_sequence_length_recovered}"" >> $config_file
	echo "remove_loci_with_less_than_this_propotion_of_samples_recovered = "${remove_loci_with_less_than_this_propotion_of_samples_recovered}"" >> $config_file
	echo "remove_loci_with_less_than_this_propotion_of_target_sequence_length_recovered = "${remove_loci_with_less_than_this_propotion_of_target_sequence_length_recovered}"" >> $config_file
	echo "" >> $config_file

	echo "# Paralogs" >> $config_file
	echo "remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs = "${remove_loci_for_all_samples_with_more_than_this_mean_proportion_of_SNPs}"" >> $config_file
	echo "file_with_putative_paralogs_to_remove_for_all_samples = \"${file_with_putative_paralogs_to_remove_for_all_samples}\"" >> $config_file
	echo "remove_outlier_loci_for_each_sample = \"${remove_outlier_loci_for_each_sample}\"" >> $config_file
fi


# Now, we change to the outdirectory and run the R script to count SNPs
cd $out_dir
Rscript ${hybphaser_path}/1a_count_snps.R >> $logfile 2>&1
wait

# concatenate the R scripts to assess the dataset (1b) and to trim the dataset (1c)
cat ${hybphaser_path}/1b_assess_dataset.R ${hybphaser_path}/1c_generate_sequence_lists.R > concat_1b1c.R

# Finally, run the concatenated R script
Rscript concat_1b1c.R >> $logfile 2>&1
wait
rm concat_1b1c.R

cd $launch_dir


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished HybPhaser runs at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" >> $logfile 2>&1
