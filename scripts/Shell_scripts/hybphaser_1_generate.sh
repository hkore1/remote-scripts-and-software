#!/bin/bash

######################
# Author: Benjamin Anderson
# Date: July 2022
# Modified: Oct 2022 (add intronerate option), 6 Nov 2022 (Harvey - modified for slurm)
# Description: run the first HybPhaser step on Gadi, using Lars' script in my singularity container
# Note: Enter input arguments below, requires: a file with a list of samples, one per line, and the path to the HybPiper output (decompressed)
#		also, set whether this is from Theo's wrapper (needs decompression) [default] or not
#		also, set whether to run with the intronerate supercontigs or not [default]
#		launch (where you want output folder).
#		The script will generate a folder called HybPhaser in the launch directory
######################

# Partition for the job:
#SBATCH --partition=physical

# Multithreaded (SMP) job: must run on one node 
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name="HybPhaser_generate"

# The project ID which this job should run under:
#SBATCH --account="punim1517"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=190G

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
#SBATCH --time=0-4:0:00

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi


##########################
### Set arguments here ###
##########################


# define the container being used
container="/data/gpfs/projects/punim1517/bin/hybphaser.sif"


# set the script for replacing NNNs when intronerate is set
nscript="/g/data/nm31/ben_anderson/scripts/fasta_Ns_sub.py"


# define the maximum number of samples to run at once (<= ncpus)
maxcontemp=24


# set up a log file
logfile="$(pwd)/hybphaser_generate_${SLURM_JOB_ID}.log"


# check args
samples="/data/gpfs/projects/punim1517/harvey_orel/projects/Rutaceae/Correa/HybPiper/done.txt"
hp2_dir="/data/gpfs/projects/punim1517/harvey_orel/projects/Rutaceae/Correa/HybPiper/results"
decompress="y"
intronerate="n"



###################################
### Leave everything below here ###
###################################

# Record and echo start time
start="$(date +%s)"
echo -e "\n**********\nStarting HybPhaser consensus generation at $(date)\n**********" >> $logfile 2>&1


# load the singularity module
module load singularity/3.5.3 >> $logfile 2>&1
#module load python3 #### HARVEY -- blanked out because doesnt work with singularity on SPARTAN and is only used for Ben's intronerate N's script

# Create a temporary directory (SCRATCH) and cd into it
export SCRATCH="/tmp/${SLURM_JOB_ID}"
mkdir -p $SCRATCH
cd $SCRATCH


# To make sure it doesn't go over size, split the samples into groups
split -l "$maxcontemp" "$samples" subset_
numbatches=$(find . -maxdepth 1 -name "subset_*" -printf '.' | wc -m)
echo -e "\nWill process $numbatches batches" >> $logfile 2>&1


# Move the HybPiper results to the working area
# If necessary, decompress the HybPiper results folders
index=1
if [ "$decompress" == "y" ]; then
	# Since there will be a lot to decompress, let's run in batches
	# Updated now to just run based on the split files (comment out old)
	#max_jobs="$PBS_NCPUS"
	for subset in subset_*
	do
		#cur_jobs=0
		echo -e "\n\nRunning decompression for batch $index" >> $logfile 2>&1
		for sample in $(cat $subset)
		do
			#((cur_jobs >= max_jobs)) && wait -n
			tar -xzf "$hp2_dir"/"$sample"/"$sample".tar.gz &
			#((++cur_jobs))
		done
		wait
		# Launch HybPhaser consensus generation
		echo -e "Running consensus generation for batch $index" >> $logfile 2>&1
		if [ "$intronerate" == "y" ]; then
			# first, we need to change the supercontigs to have 100 Ns instead of 10 when stitched
			echo -e "Substituting 100 Ns for 10 Ns in supercontig fasta files..." >> $logfile 2>&1
			python3 "$nscript" -m 100 */*/*/sequences/intron/*_supercontig.fasta
			# now, run the consensus generation
			echo -e "Generating consensus sequences...\n" >> $logfile 2>&1
			singularity exec -H "$(pwd)" "$container" 1_generate_consensus_sequences.sh -n "$subset" \
				-t "$maxcontemp" -o HybPhaser -i >> $logfile 2>&1
		else
			echo >> $logfile 2>&1
			singularity exec -H "$(pwd)" "$container" 1_generate_consensus_sequences.sh -n "$subset" \
				-t "$maxcontemp" -o HybPhaser >> $logfile 2>&1
		fi
		# manually remove read and mapping files and directories
		rm -r ./HybPhaser/*/*/reads/ ./HybPhaser/*/*/mapping_files/
		# move the HybPhaser results to the working directory
		rsync -rut HybPhaser $SLURM_SUBMIT_DIR/
		# remove the remaining files
		rm -r HybPhaser
		for sample in $(cat $subset)
		do
			rm -r "$sample"
		done
		((++index))
	done
else
	for subset in subset_*
	do
		for sample in $(cat $subset)
		do
			rsync -rt "$hp2_dir"/"$sample" .
		done
		# Launch HybPhaser consensus generation
		echo -e "\n\nRunning consensus generation for batch $index" >> $logfile 2>&1
		if [ "$intronerate" == "y" ]; then
			# first, we need to change the supercontigs to have 100 Ns instead of 10 when stitched
			echo -e "Substituting 100 Ns for 10 Ns in supercontig fasta files..." >> $logfile 2>&1
			python3 "$nscript" -m 100 */*/*/sequences/intron/*_supercontig.fasta
			# now, run the consensus generation
			echo -e "Generating consensus sequences...\n" >> $logfile 2>&1
			singularity exec -H "$(pwd)" "$container" 1_generate_consensus_sequences.sh -n "$subset" \
				-t "$maxcontemp" -o HybPhaser -i >> $logfile 2>&1
		else
			echo >> $logfile 2>&1
			singularity exec -H "$(pwd)" "$container" 1_generate_consensus_sequences.sh -n "$subset" \
				-t "$maxcontemp" -o HybPhaser >> $logfile 2>&1
		fi
		# manually remove read and mapping files and directories
		rm -r ./HybPhaser/*/*/reads/ ./HybPhaser/*/*/mapping_files/
		# move the HybPhaser results to the working directory
		rsync -rut HybPhaser $SLURM_SUBMIT_DIR/
		# remove the remaining files
		rm -r HybPhaser
		for sample in $(cat $subset)
		do
			rm -r "$sample"
		done
		((++index))
	done
fi


# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished HybPhaser consensus generation at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" >> $logfile 2>&1
