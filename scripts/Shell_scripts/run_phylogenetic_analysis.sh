#!/bin/bash

################################################################
# Run single phylogenetic analysis on a folder of alignments (IQTREE concat, gene trees, concordance factors; ASTRAL on gene trees)
# Author: Harvey Orel
# Usage: Set script parameters manually - things above the "Leave below here" comment. Run from working 
#			directory.
################################################################

# Partition for the job:
#SBATCH --partition=physical

# Multithreaded (SMP) job: must run on one node 
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name="phylo_run"

# The project ID which this job should run under:
#SBATCH --account="punim1517"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=32G

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
#SBATCH --time=1-00:0:00

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi


# Run the job from the directory where it was launched (default)

# The job command(s):

# cd to launch directory
cd $SLURM_SUBMIT_DIR

# set the input alignment directories
export INPUT="/data/gpfs/projects/punim1517/harvey_orel/projects/Rutaceae/Eriostemon_Group_analyses/7b_phasing_phylogenetics/2_phonyphased_1to1_ortho_filtered_alignments"

# set the singularity container
container="/data/projects/punim1517/bin/hybpiper-paragone_latest.sif"



########################
### Leave below here ###
########################

# set up a master log file
logfile="$SLURM_SUBMIT_DIR/phylo_$SLURM_JOB_ID.log"

# load modules
module load singularity/3.5.3 >> $logfile 2>&1
module load java/11.0.2 

# define alignment base directories
export BASE=$(basename $INPUT)

# Record and echo start time and input
start="$(date +%s)"
echo -e "****************************************************************" >> $logfile 2>&1
echo -e "*** Starting phylogenetics run at $(date) ***" >> $logfile 2>&1
echo -e "****************************************************************" >> $logfile 2>&1
echo -e "Provided path to alignments: $INPUT" >> $logfile 2>&1

# Copy script to submit directory for future reference of employed settings
\cp /data/gpfs/projects/punim1517/scripts/run_phylogenetic_analysis.sh $SLURM_SUBMIT_DIR; mv run_phylogenetic_analysis.sh logged_script_used.sh

#################### Run analyses on MO ####################
if [ ! -d "$SLURM_SUBMIT_DIR/phylogenetic_analysis" ]; then
	echo -e "****************************************************************" >> $logfile 2>&1
	echo -e "*** Beginning Phylogenetic Analyses at $(date) ***" >> $logfile 2>&1
	echo -e "****************************************************************" >> $logfile 2>&1

	rm -r phylogenetic_analysis; mkdir phylogenetic_analysis; cd phylogenetic_analysis/

	## Run a concatenation tree using IQTREE
	rm -r concatenated_IQTREE; mkdir concatenated_IQTREE; cd concatenated_IQTREE/
	cp -R $INPUT $PWD
	singularity exec "$container" iqtree --prefix concat -p $BASE --threads-max "48" -T AUTO -B 1000
	rm -r $BASE; cd ..
	echo -e "Concatenated tree complete...   $(date)" >> $logfile 2>&1

	## Run gene trees for each alignment using IQTREE
	rm -r genetrees_IQTREE; mkdir genetrees_IQTREE; cd genetrees_IQTREE/
	cp -R $INPUT $PWD
	cd $BASE
	for file in *; do singularity exec "$container" iqtree -s "$file" --threads-max "48" -T AUTO -B 1000; done
	cat *.treefile > loci.treefile; mv loci.treefile ../
	cat *.log > loci.log; mv loci.log ../
	cat *.bionj > loci.bionj; mv loci.bionj ../
	cat *.iqtree > loci.iqtree; mv loci.iqtree ../
	cat *.mldist > loci.mldist; mv loci.mldist ../
	cd ..; rm -r $BASE; cd ..
	echo -e "Individual gene trees complete...   $(date)" >> $logfile 2>&1

	## Calculate concordance factors on the concatenation tree
	rm -r concordance_IQTREE; mkdir concordance_IQTREE; cd concordance_IQTREE/
	cp -R $INPUT $PWD
	cp ../concatenated_IQTREE/concat.treefile $PWD
	cp ../genetrees_IQTREE/loci.treefile $PWD
	singularity exec "$container" iqtree -t concat.treefile --gcf loci.treefile -p $BASE --scf 100 --prefix concord
	rm -r $BASE; cd ..
	echo -e "Concordance factors complete...   $(date)" >> $logfile 2>&1

	## Run an ASTRAL tree using the gene trees
	rm -r species_tree_ASTRAL; mkdir species_tree_ASTRAL
	java -Xmx8000M -jar /data/projects/punim1517/bin/Astral/astral.5.7.8.jar -i "$(pwd)/genetrees_IQTREE/loci.treefile" -o "$(pwd)/species_tree_ASTRAL/astral.tre" --branch-annotate 2 2> species_tree_ASTRAL/astral_logfile.log
	echo -e "Astral tree complete...   $(date)" >> $logfile 2>&1

	cd ..

	echo -e "****************************************************************" >> $logfile 2>&1
	echo -e "*** Completed Phylogenetic Analyses at $(date) ***" >> $logfile 2>&1
	echo -e "****************************************************************\n" >> $logfile 2>&1
else
	echo -e "Directory 'phylogenetic_analysis' detected in launch directory, EXITING without processing..." >> $logfile 2>&1
fi


###### DONE ######
# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished phylogenetics run at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" >> $logfile 2>&1