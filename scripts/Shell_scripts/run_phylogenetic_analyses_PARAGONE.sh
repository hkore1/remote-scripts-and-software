#!/bin/bash

################################################################
# Run phylogenetic analyses on ParaGone outputs
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
#SBATCH --account="punim1533"

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
#SBATCH --time=3-00:0:00

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
export MO_INPUT="/data/gpfs/projects/punim1533/corymbia_phylogeny/6_analyses/nDNA/3_ParaGone/paragone_original_MO/26_MO_final_alignments_trimmed"
export RT_INPUT="/data/gpfs/projects/punim1533/corymbia_phylogeny/6_analyses/nDNA/3_ParaGone/paragone_original_MO/28_RT_final_alignments_trimmed"
export MI_INPUT="/data/gpfs/projects/punim1533/corymbia_phylogeny/6_analyses/nDNA/3_ParaGone/paragone_original_MO/27_MI_final_alignments_trimmed"



########################
### Leave below here ###
########################

# set up a master log file
logfile="$SLURM_SUBMIT_DIR/phylo_$SLURM_JOB_ID.log"

# load modules
module load java/11.0.2 

# define alignment base directories
export MO_BASE=$(basename $MO_INPUT)
export RT_BASE=$(basename $RT_INPUT)
export MI_BASE=$(basename $MI_INPUT)

# Record and echo start time and input
start="$(date +%s)"
echo -e "****************************************************************" >> $logfile 2>&1
echo -e "*** Starting phylogenetics run at $(date) ***" >> $logfile 2>&1
echo -e "****************************************************************" >> $logfile 2>&1
echo -e "Provided path to MO alignments: $MO_INPUT" >> $logfile 2>&1
echo -e "Provided path to RT alignments: $RT_INPUT" >> $logfile 2>&1
echo -e "Provided path to MI alignments: $MI_INPUT\n" >> $logfile 2>&1

# Copy script to submit directory for future reference of employed settings
\cp /data/gpfs/projects/punim1533/scripts/run_phylogenetic_analyses_PARAGONE.sh $SLURM_SUBMIT_DIR; mv run_phylogenetic_analyses.sh logged_script_used.sh

#################### Run analyses on MO ####################
if [ ! -d "$SLURM_SUBMIT_DIR/phylogenetics_MO" ]; then
	echo -e "****************************************************************" >> $logfile 2>&1
	echo -e "*** Beginning Monophyletic Outgroups at $(date) ***" >> $logfile 2>&1
	echo -e "****************************************************************" >> $logfile 2>&1

	rm -r phylogenetics_MO; mkdir phylogenetics_MO; cd phylogenetics_MO/

	## Run a concatenation tree using IQTREE
	rm -r concatenated_IQTREE; mkdir concatenated_IQTREE; cd concatenated_IQTREE/
	cp -R $MO_INPUT $PWD
	iqtree2 --prefix concat -p $MO_BASE --threads-max "48" -T AUTO -B 1000
	rm -r $MO_BASE; cd ..
	echo -e "Concatenated tree complete...   $(date)" >> $logfile 2>&1

	## Run gene trees for each alignment using IQTREE
	rm -r genetrees_IQTREE; mkdir genetrees_IQTREE; cd genetrees_IQTREE/
	cp -R $MO_INPUT $PWD
	run_iqtree_in_parallel.py -af $MO_BASE -s "-B 1000" -c 8 -t 6
	wait; cd $MO_BASE
	cat *.treefile > loci.treefile; mv loci.treefile ../
	cat *.log > loci.log; mv loci.log ../
	cat *.bionj > loci.bionj; mv loci.bionj ../
	cat *.iqtree > loci.iqtree; mv loci.iqtree ../
	cat *.mldist > loci.mldist; mv loci.mldist ../
	cd ..; rm -r $MO_BASE; cd ..
	echo -e "Individual gene trees complete...   $(date)" >> $logfile 2>&1

	## Calculate concordance factors on the concatenation tree
	rm -r concordance_IQTREE; mkdir concordance_IQTREE; cd concordance_IQTREE/
	cp -R $MO_INPUT $PWD
	cp ../concatenated_IQTREE/concat.treefile $PWD
	cp ../genetrees_IQTREE/loci.treefile $PWD
	iqtree2 -t concat.treefile --gcf loci.treefile -p $MO_BASE --scfl 1000 --prefix concord
	rm -r $MO_BASE; cd ..
	echo -e "Concordance factors complete...   $(date)" >> $logfile 2>&1

	## Run an ASTRAL tree using the gene trees
	rm -r species_tree_ASTRAL; mkdir species_tree_ASTRAL
	java -Xmx8000M -jar /data/projects/punim1517/bin/Astral/astral.5.7.8.jar -i "$(pwd)/genetrees_IQTREE/loci.treefile" -o "$(pwd)/species_tree_ASTRAL/astral.tre" --branch-annotate 2 2> species_tree_ASTRAL/astral_logfile.log
	echo -e "Astral tree complete...   $(date)" >> $logfile 2>&1

	cd ..

	echo -e "****************************************************************" >> $logfile 2>&1
	echo -e "*** Completed Monophyletic Outgroups at $(date) ***" >> $logfile 2>&1
	echo -e "****************************************************************\n" >> $logfile 2>&1
else
	echo -e "Directory 'phylogenetics_MO' detected in launch directory, skipping MO" >> $logfile 2>&1
fi


#################### Run analyses on RT ####################
if [ ! -d "$SLURM_SUBMIT_DIR/phylogenetics_RT" ]; then
	echo -e "****************************************************************" >> $logfile 2>&1
	echo -e "*** Beginning Rooted Subtrees at $(date) ***" >> $logfile 2>&1
	echo -e "****************************************************************" >> $logfile 2>&1

	rm -r phylogenetics_RT; mkdir phylogenetics_RT; cd phylogenetics_RT/

	## Run a concatenation tree using IQTREE
	rm -r concatenated_IQTREE; mkdir concatenated_IQTREE; cd concatenated_IQTREE/
	cp -R $RT_INPUT $PWD
	iqtree2 --prefix concat -p $RT_BASE --threads-max "48" -T AUTO -B 1000
	rm -r $RT_BASE; cd ..
	echo -e "Concatenated tree complete...   $(date)" >> $logfile 2>&1

	## Run gene trees for each alignment using IQTREE -- This is different to others due to iqtree crashes
	rm -r genetrees_IQTREE; mkdir genetrees_IQTREE; cd genetrees_IQTREE/
	cp -R $RT_INPUT $PWD
	run_iqtree_in_parallel.py -af $RT_BASE -s "-B 1000" -c 8 -t 6
	wait; cd $RT_BASE
	cat *.treefile > loci.treefile; mv loci.treefile ../
	cat *.log > loci.log; mv loci.log ../
	cat *.bionj > loci.bionj; mv loci.bionj ../
	cat *.iqtree > loci.iqtree; mv loci.iqtree ../
	cat *.mldist > loci.mldist; mv loci.mldist ../
	cd ..; rm -r $RT_BASE; cd ..
	echo -e "Individual gene trees complete...   $(date)" >> $logfile 2>&1

	## Calculate concordance factors on the concatenation tree
	rm -r concordance_IQTREE; mkdir concordance_IQTREE; cd concordance_IQTREE/
	cp -R $RT_INPUT $PWD
	cp ../concatenated_IQTREE/concat.treefile $PWD
	cp ../genetrees_IQTREE/loci.treefile $PWD
	iqtree2 -t concat.treefile --gcf loci.treefile -p $RT_BASE --scfl 1000 --prefix concord
	rm -r $RT_BASE; cd ..
	echo -e "Concordance factors complete...   $(date)" >> $logfile 2>&1

	## Run an ASTRAL tree using the gene trees
	rm -r species_tree_ASTRAL; mkdir species_tree_ASTRAL
	java -Xmx8000M -jar /data/projects/punim1517/bin/Astral/astral.5.7.8.jar -i "$(pwd)/genetrees_IQTREE/loci.treefile" -o "$(pwd)/species_tree_ASTRAL/astral.tre" --branch-annotate 2 2> species_tree_ASTRAL/astral_logfile.log
	echo -e "Astral tree complete...   $(date)" >> $logfile 2>&1

	cd ..

	echo -e "****************************************************************" >> $logfile 2>&1
	echo -e "*** Completed Rooted Subtrees at $(date) ***" >> $logfile 2>&1
	echo -e "****************************************************************\n" >> $logfile 2>&1
else
	echo -e "Directory 'phylogenetics_RT' detected in launch directory, skipping RT" >> $logfile 2>&1
fi


#################### Run analyses on MI ####################
if [ ! -d "$SLURM_SUBMIT_DIR/phylogenetics_MI" ]; then
	echo -e "****************************************************************" >> $logfile 2>&1
	echo -e "*** Beginning Monophyletic Ingroups at $(date) ***" >> $logfile 2>&1
	echo -e "****************************************************************" >> $logfile 2>&1

	rm -r phylogenetics_MI; mkdir phylogenetics_MI; cd phylogenetics_MI/

	## Run a concatenation tree using IQTREE
	rm -r concatenated_IQTREE; mkdir concatenated_IQTREE; cd concatenated_IQTREE/
	cp -R $MI_INPUT $PWD
	iqtree2 --prefix concat -p $MI_BASE --threads-max "48" -T AUTO -B 1000
	rm -r $MI_BASE; cd ..
	echo -e "Concatenated tree complete...   $(date)" >> $logfile 2>&1

	## Run gene trees for each alignment using IQTREE -- This is different to others due to iqtree crashes
	rm -r genetrees_IQTREE; mkdir genetrees_IQTREE; cd genetrees_IQTREE/
	cp -R $MI_INPUT $PWD
	run_iqtree_in_parallel.py -af $MI_BASE -s "-B 1000" -c 8 -t 6
	wait; cd $MI_BASE
	cat *.treefile > loci.treefile; mv loci.treefile ../
	cat *.log > loci.log; mv loci.log ../
	cat *.bionj > loci.bionj; mv loci.bionj ../
	cat *.iqtree > loci.iqtree; mv loci.iqtree ../
	cat *.mldist > loci.mldist; mv loci.mldist ../
	cd ..; rm -r $MI_BASE; cd ..
	echo -e "Individual gene trees complete...   $(date)" >> $logfile 2>&1

	## Calculate concordance factors on the concatenation tree
	rm -r concordance_IQTREE; mkdir concordance_IQTREE; cd concordance_IQTREE/
	cp -R $MI_INPUT $PWD
	cp ../concatenated_IQTREE/concat.treefile $PWD
	cp ../genetrees_IQTREE/loci.treefile $PWD
	iqtree2 -t concat.treefile --gcf loci.treefile -p $MI_BASE --scf 100 --prefix concord
	rm -r $MI_BASE; cd ..
	echo -e "Concordance factors complete...   $(date)" >> $logfile 2>&1

	## Run an ASTRAL tree using the gene trees
	rm -r species_tree_ASTRAL; mkdir species_tree_ASTRAL
	java -Xmx8000M -jar /data/projects/punim1517/bin/Astral/astral.5.7.8.jar -i "$(pwd)/genetrees_IQTREE/loci.treefile" -o "$(pwd)/species_tree_ASTRAL/astral.tre" --branch-annotate 2 2> species_tree_ASTRAL/astral_logfile.log
	echo -e "Astral tree complete...   $(date)" >> $logfile 2>&1

	cd ..

	echo -e "****************************************************************" >> $logfile 2>&1
	echo -e "*** Completed Monophyletic Ingroups at $(date) ***" >> $logfile 2>&1
	echo -e "****************************************************************" >> $logfile 2>&1
else
	echo -e "Directory 'phylogenetics_MI' detected in launch directory, skipping MI" >> $logfile 2>&1
fi


###### DONE ######
# Record and echo end time and duration
end="$(date +%s)"
duration="$(( $end - $start ))"
duration_mins=$(echo "scale=2; ${duration}/60" | bc)
duration_hours=$(echo "scale=2; ${duration}/3600" | bc)

echo -e "\nFinished phylogenetics run at $(date) after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours" >> $logfile 2>&1