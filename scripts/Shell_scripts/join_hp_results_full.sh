#!/bin/bash

######################
# Author: Harvey Orel
# Date: 22 Nov 2022
# Description: Join Hybpiper2 results using a modified version of Theo's script, and then perform some
# 				additional steps to neaten joined files.
# Note: Actions - moves joined results to a new folder, plots gene recovery and paralog count heatmaps.
# Usage: join_hp_results_full.sh -r results/ -n 6 -i y
######################

### Check args

helpFunction()
{
   echo ""
   echo "Usage: $0 -r hybpiper_results_dir -n numCPUs -i join_supercontigs"
   echo -e "\t-r Path to hybpiper results/"
   echo -e "\t-n Number of cores to use"
   echo -e "\t-i Join intronerated supercontigs ('y' or 'n')"
   exit 1 # Exit script after printing help
}

while getopts "r:n:i:" opt
do
   case "$opt" in
      r ) hybpiper_results_dir="$OPTARG" ;;
      n ) numCPUs="$OPTARG" ;;
      i ) join_supercontigs="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$hybpiper_results_dir" ] || [ -z "$numCPUs" ] || [ -z "$join_supercontigs" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Begin script in case all parameters are correct
echo "Commencing with following settings"
echo "Results Directory: $hybpiper_results_dir"
echo "Number of cores: $numCPUs"
echo -e "Joining intronerated supercontigs: $join_supercontigs\n"



# Load modules
#module load gcccore/8.3.0
module load GCC/11.3.0
#module load r/4.0.0
module load R/4.2.2

# Run the modified join function
join_hp_results_mod.py $hybpiper_results_dir $numCPUs $join_supercontigs

# create dir for joined results
mkdir combined_results

# move joined results to new dir
mv results/all_lengths.txt combined_results/
mv results/all_paralogs_joined combined_results/
mv results/compiled_paralogs.txt combined_results/
mv results/paralogs_no_chimeras_joined combined_results/
mv results/stats_joined.txt combined_results/
mv results/supercontig_seqs_joined combined_results/
mv results/intron_seqs_joined combined_results/

# Run scripts to create heatmaps
Rscript /data/projects/punim1533/scripts/R_CJJ_gene_recovery_heatmap_ggplot.R combined_results/all_lengths.txt
Rscript /data/projects/punim1533/scripts/R_sample_histogram_ggplot.R combined_results/all_lengths.txt
Rscript /data/projects/punim1533/scripts/R_plot_compiled_paralogs.R combined_results/compiled_paralogs.txt

rm Rplots.pdf
