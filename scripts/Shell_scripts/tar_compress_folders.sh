#!/bin/bash


# Set Slurm directives

#SBATCH --partition=cascade
#SBATCH --nodes=1
#SBATCH --job-name="compress_folder"
#SBATCH --account="punim1533"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=96G
#SBATCH --time=2-00:0:00

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

module load GCCcore/11.3.0
module load pigz/2.7

tar --remove-files --use-compress-program=pigz -cvf 1_New_References.tar.gz 1_New_References/

echo -e "\nFinished compressing folder after running for $duration seconds, or" \
	"$duration_mins minutes, or $duration_hours hours"
