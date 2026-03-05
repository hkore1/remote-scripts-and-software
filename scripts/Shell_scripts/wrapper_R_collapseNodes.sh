# Load required modules
module load gcccore/8.3.0
module load r/4.1.0

while getopts "i:o:" opt
do
   case "$opt" in
      i ) input_trees="$OPTARG" ;;
      o ) output_trees="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

Rscript /data/gpfs/projects/punim1517/scripts/R_collapseNodes.R $input_trees $output_trees