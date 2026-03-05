#!/bin/bash

# Author: Harvey Orel
# Date: Feb 19 2023

# Concatenates hp2 supercontig and intronerate files output from Theo's hp2.py wrapper script in preparation for
# allele phasing as outlined in Kates et al. 2018 (https://github.com/mossmatters/phyloscripts/tree/master/alleles_workflow)

# Usage: 'bash concat_hp2_supercontigs_and_introns_rename.sh <hp2_results_directory>'
#           where the directory is the intronerated results/ output of hp2.py containing tar files.

# The script concatenates files and changes names in the supercontig.fasta output to the format '>GENE_SAMPLE'.
# Intronerate files are renamed so that the seqname column of the .gff is the gene name.
# Concatenated files are deposited in 'supercontigs' and 'intronerate' folders in the working directory.

# Check for input directory argument
if [ $# -eq 0 ]
then
  echo "Usage: $0 <input_directory>"
  exit 1
fi

# Create output directories
mkdir -p supercontigs
mkdir -p intronerate

# Iterate over all directories in input directory
for dir in "$1"/*/
do
  # Extract prefix from directory name
  prefix=$(basename "${dir}")

  # Print progress
  echo -e "Working on sample: $prefix"

  # Extract tar file
  tar -xf "${dir%/}/$prefix.tar.gz"

  # Concatenate supercontig sequences into single fasta file
  for file in "$prefix"/*/"$prefix"/sequences/intron/*_supercontig.fasta
  do
    filename=$(basename "$file")
    filename=${filename/_supercontig.fasta/_$prefix}
    # Pipe the output to grep to remove lines containing the phrase "Joins between unique SPAdes contigs are separated by 10 \"N\" characters"
    cat "$file" | grep -v "Joins between unique SPAdes contigs are separated by 10 \"N\" characters" | awk -v fname="$filename" 'BEGIN {print ">"fname} {print}' >> "$prefix".supercontigs.fasta
  done

  # Concatenate intronerate files into single fasta file
  for file in "$prefix"/*/"$prefix"/intronerate/intronerate.gff
  do
    filename=$(dirname "$file")
    filename=$(echo $filename | cut -d'/' -f2)
    cat "$file" | sed "s/$prefix/$filename/g" >> "$prefix".intronerate.fasta
  done

  # Rename sequences in supercontig files so that the locus name and sample name are separated by "---"
  prefix_init_sep="_$prefix"
  prefix_new_sep="---$prefix"

  if [[ "$OSTYPE" == "darwin"* ]]; then
  sed -i '' -e "s/${prefix_init_sep}/${prefix_new_sep}/g" "$prefix".supercontigs.fasta
  else
  sed -i -e "s/${prefix_init_sep}/${prefix_new_sep}/g" "$prefix".supercontigs.fasta
  fi

  # Move output files to respective output directories
  mv "$prefix".supercontigs.fasta supercontigs/
  mv "$prefix".intronerate.fasta intronerate/

  # Remove extracted directory
  rm -r "$prefix"

  # Increment counter
  ((count++))
done
