#!/bin/bash


# Author: Harvey Orel (modified from Matt Johnson - https://github.com/mossmatters/phyloscripts/blob/master/alleles_workflow/create_alleles_alignments.sh)
# Description: Script to generate cleaned/trimmed phased exon alignments from allele sequences generated in scripts 1-3.

# Usage: Run from directory with '2_Phased_Sequences/' folder in it.


### Parameters to set ###
# Set output folder name
outfolder="phased_allele_exon_alignments"

# List of genes and samples
genelist=/data/gpfs/projects/punim1533/corymbia_phylogeny/6_analyses/nDNA/4_allele_phasing/genelist.txt
namelist=/data/gpfs/projects/punim1533/corymbia_phylogeny/6_analyses/nDNA/4_allele_phasing/_sample_list.txt
		
# Check to make sure script is being run from correct directory
if [ ! -d 2_Phased_Sequences/ ]; then
  echo "Directory '2_Phased_Sequences' does not exist! Exiting..."
  exit 1
fi

#set -eo pipefail
module load java/17.0.1
module load parallel/20210322

# Create new directory.
if [ ! -d 3_exon_alignments/ ]; then
	mkdir -p 3_exon_alignments
	rm -f 3_exon_alignments/*
fi

#cd /data/gpfs/projects/punim1517/harvey_orel/projects/Rutaceae/Eriostemon_Group_analyses/7_allele_phasing/3_exon_alignments
cd 3_exon_alignments # This should work but is untested, if failing here provide full path as above.

##########EXONS############

###### Exon sequences generated from HybPiper output:

mkdir -p exon
rm -f exon/* 
parallel --colsep='\n' "cat /data/gpfs/projects/punim1533/corymbia_phylogeny/6_analyses/nDNA/4_allele_phasing/2_Phased_Sequences/alleles/{1}_ALLELES/{1}_exon/{2}.alleles.FNA >> exon/{2}.alleles.FNA" :::: $namelist :::: $genelist
wait

##### Alignments with MACSE

parallel --eta java -jar /data/gpfs/projects/punim1533/bin/macse_v2.06.jar -prog alignSequences -seq exon/{}.alleles.FNA :::: $genelist
wait

##### Replace frame shifts ! with gaps -

mkdir -p macse
rm -f macse/*
mv exon/*_NT* macse

parallel sed -i "s/\!/-/g" macse/{}.alleles_NT.FNA :::: $genelist
wait

##### Trim alignments to retain only sites present in 75% of taxa

mkdir -p exon_trimmed
rm -f exon_trimmed/*

parallel "trimal -gt 0.75 -in macse/{}.alleles_NT.FNA -out exon_trimmed/{}.alleles.macse.trimmed.FNA" :::: $genelist

# Move files
mkdir $outfolder
mv exon $outfolder
mv macse $outfolder
mv exon_trimmed $outfolder

# Generate mapping file for Astral-Pro (relies on the python script 'generate_AstralPRO_mapping_file.py', from https://github.com/hkore1/python_scripts/blob/main/generate_AstralPRO_mapping_file.py)
for file in $outfolder/exon_trimmed/*.FNA; do generate_AstralPRO_mapping_file.py $file; done
cat $outfolder/exon_trimmed/*.txt > $outfolder/AstralPro_mapping_duplicates.txt
rm $outfolder/exon_trimmed/*_apro_mapping.txt
sort $outfolder/AstralPro_mapping_duplicates.txt | uniq > $outfolder/AstralPro_mapping_unique.txt


echo -e "\n\nDone!!!\n\n"
