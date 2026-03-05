#!/bin/bash

### Parameters to set ###
# Set output folder name
outfolder="iupac_exon_alignments"
		

# Check to make sure script is being run from correct directory
if [ ! -d 2_Phased_Sequences/ ]; then
  echo "Directory '2_Phased_Sequences' does not exist! Exiting..."
  exit 1
fi

#set -eo pipefail
module load java/17.0.1
module load parallel/20210322

# Shell script to recreate the IUPAC ambiguity coded alignments for Artocarpus.
if [ ! -d 3_exon_alignments/ ]; then
	mkdir -p 3_exon_alignments
	rm -f 3_exon_alignments/*
fi

cd /data/gpfs/projects/punim1517/harvey_orel/projects/Rutaceae/Eriostemon_Group_analyses/7_allele_phasing/3_exon_alignments

genelist=/data/gpfs/projects/punim1517/harvey_orel/projects/Rutaceae/Eriostemon_Group_analyses/7_allele_phasing/genelist.txt
namelist=/data/gpfs/projects/punim1517/harvey_orel/projects/Rutaceae/Eriostemon_Group_analyses/7_allele_phasing/sample_list.txt

##########EXONS############

###### Exon sequences generated from HybPiper output:

mkdir -p exon
rm -f exon/* 
parallel "cat /data/gpfs/projects/punim1517/harvey_orel/projects/Rutaceae/Eriostemon_Group_analyses/7_allele_phasing/2_Phased_Sequences/iupac/{1}_IUPAC/{1}_exon/{2}.iupac.FNA >> exon/{2}.iupac.FNA" :::: $namelist :::: $genelist
wait

##### Alignments with MACSE

parallel --eta java -jar /data/gpfs/projects/punim1517/bin/macse_v2.06.jar -prog alignSequences -seq exon/{}.iupac.FNA :::: $genelist
wait

##### Replace frame shifts ! with gaps -

mkdir -p macse
rm -f macse/*
mv exon/*_NT* macse

parallel sed -i "s/\!/-/g" macse/{}.iupac_NT.FNA :::: $genelist
wait

##### Trim alignments to retain only sites present in 75% of taxa

mkdir -p exon_trimmed
rm -f exon_trimmed/*

parallel "trimal -gt 0.75 -in macse/{}.iupac_NT.FNA -out exon_trimmed/{}.iupac.macse.trimmed.FNA" :::: $genelist

# Move files
mkdir $outfolder
mv exon $outfolder
mv macse $outfolder
mv exon_trimmed $outfolder

for file in $outfolder/exon_trimmed/*.FNA; do generate_AstralPRO_mapping_file.py $file; done
cat $outfolder/exon_trimmed/*.txt > $outfolder/AstralPro_mapping_duplicates.txt
rm $outfolder/exon_trimmed/*_apro_mapping.txt
sort $outfolder/AstralPro_mapping_duplicates.txt | uniq > $outfolder/AstralPro_mapping_unique.txt


echo -e "\n\nDone!!!\n\n"
