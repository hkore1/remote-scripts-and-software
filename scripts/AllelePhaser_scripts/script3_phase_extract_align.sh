#!/bin/bash

# A wrapper script to perform functions in 4 subscripts that:
#		1. Phase alleles from IUPAC reference supercontigs (with WhatsHap)
#		2. Generate separate fasta files for phased sequences (with bcftools)
#		3. Replace variant sites outside the longest phase block with ambiguity codes (with haplonerate.py)
#		4. Generate separate files for the intron, exon, and supercontig sequences for samples for default HybPiper, 
#			IUPAC-coded, and phased-allele data (with intron_exon_extractor.py)

# Assumes 'script1_concat_hp2_supercontigs_and_introns_rename.sh' and 'script2_map_to_supercontigs_HO.sh' have already been run,
# and that there is a '1_New_references' folder containing IUPAC consensus/reference sequences for each sample in the launch directory.


# Check to make sure script is being run from correct directory
if [ ! -d 1_New_References/ ]; then
  echo "Directory '1_New_References' does not exist! Exiting..."
  exit 1
fi

# Specify some initial directories and variables
start_directory=$(pwd)

rm -rf 2_Phased_Sequences
mkdir 2_Phased_Sequences

## Make sure to set these!! (Path to subscript3.4 script)
Step4_subscript_path="/data/projects/punim1533/scripts/AllelePhaser_scripts/subscript3.4_intron_exon_extractor.py"



##############################
##### Leave below here #######
##############################


####### Perform Step 1 #######

echo -e "\n\n\n***********************************"
echo -e "******** Running WhatsHap *********"
echo -e "***********************************\n\n\n"

cd 1_New_References/

module load anaconda3/2021.11
eval "$(conda shell.bash hook)"
conda activate whatshap-env

for i in *; do cd $i; whatshap phase \
-o $i.supercontigs.fasta.snps.whatshap.vcf \
--no-reference \
$i.supercontigs.fasta.snps.vcf \
$i.supercontigs.fasta.marked.bam; 
whatshap stats \
--gtf $i.whatshap.gtf \
--tsv $i.whatshap.stats.tsv \
$i.supercontigs.fasta.snps.whatshap.vcf; 
cd ..;
done

conda deactivate



####### Perform Step 2 #######

echo -e "\n\n\n***********************************"
echo -e "******** Running bcftools *********"
echo -e "***********************************\n\n\n"

module load samtools/1.16.1
module load bcftools/1.15
module load parallel/20210322

set -eo pipefail

genelist_file=$1

for directory in */ ; do
    prefix="${directory///}"
    bash subscript3.2_extract_phase_bcftools.sh $prefix $genelist_file
done



####### Perform Step 3 #######

echo -e "\n\n\n***********************************"
echo -e "******* Running Haplonerate *******"
echo -e "***********************************\n\n\n"

for directory in */ ; do
	cd $directory
    prefix="${directory///}"
    echo "Haplonerating: ${prefix}"
    subscript3.3_haplonerate.py ${prefix}.whatshap.gtf phased_bcftools/${prefix}_h1.fasta phased_bcftools/${prefix}_h2.fasta --edit mask --reference ${prefix}.supercontigs.fasta.iupac.formatted.fasta --output ${prefix}.supercontigs.alleles.fasta --block ${prefix}_haplonerate_block_info.txt
    cd ..
done



####### Perform Step 4 #######

echo -e "\n\n\n*****************************************"
echo -e "***** Running Exon/Intron Extractor *****"
echo -e "*****************************************\n\n\n"

for directory in */ ; do
    prefix="${directory///}"
    mkdir ../2_Phased_Sequences/${prefix}_IUPAC
    python $Step4_subscript_path $prefix "iupac" ${start_directory}/intronerate
    mv ${prefix}/${prefix}_exon ../2_Phased_Sequences/${prefix}_IUPAC; mv ${prefix}/${prefix}_intron ../2_Phased_Sequences/${prefix}_IUPAC; mv ${prefix}/${prefix}_supercontig ../2_Phased_Sequences/${prefix}_IUPAC

    mkdir ../2_Phased_Sequences/${prefix}_ALLELES
    python $Step4_subscript_path $prefix "alleles" ${start_directory}/intronerate
    mv ${prefix}/${prefix}_exon ../2_Phased_Sequences/${prefix}_ALLELES; mv ${prefix}/${prefix}_intron ../2_Phased_Sequences/${prefix}_ALLELES; mv ${prefix}/${prefix}_supercontig ../2_Phased_Sequences/${prefix}_ALLELES

    mkdir ../2_Phased_Sequences/${prefix}_DEFAULT
    python $Step4_subscript_path $prefix "default" ${start_directory}/intronerate
    mv ${prefix}/${prefix}_exon ../2_Phased_Sequences/${prefix}_DEFAULT; mv ${prefix}/${prefix}_intron ../2_Phased_Sequences/${prefix}_DEFAULT; mv ${prefix}/${prefix}_supercontig ../2_Phased_Sequences/${prefix}_DEFAULT
done

mkdir ../2_Phased_Sequences/alleles
mkdir ../2_Phased_Sequences/default
mkdir ../2_Phased_Sequences/iupac
mv ../2_Phased_Sequences/*_ALLELES/ ../2_Phased_Sequences/alleles/
mv ../2_Phased_Sequences/*_DEFAULT/ ../2_Phased_Sequences/default/
mv ../2_Phased_Sequences/*_IUPAC/ ../2_Phased_Sequences/iupac/





