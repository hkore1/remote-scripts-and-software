#!/bin/bash

# Author: Harvey Orel
# Date: Feb 19 2023

# This workflow will take the concatenated (using 'concat_hp2_supercontigs_and_introns.sh') supercontig 
# output of HybPiper and return a supercontig that contains heterozygous positions as 
# ambiguity bases. Uses paired reads. Modified from 'map_to_supercontigs.sh' script from 
# https://github.com/mossmatters/phyloscripts/tree/master/alleles_workflow, and employed in Kates et al. 2018

# The script should be run on a directory containing the supercontig files and a directory containing 
# the corresponding read files.

if [[ $# -ne 3 ]] ; then
    echo 'usage: script2_map_to_supercontigs.sh supercontig_directory read_directory number_of_CPUs'
    exit 1
fi

supercontig_dir=$1
read_dir=$2
num_CPUs=$3
CPUsMinusOne=$((num_CPUs-1))

#########CHANGE THESE PATHS AS NEEDED###########

gatkpath=/data/gpfs/projects/punim1517/bin/gatk/gatk
picardpath=/data/gpfs/projects/punim1517/bin/picard/build/libs/picard.jar

#############COMMAND LINE ARGUMENTS############

module load samtools/1.16.1
module load java/17.0.1

output_dir="1_New_References"
mkdir -p $output_dir

# Create stored copy of the supercontigs folder, as these files are moved out of the folder during the below for loop.
cp -r $supercontig_dir/* supercontigs_stored/

# Initialise for loop and specify input variables
for supercontig_file in "$supercontig_dir"/*.supercontigs.fasta; do
    filename=$(basename "$supercontig_file" .supercontig.fasta)
    prefix=${filename%%.*}
    read1suffix="R1.fastq.gz"
    read2suffix="R2.fastq.gz"
    read1fq="$(find $read_dir -name "${prefix}_*${read1suffix}")"
    read2fq="$(find $read_dir -name "${prefix}_*${read2suffix}")"

    supercontig_filename=$prefix.supercontigs.fasta
    supercontig_path=$(dirname "$supercontig_file")
    supercontig="$(pwd)/${supercontig_path}/${supercontig_filename}"

    echo -e "\n\n*************************"
    echo -e "Running $prefix"
    echo -e "*************************\n"

    echo -e "Using supercontig file: $supercontig"
    echo -e "Using read file 1: $read1fq"
    echo -e "Using read file 2: $read2fq\n"

    mkdir $prefix
    SECONDS=0

    #####STEP ZERO: Make Reference Databases

    java -jar $picardpath CreateSequenceDictionary \
    R=$supercontig
    bwa index $supercontig 
    samtools faidx $supercontig
    wait


    #####STEP ONE: Map reads

    echo "Mapping Reads"

    bwa mem $supercontig $read1fq $read2fq -t num_CPUs| samtools view -bS - -@ CPUsMinusOne| samtools sort - -o $supercontig.sorted.bam -@ CPUsMinusOne

    java -jar $picardpath FastqToSam  \
    F1=$read1fq \
    F2=$read2fq \
    O=$supercontig.unmapped.bam \
    SM=$supercontig

    java -jar $picardpath MergeBamAlignment \
    ALIGNED=$supercontig.sorted.bam \
    UNMAPPED=$supercontig.unmapped.bam \
    O=$supercontig.merged.bam \
    R=$supercontig


    #####STEP TWO: Mark duplicates

    echo "Marking Duplicates"
    java -jar $picardpath MarkDuplicates \
    I=$supercontig.merged.bam \
    O=$supercontig.marked.bam \
    M=$supercontig.metrics.txt


    #######STEP THREE: Identify variants, select only SNPs

    echo "Identifying variants"

    samtools index $supercontig.marked.bam

    $gatkpath HaplotypeCaller \
    -R $supercontig \
    -I $supercontig.marked.bam \
    -O $supercontig.vcf

    time $gatkpath SelectVariants \
    -R $supercontig \
    -V $supercontig.vcf \
    --select-type-to-include SNP \
    -O $supercontig.snps.vcf 


    ######STEP FOUR: Output new supercontig FASTA with ambiguity codes

    echo "Generating IUPAC FASTA file"

    $gatkpath FastaAlternateReferenceMaker \
    -R $supercontig \
    -O $supercontig.iupac.fasta \
    -V $supercontig.snps.vcf \
    --use-iupac-sample $supercontig


    # Generate another IUPAC fasta file in the correct format for downstream processes
    cp $supercontig.iupac.fasta $supercontig.iupac.formatted.fasta
    sed -i 's/^>[^ ]*/>/g' $supercontig.iupac.formatted.fasta
    sed -i 's/^>[ ]/>/g' $supercontig.iupac.formatted.fasta
    sed -i 's/:.*//g' $supercontig.iupac.formatted.fasta

    wait
    mv "$(pwd)"/"$supercontig_path"/"${prefix}"* $prefix
    mv $prefix $output_dir
    t=$SECONDS
    echo -e '\n\n______________________'
    echo -e 'Time taken:' "$(( t/60 - 1440*(t/86400) ))" 'minutes'
    echo -e '______________________\n'

done

