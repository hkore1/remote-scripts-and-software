# Run to activate conda whatshap environment on SPARTAN (can be run on an interactive node 
#for an entire runthrough)

module load anaconda3/2021.11
eval "$(conda shell.bash hook)"
conda activate whatshap-env

# Main
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

# This is a Bash shell script that iterates over a list of directories specified 
# in a file named "namelist.txt". For each directory in the list, it changes the current 
# working directory to that directory using the "cd" command. Then, it runs two "whatshap" commands.

# The first "whatshap" command is "phase", which takes as input a VCF file containing SNPs, 
# and a BAM file containing aligned reads, and produces a phased VCF file as output. The phased VCF 
# file is named after the input file, with the suffix ".supercontigs.fasta.snps.whatshap.vcf".

# The second "whatshap" command is "stats", which takes as input the phased VCF file produced 
# by the previous command, and produces a statistics file in both TSV and GTF formats. The statistics 
# files are named after the input file, with the suffixes ".whatshap.stats.tsv" and ".whatshap.gtf", respectively.

# After running the two "whatshap" commands, the script changes the current working 
# directory back to the parent directory using the "cd" command.

# From https://github.com/mossmatters/phyloscripts/tree/master/alleles_workflow