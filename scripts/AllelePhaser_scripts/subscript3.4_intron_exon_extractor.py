#!/usr/bin/env python

# Script to use the GFF files from intronerate and the ambiguity-encoded FASTA files to generate separate intron and exon files for each gene.

import re,sys,os,errno,shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

my_re = re.compile(r"([0-9]+ )(.+):1")

#Fix the names in IUPAC file

prefix=sys.argv[1]
fileType=sys.argv[2]
intronerate_folder=sys.argv[3]

if os.path.isdir(prefix):
    os.chdir(prefix)
        
print("Working on sample: " + prefix)
#Parse GFF into dictionaries for each gene (one for introns, one for exons)
# ASSUMES THE GFF IS SORTED WITHIN EACH GENE!!!!

intron_dict = {}
exon_dict = {}

gff_fn = intronerate_folder+"/"+prefix+".intronerate.fasta" #sys.argv[3]
for line in open(gff_fn):
    line=line.split()
    try:
        if line[2] == "exon":
            try:
                exon_dict[line[0]].append((int(line[3])-1,int(line[4])))
            except KeyError:
                exon_dict[line[0]] = [(int(line[3])-1,int(line[4]))]
        elif line[2] == "intron":
            try:
                intron_dict[line[0]].append((int(line[3])-1,int(line[4])))
            except KeyError:
                intron_dict[line[0]] = [(int(line[3])-1,int(line[4]))]
    except IndexError:
        pass
exon_dict = {f'{k}---'+prefix: v for k, v in exon_dict.items()} # reformat dict to same style as supercontig dict



if fileType == "iupac":
    try:
        supercontig_dict = SeqIO.to_dict(SeqIO.parse("{}.supercontigs.fasta.iupac.formatted.fasta".format(prefix),'fasta'))
        dataType = "iupac"
        print("Generating IUPAC sequences...")
    except FileNotFoundError:
        print("Error in 'iupac' option: SAMPLE.supercontigs.fasta.iupac.formatted.fasta file not found --- EXITING")
        sys.exit()
elif fileType == "alleles":
    try:
        supercontig_dict = SeqIO.to_dict(SeqIO.parse("{}.supercontigs.alleles.fasta".format(prefix),'fasta'))
        dataType = "alleles"
        print("Generating phased allele sequences...")
    except FileNotFoundError:
        print("Error in 'alleles' option: SAMPLE.supercontigs.alleles.fasta file not found --- EXITING")
        sys.exit()
elif fileType == "default":
    try:
        supercontig_dict = SeqIO.to_dict(SeqIO.parse("{}.supercontigs.fasta".format(prefix),'fasta'))
        dataType = "default"
        print("Generating default sequences...")
    except FileNotFoundError:
        print("Error in 'default' option: SAMPLE.supercontigs.fasta file not found --- EXITING")
        sys.exit()
else:
    print("Error: No valid option selected for input data type --- EXITING")
    sys.exit()



for gene in exon_dict:
    try:
        geneLength = len(supercontig_dict[gene])
    except KeyError:
        haploGeneName = gene+"_h1"
        geneLength = len(supercontig_dict[haploGeneName])
    exon_ranges = exon_dict[gene]
    for exon_interval in range(len(exon_ranges)+1):
        if exon_interval == 0:
            intron_dict[gene] = [(0,exon_ranges[exon_interval][0]-1)]
        elif exon_interval == len(exon_ranges)  :
            intron_dict[gene].append((exon_ranges[-1][1],geneLength))

        else:
            start = exon_ranges[exon_interval - 1][1] 
            stop = exon_ranges[exon_interval][0] - 1
            intron_dict[gene].append((start,stop))



newseq = ''

for seqType in [prefix+"_exon",prefix+"_intron",prefix+"_supercontig"]:
    if os.path.exists(seqType):
        shutil.rmtree(seqType)
    os.makedirs(seqType)

for gene in supercontig_dict:
    geneName = gene.split("---")[0]
    sampleName = gene.split("---")[1]
            
    with open(prefix+"_exon/{}.{}.FNA".format(geneName,dataType),'a') as exonout:
        newseq = ''
        exonLookupName = supercontig_dict[gene].id.replace("_h1",'')
        exonLookupName = exonLookupName.replace("_h2",'')
        if exonLookupName not in exon_dict:
            continue
        for gff_interval in exon_dict[exonLookupName]:
            newseq += supercontig_dict[gene].seq[gff_interval[0]:gff_interval[1]]
        exonout.write(">{}\n{}\n".format(sampleName,newseq))
    ### Block to add '_h1' and '_h2' suffix to sequences in the fasta if dataType == 'alleles'
    # if dataType == "alleles":    
    #     with open(prefix+"_exon/{}.{}.FNA".format(geneName,dataType), 'r+') as fp:
    #         content=fp.readlines()
    #         line_string=[]
    #         x = len(content)
    #         if x == 4:
    #             for line in content:
    #                 line_string.append(line)
    #             line_string[0] = line_string[0].replace("\n","_h1\n")
    #             line_string[2] = line_string[2].replace("\n","_h2\n")
    #             content[0] = line_string[0]
    #             content[2] = line_string[2]
    #             with open(prefix+"_exon/{}.{}.FNA".format(geneName,dataType), 'w') as fp:
    #                 for line in content:
    #                     fp.write(line)

    with open(prefix+"_intron/{}.intron.{}.fasta".format(geneName,dataType),'a') as intronout:
        newseq=''
        for gff_interval in intron_dict[exonLookupName]:
            newseq += supercontig_dict[gene].seq[gff_interval[0]:gff_interval[1]]
        intronout.write(">{}\n{}\n".format(sampleName,newseq))
    ### Block to add '_h1' and '_h2' suffix to sequences in the fasta if dataType == 'alleles'
    # if dataType == "alleles":    
    #     with open(prefix+"_intron/{}.intron.{}.fasta".format(geneName,dataType), 'r+') as fp:
    #         content=fp.readlines()
    #         line_string=[]
    #         x = len(content)
    #         if x == 4:
    #             for line in content:
    #                 line_string.append(line)
    #             line_string[0] = line_string[0].replace("\n","_h1\n")
    #             line_string[2] = line_string[2].replace("\n","_h2\n")
    #             content[0] = line_string[0]
    #             content[2] = line_string[2]
    #             with open(prefix+"_intron/{}.intron.{}.fasta".format(geneName,dataType), 'w') as fp:
    #                 for line in content:
    #                     fp.write(line)
    
    with open(prefix+"_supercontig/{}.supercontig.{}.fasta".format(geneName,dataType),'a') as supercontigout:
        supercontigout.write(">{}\n{}\n".format(sampleName,supercontig_dict[gene].seq))
    ### Block to add '_h1' and '_h2' suffix to sequences in the fasta if dataType == 'alleles'
    # if dataType == "alleles":    
    #     with open(prefix+"_supercontig/{}.supercontig.{}.fasta".format(geneName,dataType), 'r+') as fp:
    #         content=fp.readlines()
    #         line_string=[]
    #         x = len(content)
    #         if x == 4:
    #             for line in content:
    #                 line_string.append(line)
    #             line_string[0] = line_string[0].replace("\n","_h1\n")
    #             line_string[2] = line_string[2].replace("\n","_h2\n")
    #             content[0] = line_string[0]
    #             content[2] = line_string[2]
    #             with open(prefix+"_supercontig/{}.supercontig.{}.fasta".format(geneName,dataType), 'w') as fp:
    #                 for line in content:
    #                     fp.write(line)
        
