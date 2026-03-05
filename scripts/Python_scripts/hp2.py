#!/usr/bin/env python3
import sys
import re
import glob
import subprocess as sp
import os
import concurrent.futures

# Theo Allnutt 2022
# Runs HybPiper 2.0 on master read directory - reads must already be filtered and concatenated if multiple runs exist for each file
# Usage:
# hp2.py reads_dir/ result_dir/ baitfile aa username "--diamond"
# where aa=bait type (or dna); username=your nci username; in quotes=other hybpiper assemble options, to use blastx and defaults leave empty:""

# Modified: Harvey Orel 8 May 2023, line 127 - to take prefix from start of index rather than first underscore.

digits = re.compile(r'(\d+)')
def tokenize(filename):
    return tuple(int(token) if match else token
                 for token, match in
                 ((fragment, digits.search(fragment))
                  for fragment in digits.split(filename)))

def run_hp(i):

	#runs all hp commands, 8 samples at a time (optimum for Gadi)
	
	print("Assembling sample",i)
	
	p0=sp.Popen("rm -f logs/%s.out;rm -f logs/%s.err" %(i,i),shell=True).wait()
	
	if bait_type=="aa":
		p1=sp.Popen("%s/.local/bin/hybpiper assemble --targetfile_aa %s -r %s/%s_* --prefix %s --cpu %s %s >>logs/%s.out 2>>logs/%s.err" %(username,baits,reads_dir,i,i,cputhreads,ass_options,i,i),shell=True).wait()
	else:
		p1=sp.Popen("%s/.local/bin/hybpiper assemble --targetfile_dna %s -r %s/%s_* --prefix %s --cpu %s %s >>logs/%s.out 2>>logs/%s.err" %(username,baits,reads_dir,i,i,cputhreads,ass_options,i,i),shell=True).wait()		
	
	
	#make namelist for stats
	g1=open(i+".name",'w')
	g1.write(i+"\n")
	g1.close()
	
	#make dir for results
	p2=sp.Popen("rm -rf results/%s;mkdir results/%s" %(i,i),shell=True).wait()
	
	print("Stats for sample",i)
	
	if bait_type=="aa":
		p3=sp.Popen("%s/.local/bin/hybpiper stats --targetfile_aa %s gene %s.name --seq_lengths_filename results/%s/%s_lengths --stats_filename results/%s/%s_stats >>logs/%s.out 2>>logs/%s.err" %(username,baits,i,i,i,i,i,i,i),shell=True).wait()
	else:
		p3=sp.Popen("%s/.local/bin/hybpiper stats --targetfile_dna %s gene %s.name --seq_lengths_filename results/%s/%s_lengths --stats_filename results/%s/%s_stats >>logs/%s.out 2>>logs/%s.err" %(username,baits,i,i,i,i,i,i,i),shell=True).wait()
	
	
	print("Retrieving sequences for sample",i)
	
	if bait_type=="aa":
	
		p4=sp.Popen("%s/.local/bin/hybpiper retrieve_sequences  --targetfile_aa %s dna --sample_names %s.name --fasta_dir results/%s/%s_dna_seqs >>logs/%s.out 2>>logs/%s.err" %(username,baits,i,i,i,i,i),shell=True).wait()
		
		p5=sp.Popen("%s/.local/bin/hybpiper retrieve_sequences --targetfile_aa %s aa --sample_names %s.name --fasta_dir results/%s/%s_aa_seqs >>logs/%s.out 2>>logs/%s.err" %(username,baits,i,i,i,i,i),shell=True).wait()
	
	else:
	
		p4=sp.Popen("%s/.local/bin/hybpiper retrieve_sequences  --targetfile_dna %s dna --sample_names %s.name --fasta_dir results/%s/%s_dna_seqs >>logs/%s.out 2>>logs/%s.err" %(username,baits,i,i,i,i,i),shell=True).wait()
		
		p5=sp.Popen("%s/.local/bin/hybpiper retrieve_sequences --targetfile_dna %s aa --sample_names %s.name --fasta_dir results/%s/%s_aa_seqs >>logs/%s.out 2>>logs/%s.err" %(username,baits,i,i,i,i,i),shell=True).wait()
				
	print("Retrieving paralogs for sample",i)
	
	if bait_type=="aa":
		p6=sp.Popen("%s/.local/bin/hybpiper paralog_retriever %s.name --targetfile_aa %s --fasta_dir_all results/%s/paralogs_all --fasta_dir_no_chimeras results/%s/paralogs_no_chimeras --paralog_report_filename results/%s/paralog_report --paralogs_above_threshold_report_filename results/%s/paralogs_above_threshold_report --heatmap_filename results/%s/paralog_heatmap >>logs/%s.out 2>>logs/%s.err" %(username,i,baits,i,i,i,i,i,i,i),shell=True).wait()
	
	else:
		p6=sp.Popen("%s/.local/bin/hybpiper paralog_retriever %s.name --targetfile_dna %s --fasta_dir_all results/%s/paralogs_all --fasta_dir_no_chimeras results/%s/paralogs_no_chimeras --paralog_report_filename results/%s/paralog_report --paralogs_above_threshold_report_filename results/%s/paralogs_above_threshold_report --heatmap_filename results/%s/paralog_heatmap >>logs/%s.out 2>>logs/%s.err" %(username,i,baits,i,i,i,i,i,i,i),shell=True).wait()
		
		
	#p7=sp.Popen("rm results/%s/paralog_heatmap.png" %i,shell=True).wait()
	
	p8=sp.Popen("tar -czf results/%s/%s.tar.gz %s/" %(i,i,i),shell=True).wait()
	
	p9=sp.Popen("rm -r %s/" %i,shell=True).wait()
	
	p10=sp.Popen("rm %s.name" %i,shell=True).wait()
	
	#remove seq files - not needed ### HARVEY - can change to retain supercontig & intron seqs 
	#p11=sp.Popen("rm -r results/%s/%s_intron_seqs/ results/%s/%s_supercontig_seqs/ results/%s/%s_aa_seqs/ results/%s/%s_dna_seqs/" %(i,i,i,i,i,i),shell=True).wait()
	p11=sp.Popen("rm -r results/%s/%s_aa_seqs/ results/%s/%s_dna_seqs/" %(i,i,i,i),shell=True).wait()

	
	g=open("done.txt",'a')
	g.write(i+"\n")
	g.close()
	#remove log lines with 'ETA' which fill up files
	p11=sp.Popen("sed -i '/ETA:/d' logs/%s.err" %i,shell=True).wait()
	
	#p12=sp.Popen("rsync -rut ./results $PBS_O_WORKDIR/",shell=True).wait() ##HARVEY
	#p12=sp.Popen("cp ./done.txt $PBS_O_WORKDIR/",shell=True).wait() ##HARVEY
	p12=sp.Popen("rsync -rut ./results $SLURM_SUBMIT_DIR/",shell=True).wait()
	p12=sp.Popen("cp ./done.txt $SLURM_SUBMIT_DIR/",shell=True).wait()
	
reads_dir=sys.argv[1]

baits=sys.argv[2]

bait_type=sys.argv[3]

username=sys.argv[4]

ass_options = sys.argv[5]

samplethreads=int(sys.argv[6])

cputhreads=int(sys.argv[7])
	
#get list of files, files done and list of 8 to run
filelist=os.listdir(reads_dir)

full_list=[]

p0=sp.Popen("mkdir -p results/ logs/", shell=True).wait()

for i in filelist:
	#k=i.split("_")[0] ### Old
	k=re.split(r"_H.{8}_", i) ### Harvey modification 8 May
	if k not in full_list:
		full_list.append(k)

if "done.txt" not in os.listdir():
	p0=sp.Popen("touch done.txt",shell=True).wait()

done_list=[]
g=open("done.txt",'r')
for i in g:
	done_list.append(i.rstrip("\n"))
g.close()

todo=[]
for i in full_list:
	if i not in done_list:
		todo.append(i)

print("samples to run",todo)

if __name__ == '__main__': 

	executor1 = concurrent.futures.ProcessPoolExecutor(samplethreads)
	futures1 = [executor1.submit(run_hp,i) for i in todo]
	concurrent.futures.wait(futures1)

	#debug one sample at a time:
	#for i in todo:
		#run_hp(i)

