#!/usr/bin/env python3

######################
# Author: Theo Allnutt
# Modified: Harvey Orel (21 Nov 2022 - added join 'paralog_report.tsv')
#						(22 Nov 2022 - added join supercontigs, for intronerated assemblies)
#						(9 May 2023  - added join introns, for intronerated assemblies)
#		Modification spits out an extra txt file 'compiled_paralogs.txt' that can be used to plot 
#		a heatmap for viewing distribution of paralogs across samples.
######################

import sys
import re
import glob
import subprocess as sp
import os
import concurrent.futures

digits = re.compile(r'(\d+)')
def tokenize(filename):
    return tuple(int(token) if match else token
                 for token, match in
                 ((fragment, digits.search(fragment))
                  for fragment in digits.split(filename)))

#geneName.paralogs.fasta

def join_paralogs(i,s):

	files1=os.listdir(resdir+i+"/paralogs_all/")
	
	print("Joining stats",i)
	p6=sp.Popen("head -2 %s/%s/%s_stats.tsv |tail -1 >> %s/stats_joined.txt" %(resdir,i,i,resdir),shell=True).wait

	for j in files1:
		j1=j.replace("_paralogs_all",".paralogs")
		#print("Joining %s/paralogs_all/" %i,j)
		p2=sp.Popen("cat %s/%s/paralogs_all/%s >> %s/all_paralogs_joined/%s" %(resdir,i,j,resdir,j1),shell=True).wait()
		
	files2=os.listdir(resdir+i+"/paralogs_no_chimeras/")
	
	for j in files2:
		j1=j.replace("_paralogs_no_chimeras",".paralogs")
		#print("Joining %s/paralogs_no_chimeras/" %i,j)
		p3=sp.Popen("cat %s/%s/paralogs_no_chimeras/%s >> %s/paralogs_no_chimeras_joined/%s" %(resdir,i,j,resdir,j1),shell=True).wait()

	if s == "y": ### Harvey additions 22 Nov 2022
		files3=os.listdir(resdir+i+"/"+i+"_supercontig_seqs/")
		for j in files3:
			p4=sp.Popen("cat %s/%s/%s_supercontig_seqs/%s >> %s/supercontig_seqs_joined/%s" %(resdir,i,i,j,resdir,j),shell=True).wait()

		### Harvey additions 9 May 2023
		files4=os.listdir(resdir+i+"/"+i+"_intron_seqs/")
		for j in files4:
			p5=sp.Popen("cat %s/%s/%s_intron_seqs/%s >> %s/intron_seqs_joined/%s" %(resdir,i,i,j,resdir,j),shell=True).wait()
	else:
		pass

	
def get_lengths():
	
	for i in dirs:
	
		g2=open(resdir+"/"+i+"/"+i+"_lengths.tsv",'r')
		#each file genes are in different order
		
		#genes
		genes=[]
		
		k=g2.readline().rstrip("\n").split("\t")[1:]
		
		for x in k:
		
			genes.append(x)
		
			if x not in allgenes:
				allgenes.append(x)
				
		#means only for first file, should all be the same
		if len(lens.keys())<2:
			
			k=g2.readline().rstrip("\n").split("\t")[1:]
			c=-1
			for x in k:
				if x not in means.keys():
					c=c+1
					means[genes[c]]=x
					
		else:
			g2.readline() #skip means other than first file
				
		#lengths
		lens[i]={}
		k=g2.readline().rstrip("\n").split("\t")[1:]
		c=-1
		if i not in species:
			species.append(i)
			
		for x in k:
			c=c+1
			lens[i][genes[c]]=x
				
		g2.close()

def count_paralogs():
	
	for i in dirs:
	
		try:
			g4=open(resdir+"/"+i+"/"+"paralog_report.tsv",'r')

			#genes
			genes=[]
			
			k=g4.readline().rstrip("\n").split("\t")[1:]
			
			for x in k:
			
				genes.append(x)
			
				if x not in allgenes:
					allgenes.append(x)
					
			#lengths
			lens[i]={}
			k=g4.readline().rstrip("\n").split("\t")[1:]
			c=-1
			if i not in species:
				species.append(i)
				
			for x in k:
				c=c+1
				lens[i][genes[c]]=x
					
			g4.close()

		except FileNotFoundError:
			pass
		#each file genes are in different order
		
		
	#return genes,lens


################global

resdir=sys.argv[1]
threads=int(sys.argv[2])
dirs=next(os.walk(resdir))[1]
join_intronerate=sys.argv[3]

if "paralogs_no_chimeras_joined" in dirs:
	dirs.remove("paralogs_no_chimeras_joined")
if "all_paralogs_joined" in dirs:
	dirs.remove("all_paralogs_joined")

print("N.B. deleting any existing directories of joined files")

p0=sp.Popen("rm -rf %s/all_paralogs_joined/" %resdir,shell=True).wait()
p1=sp.Popen("rm -rf %s/paralogs_no_chimeras_joined/" %resdir,shell=True).wait()

p0=sp.Popen("mkdir -p %s/all_paralogs_joined/" %resdir,shell=True).wait()
p1=sp.Popen("mkdir -p %s/paralogs_no_chimeras_joined/" %resdir,shell=True).wait()

if join_intronerate == "y": ### Harvey additions
	print("\nJoining intronerate files (supercontigs and introns)\n")
	if "supercontig_seqs_joined" in dirs:
		dirs.remove("supercontig_seqs_joined")

	if "intron_seqs_joined" in dirs:
		dirs.remove("intron_seqs_joined")

	p7=sp.Popen("rm -rf %s/supercontig_seqs_joined/" %resdir,shell=True).wait()
	p8=sp.Popen("rm -rf %s/intron_seqs_joined/" %resdir,shell=True).wait()

	p7=sp.Popen("mkdir -p %s/supercontig_seqs_joined/" %resdir,shell=True).wait()
	p8=sp.Popen("mkdir -p %s/intron_seqs_joined/" %resdir,shell=True).wait()
else:
	print("\nNot joining intronerate files (supercontigs and introns)\n")
	pass

g1=open(resdir+"/stats_joined.txt",'w')
g1.write("Name\tNumReads\tReadsMapped\tPctOnTarget\tGenesMapped\tGenesWithContigs\tGenesWithSeqs\tGenesAt25pct\tGenesAt50pct\tGenesAt75pct\tGenesAt150pct\tParalogWarningsLong\tParalogWarningsDepth\tGenesWithoutStitchedContigs\tGenesWithStitchedContigs\tGenesWithStitchedContigsSkipped\tGenesWithChimeraWarning\n")
g1.close()

lens={}
species=[]
genes=[]
means={}
allgenes=[]

if __name__ == '__main__':
	
	executor1 = concurrent.futures.ProcessPoolExecutor(threads)
	futures1 = [executor1.submit(join_paralogs,i,join_intronerate) for i in dirs] ### Harvey additions
	concurrent.futures.wait(futures1)
	
	#get lengths
	print("Getting lengths")
	get_lengths()
	
	
	#print(genes)
	#print(lens)
	
	allgenes.sort(key=tokenize)

	g3=open(resdir+"/"+"all_lengths.txt",'w')

	g3.write("Species\t"+"\t".join(str(p) for p in allgenes)+"\n")

	g3.write("MeanLength\t")
	#print(lens)
	#print('genes',genes)
	for x in allgenes:
		g3.write(means[x]+"\t")
	g3.write("\n")

	species.sort(key=tokenize)

	for x in species:
		g3.write(x)
		for y in allgenes:
		
			if y in lens[x].keys():
				g3.write("\t"+lens[x][y])
				
			else:
				g3.write("\t0")
			
		g3.write("\n")
		
	g3.close()

	#get lengths
	print("Counting paralogs")
	count_paralogs()
	
	
	#print(genes)
	#print(lens)
	
	allgenes.sort(key=tokenize)

	g5=open(resdir+"/"+"compiled_paralogs.txt",'w')

	g5.write("Species\t"+"\t".join(str(p) for p in allgenes)+"\n")

	#print(lens)
	#print('genes',genes)

	species.sort(key=tokenize)

	for x in species:
		g5.write(x)
		for y in allgenes:
		
			if y in lens[x].keys():
				g5.write("\t"+lens[x][y])
				
			else:
				g5.write("\t0")
			
		g5.write("\n")
		
	g5.close()

	print("Done")