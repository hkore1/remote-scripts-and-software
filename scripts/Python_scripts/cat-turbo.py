#!/usr/bin/env python3
import sys
import subprocess as sp
import multiprocessing
import concurrent.futures
import re
import glob
import os

#cat-turbo.py infolder/ catfiles/ 24

def tokenize(filename):
	digits = re.compile(r'(\d+)')
	return tuple(int(token) if match else token for token, match in ((fragment, digits.search(fragment)) for fragment in digits.split(filename)))


def catturbo(i):
	
	for x in filelist:
		if x.split("/")[-1].split("_")[0]==i:
			ext=x.split("/")[-1].split(".")[-1]
			break
	
	if ext=="gz":
	
		print('cat',i)
		name1=i+"_*R1*.gz"
		name2=i+"_R1.fastq.gz"
		
		p0=sp.Popen("cat  %s/%s >%s/%s" %(infolder,name1,outfolder,name2),shell=True).wait()

		name3=i+"_*R2*.gz"
		name4=i+"_R2.fastq.gz"
		
		p2=sp.Popen("cat  %s/%s >%s/%s" %(infolder,name3,outfolder,name4),shell=True).wait()

	else:
	
		print('cat',i)
		name1=i+"_*R1*"
		name2=i+"_R1.fastq"
		
		p0=sp.Popen("cat  %s/%s >%s/%s" %(infolder,name1,outfolder,name2),shell=True).wait()

		name3=i+"_*R2*"
		name4=i+"_R2.fastq"
		
		p2=sp.Popen("cat  %s/%s >%s/%s" %(infolder,name3,outfolder,name4),shell=True).wait()

infolder=sys.argv[1]

filelist=os.listdir(infolder)

namelist=[]

for i in filelist:
	k=i.split("_")[0]
	if k not in namelist:
		namelist.append(k)

outfolder = sys.argv[2] #outfolder
threads=int(sys.argv[3])

if outfolder == infolder:
	print("You cannot use the same infolder and outfolder, please change outfolder. Exiting.")
	quit()

p1=sp.Popen("mkdir -p %s" %outfolder,shell=True).wait()
			

if __name__ == '__main__':

	print("cat..",threads,"threads")
	executor1 = concurrent.futures.ProcessPoolExecutor(threads)
	futures1 = [executor1.submit(catturbo,i) for i in namelist]
	concurrent.futures.wait(futures1)
	
	#for x in namelist:
		#catturbo(x)
	
	print("cat-turbo completed")