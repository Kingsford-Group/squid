#!/bin/python

import sys
import subprocess

def UnmapList(bamfile):
	commandbam="samtools view -f4 "+bamfile+" | cut -f1"
	p=subprocess.Popen(commandbam, shell=True, stdout=subprocess.PIPE)
	namebam=p.communicate()[0]
	namebam=namebam.decode("utf-8").strip()
	namebam=namebam.split("\n")
	namebam=set(namebam)

	print(len(namebam))
	return namebam

def UnmapListfromStar(unmapprefix):
	fp=open(unmapprefix+"1", 'r')
	namestar=[]
	linecount=0
	for line in fp:
		linecount+=1
		if linecount%4==1:
			strs=line.strip().split()
			namestar.append(strs[0][1:])
	fp.close()
	fp=open(unmapprefix+"2", 'r')
	linecount=0
	for line in fp:
		linecount+=1
		if linecount%4==1:
			strs=line.strip().split()
			namestar.append(strs[0][1:])
	fp.close()
	namestar=set(namestar)
	return namestar

def ReadBEDList(bamfile, bedfile):
	command="samtools view -L "+bedfile+" "+bamfile+" | cut -f1"
	p=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
	namebed=p.communicate()[0]
	namebed=namebed.decode("utf-8").strip()
	namebed=namebed.split("\n")
	namebed=set(namebed)

	return namebed

def Reads2Remap(in_prefix, in_suffix, out_prefix, NameList):
	fpout=open(out_prefix+"_1.fastq", 'w')
	command=""
	if in_suffix[-2:]=="gz":
		command="gunzip -c "+in_prefix+"_1."+in_suffix
	else:
		command="cat "+in_prefix+"_1."+in_suffix
	p=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
	reads=p.communicate()[0]
	reads=reads.decode("utf-8").strip()
	reads=reads.split("\n")
	i=0
	while i<len(reads):
		if reads[i].split()[0][1:] in NameList: 
			for j in range(4):
				fpout.write(reads[i+j]+"\n")
		i+=4
	fpout.close()

	fpout=open(out_prefix+"_2.fastq", 'w')
	if in_suffix[-2:]=="gz":
		command="gunzip -c "+in_prefix+"_2."+in_suffix
	else:
		command="cat "+in_prefix+"_2."+in_suffix
	p=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
	reads=p.communicate()[0]
	reads=reads.decode("utf-8").strip()
	reads=reads.split("\n")
	i=0
	while i<len(reads):
		if reads[i].split()[0][1:] in NameList: 
			for j in range(4):
				fpout.write(reads[i+j]+"\n")
		i+=4
	fpout.close()

if __name__=="__main__":
	if len(sys.argv)<5:
		print("python3 Read2Remap.py 1 <Input_BAM> <Input_BED> <FirstRead.fastq> <Out_Prefix>")
		print("python3 Read2Remap.py 2 <StarUnmap> <Input_BAM> <Input_BED> <FirstRead.fastq> <Out_Prefix>")
	else:
		if sys.argv[1]=='1':
			Input_BAM=sys.argv[2]
			Input_BED=sys.argv[3]
			In_Prefix=sys.argv[4][:(sys.argv[4].rfind("_"))]
			In_Suffix=sys.argv[4][(sys.argv[4].rfind("_")+3):]
			Out_Prefix=sys.argv[5]
		elif sys.argv[1]=='2':
			Input_StarUnmap=sys.argv[2]
			Input_BAM=sys.argv[3]
			Input_BED=sys.argv[4]
			In_Prefix=sys.argv[5][:(sys.argv[5].rfind("_"))]
			In_Suffix=sys.argv[5][(sys.argv[5].rfind("_")+3):]
			Out_Prefix=sys.argv[6]

		if sys.argv[1]=='1':
			NameList1=UnmapList(Input_BAM)
		elif sys.argv[1]=='2':
			NameList1=UnmapListfromStar(Input_StarUnmap)
		NameList2=ReadBEDList(Input_BAM, Input_BED)
		NameList=NameList1|NameList2
		Reads2Remap(In_Prefix, In_Suffix, Out_Prefix, NameList)