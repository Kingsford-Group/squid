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
	if len(sys.argv)!=5:
		print("python3 Read2Remap.py <Input_BAM> <Input_BED> <FirstRead.fastq> <Out_Prefix>")
	else:
		Input_BAM=sys.argv[1]
		Input_BED=sys.argv[2]
		In_Prefix=sys.argv[3][:(sys.argv[3].rfind("_"))]
		In_Suffix=sys.argv[3][(sys.argv[3].rfind("_")+3):]
		Out_Prefix=sys.argv[4]

		NameList1=UnmapList(Input_BAM)
		NameList2=ReadBEDList(Input_BAM, Input_BED)
		NameList=NameList1|NameList2
		Reads2Remap(In_Prefix, In_Suffix, Out_Prefix, NameList)