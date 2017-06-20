#!/bin/python

import sys
import pysam

def IsReadConcordant(read):
	thresh=3
	flag=True
	if read.is_unmapped or read.mate_is_unmapped:
		flag=False
	elif read.reference_id!=read.next_reference_id:
		flag=False
	elif read.is_reverse==read.mate_is_reverse:
		flag=False
	elif read.is_reverse and read.reference_start<read.next_reference_start-thresh:
		flag=False
	elif not read.is_reverse and read.reference_start>read.next_reference_start+thresh:
		flag=False
	return flag

def TrueDiscbyBWA(speedseqbam):
	DisNames=[]
	fp=pysam.AlignmentFile(speedseqbam)
	for read in fp:
		if read.is_secondary or read.is_supplementary:
			continue
		if not IsReadConcordant(read) and read.mapping_quality>0 and not read.has_tag("XA"):
			DisNames.append(read.query_name)
	fp.close()
	DisNames.sort()

	FilteredDisName=[]
	i=1
	while i<len(DisNames):
		if DisNames[i-1]==DisNames[i]:
			FilteredDisName.append(DisNames[i])
		i+=1
	print("There are {} valid discordant read pairs".format(len(FilteredDisName)))
	return set(FilteredDisName)

def FilterChimBam(FilteredDisName, chimbam, outbam):
	fpin=pysam.AlignmentFile(chimbam)
	fpout=pysam.AlignmentFile(outbam, "wb", template=fpin)
	for read in fpin:
		if read.query_name in FilteredDisName:
			fpout.write(read)
	fpin.close()
	fpout.close()

if __name__=="__main__":
	if len(sys.argv)<4:
		print("python FilterChimOut.py <SpeedSeqBam> <ChimericBam> <OutBam>")
	else:
		SpeedSeqBam=sys.argv[1]
		ChimericBam=sys.argv[2]
		OutBam=sys.argv[3]

		FilteredDisName=TrueDiscbyBWA(SpeedSeqBam)
		FilterChimBam(FilteredDisName, ChimericBam, OutBam)