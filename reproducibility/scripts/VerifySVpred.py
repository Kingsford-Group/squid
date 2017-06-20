#/bin/python

import sys
import pysam
sys.path.append('../HCC1954')
from BPSVclass import *

def ReadTruth(filename):
	SVs=[]
	IDs=[]
	fp=open(filename, 'r')
	for line in fp:
		if line[0]=='#':
			continue
		strs=line.strip().split("\t")
		bp1=BP_t(strs[0], int(strs[1]), int(strs[2]), (strs[8]=="-"))
		bp2=BP_t(strs[3], int(strs[4]), int(strs[5]), (strs[9]=="-"))
		sv=SV_t(bp1, bp2)
		SVs.append(sv)
		IDs.append(int(strs[6]))
	fp.close()
	return [SVs, IDs]

def BuildReference(RNAbam):
	RefTable={}
	RefLength=[]
	RefName=[]
	count=0
	fp=pysam.pysam.AlignmentFile(RNAbam)
	for e in fp.header["SQ"]:
		RefTable[e["SN"]]=count
		RefLength.append(e["LN"])
		RefName.append(e["SN"])
		count+=1
	fp.close()
	return [RefTable, RefLength, RefName]

def ReadSQUID(filename, RefName):
	SVs=[]
	fp=open(filename, 'r')
	for line in fp:
		if line[0]=='#':
			continue
		strs=line.strip().split("\t")
		bp1=BP_t(strs[0], int(strs[1]), int(strs[2]), (strs[8]=="-"))
		bp2=BP_t(strs[3], int(seg2[4]), int(seg2[5]), (strs[9]=="-"))
		sv=SV_t(bp1, bp2)
		SVs.append(sv)
	fp.close()
	return SVs

def ReadLumpy(filename):
	fp=open(filename, 'r')
	SVs=[]
	correctedcount=0
	for line in fp:
		if line[0]=='#':
			continue
		strs=line.strip().split()
		if strs[4][0]=='A' or strs[4][0]=='T' or strs[4][0]=='C' or strs[4][0]=='G':
			continue
		if strs[4]=="<DUP>":
			bp1=BP_t(strs[0], int(strs[1]), int(strs[1])+1, True)
			s=strs[7].index(";END")
			t=strs[7].index(";", s+1)
			bp2=BP_t(strs[0], int(strs[7][s+5:t])-1, int(strs[7][s+5:t]), False)
			sv=SV_t(bp1, bp2)
			correctedcount+=1
			SVs.append(sv)
		elif strs[4]=="<INV>":
			bp1=BP_t(strs[0], int(strs[1])-1, int(strs[1]), False)
			s=strs[7].index(";END")
			t=strs[7].index(";", s+1)
			bp2=BP_t(strs[0], int(strs[7][s+5:t])-1, int(strs[7][s+5:t]), False)
			sv=SV_t(bp1, bp2)
			SVs.append(sv)
			correctedcount+=1

			bp1=BP_t(strs[0], int(strs[1]), int(strs[1])+1, True)
			bp2=BP_t(strs[0], int(strs[7][s+5:t]), int(strs[7][s+5:t])+1, True)
			sv=SV_t(bp1, bp2)
			SVs.append(sv)
		elif strs[4]=="<DEL>":
			bp1=BP_t(strs[0], int(strs[1])-1, int(strs[1]), False)
			s=strs[7].index(";END")
			t=strs[7].index(";", s+1)
			bp2=BP_t(strs[0], int(strs[7][s+5:t]), int(strs[7][s+5:t])+1, True)
			sv=SV_t(bp1, bp2)
			SVs.append(sv)
		elif strs[4]=="<TRA>":
			bp1=BP_t(strs[0], int(strs[1]), int(strs[1]), True)
			s1=strs[7].index(";CHR2")
			t1=strs[7].index(";", s1+1)
			s2=strs[7].index(";END")
			t2=strs[7].index(";", s2+1)
			bp2=BP_t(strs[7][s1+6:t1], int(strs[7][s2+5:t2]), int(strs[7][s2+5:t2]), True)
			s=strs[7].index(";CT")
			if strs[7][s+4]=='3' and strs[7][s+7]=='3':
				bp1.StartPos-=1
				bp1.IsLeft=False
				bp2.StartPos-=1
				bp2.IsLeft=False
			elif strs[7][s+4]=='3' and strs[7][s+7]=='5':
				bp1.StartPos-=1
				bp1.IsLeft=False
				bp2.EndPos+=1
				bp2.IsLeft=True
			elif strs[7][s+4]=='5' and strs[7][s+7]=='3':
				bp1.EndPos+=1
				bp1.IsLeft=True
				bp2.StartPos-=1
				bp2.IsLeft=False
			elif strs[7][s+4]=='5' and strs[7][s+7]=='5':
				bp1.EndPos+=1
				bp1.IsLeft=True
				bp2.EndPos+=1
				bp2.IsLeft=True
			sv=SV_t(bp1, bp2)
			SVs.append(sv)
			correctedcount+=1
		else:
			if strs[2][-2:]=="_2":
				continue
			if strs[4][0]=='N' and strs[4][1]=='[':
				bp1=BP_t(strs[0], int(strs[1])-1, int(strs[1]), False)
				info=strs[4][2:-1].split(":")
				bp2=BP_t(info[0], int(info[1]), int(info[1])+1, True)
			elif strs[4][0]=='N' and strs[4][1]==']':
				bp1=BP_t(strs[0], int(strs[1])-1, int(strs[1]), False)
				info=strs[4][2:-1].split(":")
				bp2=BP_t(info[0], int(info[1])-1, int(info[1]), False)
			elif strs[4][-1]=='N' and strs[4][-2]=='[':
				bp1=BP_t(strs[0], int(strs[1]), int(strs[1])+1, True)
				info=strs[4][1:-2].split(":")
				bp2=BP_t(info[0], int(info[1]), int(info[1])+1, True)
			else:
				bp1=BP_t(strs[0], int(strs[1]), int(strs[1])+1, True)
				info=strs[4][1:-2].split(":")
				bp2=BP_t(info[0], int(info[1])-1, int(info[1]), False)
			sv=SV_t(bp1, bp2)
			correctedcount+=1
			SVs.append(sv)
	fp.close()
	return [SVs, correctedcount]

def ReadBreakDancer(filename):
	fp=open(filename, 'r')
	SVs=[]
	correctedcount=0
	for line in fp:
		if line[0]=='#':
			continue
		strs=line.strip().split("\t")
		if strs[6]=="DEL":
			bp1=BP_t(strs[0], int(strs[1])-1, int(strs[1]), False)
			bp2=BP_t(strs[3], int(strs[4]), int(strs[4])+1, True)
			sv=SV_t(bp1, bp2)
			SVs.append(sv)
		elif strs[6]=="INV":
			bp1=BP_t(strs[0], int(strs[1])-1, int(strs[1]), False)
			bp2=BP_t(strs[3], int(strs[4])-1, int(strs[4]), False)
			sv=SV_t(bp1, bp2)
			correctedcount+=1
			SVs.append(sv)
			bp1=BP_t(strs[0], int(strs[1]), int(strs[1])+1, True)
			bp2=BP_t(strs[3], int(strs[4]), int(strs[4])+1, True)
			sv=SV_t(bp1, bp2)
			SVs.append(sv)
		elif strs[6]=="ITX":
			bp1=BP_t(strs[0], int(strs[1])-1, int(strs[1]), False)
			bp2=BP_t(strs[3], int(strs[4])-1, int(strs[4]), False)
			sv=SV_t(bp1, bp2)
			correctedcount+=1
			SVs.append(sv)
			bp1=BP_t(strs[0], int(strs[1]), int(strs[1])+1, True)
			bp2=BP_t(strs[3], int(strs[4]), int(strs[4])+1, True)
			sv=SV_t(bp1, bp2)
			SVs.append(sv)
			bp1=BP_t(strs[0], int(strs[1])-1, int(strs[1]), False)
			bp2=BP_t(strs[3], int(strs[4]), int(strs[4])+1, True)
			sv=SV_t(bp1, bp2)
			SVs.append(sv)
			bp1=BP_t(strs[0], int(strs[1]), int(strs[1])+1, True)
			bp2=BP_t(strs[3], int(strs[4])-1, int(strs[4]), False)
			sv=SV_t(bp1, bp2)
			SVs.append(sv)
		elif strs[6]=="CTX":
			bp1=BP_t(strs[0], int(strs[1])-1, int(strs[1]), False)
			bp2=BP_t(strs[3], int(strs[4])-1, int(strs[4]), False)
			sv=SV_t(bp1, bp2)
			correctedcount+=1
			SVs.append(sv)
			bp1=BP_t(strs[0], int(strs[1]), int(strs[1])+1, True)
			bp2=BP_t(strs[3], int(strs[4]), int(strs[4])+1, True)
			sv=SV_t(bp1, bp2)
			SVs.append(sv)
			bp1=BP_t(strs[0], int(strs[1])-1, int(strs[1]), False)
			bp2=BP_t(strs[3], int(strs[4]), int(strs[4])+1, True)
			sv=SV_t(bp1, bp2)
			SVs.append(sv)
			bp1=BP_t(strs[0], int(strs[1]), int(strs[1])+1, True)
			bp2=BP_t(strs[3], int(strs[4])-1, int(strs[4]), False)
			sv=SV_t(bp1, bp2)
			SVs.append(sv)
	fp.close()
	return [SVs,correctedcount]

def ReadHYDRA(filename, FilterLowQual):
	fp=open(filename, 'r')
	SVs=[]
	correctedcount=0
	for line in fp:
		if line[0]=='#':
			continue
		strs=line.strip().split("\t")
		if int(strs[7])<FilterLowQual:
			continue;
		if strs[10]!="local_deletion":
			correctedcount+=1
		bp1=BP_t(strs[0], int(strs[1]), int(strs[2]), True)
		if strs[8]=="+":
			bp1.IsLeft=False
		bp2=BP_t(strs[3], int(strs[4]), int(strs[5]), True)
		if strs[9]=="+":
			bp2.IsLeft=False
		sv=SV_t(bp1, bp2)
		SVs.append(sv)
	fp.close()
	return [SVs, correctedcount]

def ReadBedpe(filename):
	fp=open(filename, 'r')
	SVs=[]
	for line in fp:
		if line[0]=='#':
			continue
		strs=line.strip().split("\t")
		bp1=BP_t(strs[0], int(strs[1]), int(strs[2]), True)
		if strs[8]=="+":
			bp1.IsLeft=False
		bp2=BP_t(strs[3], int(strs[4]), int(strs[5]), True)
		if strs[9]=="+":
			bp2.IsLeft=False
		sv=SV_t(bp1, bp2)
		SVs.append(sv)
	fp.close()
	return SVs

def CompareSV(trueSV, predSV):
	global thresh
	Result=[-1]*len(predSV)
	for i in range(len(predSV)):
		psv=predSV[i]
		for j in range(len(trueSV)):
			sv=trueSV[j]
			if psv.BP1.Chr==sv.BP1.Chr and psv.BP2.Chr==sv.BP2.Chr and psv.BP1.IsLeft==sv.BP1.IsLeft and psv.BP2.IsLeft==sv.BP2.IsLeft:
				if psv.BP1.IsLeft and psv.BP2.IsLeft and abs(psv.BP1.StartPos-sv.BP1.StartPos)<thresh and abs(psv.BP2.StartPos-sv.BP2.StartPos)<thresh:
					Result[i]=j
				elif psv.BP1.IsLeft and (not psv.BP2.IsLeft) and abs(psv.BP1.StartPos-sv.BP1.StartPos)<thresh and abs(psv.BP2.EndPos-sv.BP2.EndPos)<thresh:
					Result[i]=j
				elif (not psv.BP1.IsLeft) and psv.BP2.IsLeft and abs(psv.BP1.EndPos-sv.BP1.EndPos)<thresh and abs(psv.BP2.StartPos-sv.BP2.StartPos)<thresh:
					Result[i]=j
				elif (not psv.BP1.IsLeft) and (not psv.BP2.IsLeft) and abs(psv.BP1.EndPos-sv.BP1.EndPos)<thresh and abs(psv.BP2.EndPos-sv.BP2.EndPos)<thresh:
					Result[i]=j
	return Result

def OutputResult(accuracy, sensitivity, predSV, Result, IDs, outputfile):
	fp=open(outputfile, 'w')
	fp.write("#accuracy={:.6f}\tsensitivity={:.6f}\n".format(accuracy, sensitivity))
	fp.write("#chrom1\tstart1\tend1\tstrand1\tchrom2\tstart2\tend2\tstrand2\tis_hit\tsvID\n")
	for i in range(len(predSV)):
		sv=predSV[i]
		isleft1="+"
		if sv.BP1.IsLeft:
			isleft1="-"
		isleft2="+"
		if sv.BP2.IsLeft:
			isleft2="-"
		ishit=0
		svID=-1
		if Result[i]!=-1:
			ishit=1
			svID=IDs[Result[i]]
		fp.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sv.BP1.Chr, sv.BP1.StartPos, sv.BP1.EndPos, isleft1, sv.BP2.Chr, sv.BP2.StartPos, sv.BP2.EndPos, isleft2, ishit, svID))
	fp.close()

if __name__=="__main__":
	if len(sys.argv)<3:
		print("python VerifySQUID.py <option of 1-5> <TrueSV> <PredSV> <Outputfile> [<RNAbam> | <qual_threshold>]")
		print("\t1: SQUIDdiscordant, RNAbam needed")
		print("\t2: LUMPY")
		print("\t3: BreakDancer")
		print("\t4: HYDRA, qual_threshold needed")
		print("\t5: de novo + mummer/gmap")
	else:
		thresh=10
		TrueSVFile=sys.argv[2]
		PredSVFile=sys.argv[3]
		OutputFile=sys.argv[4]

		[trueSV, IDs]=ReadTruth(TrueSVFile)
		if sys.argv[1]=="1":
			RNAbam=sys.argv[5]
			[RefTable, RefLength, RefName]=BuildReference(RNAbam)
			predSV=ReadSQUID(PredSVFile, RefName)
			correctedcount=len(predSV)
		elif sys.argv[1]=="2":
			[predSV, correctedcount]=ReadLumpy(PredSVFile)
		elif sys.argv[1]=='3':
			[predSV, correctedcount]=ReadBreakDancer(PredSVFile)
		elif sys.argv[1]=='4':
			FilterLowQual=int(sys.argv[5])
			[predSV, correctedcount]=ReadHYDRA(PredSVFile, FilterLowQual)
		elif sys.argv[1]=='5':
			predSV=ReadBedpe(PredSVFile)
			correctedcount=len(predSV)
		Result=CompareSV(trueSV, predSV)

		accuracy=sum([(x!=-1) for x in Result])
		accuracy=1.0*accuracy/correctedcount
		sensitivity=len(set(Result))-1
		sensitivity=1.0*sensitivity/len(trueSV)
		print("Accuracy = {}\tsensitivity = {}".format(accuracy, sensitivity))
		OutputResult(accuracy, sensitivity, predSV, Result, IDs, OutputFile)
