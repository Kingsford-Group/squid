#!/bin/python

import sys
import pysam

class BP_t(object):
	def __init__(self, Chr_, Position_, IsLeft_):
		self.Chr=Chr_
		self.Position=Position_
		self.IsLeft=IsLeft_
	def __eq__(self, other):
		if isinstance(other, BP_t):
			return (self.Chr==other.Chr and self.Position==other.Position and self.IsLeft==other.IsLeft)
		return NotImplemented
	def __ne__(self, other):
		result=self.__eq__(other)
		if result is NotImplemented:
			return result
		return not result
	def __lt__(self, other):
		if isinstance(other, BP_t):
			if self.Chr!=other.Chr:
				return self.Chr<other.Chr
			elif self.Position!=other.Position:
				return self.Position<other.Position
			else:
				return self.IsLeft<other.IsLeft
		return NotImplemented
	def __gt__(self, other):
		if isinstance(other, BP_t):
			if self.Chr!=other.Chr:
				return self.Chr>other.Chr
			elif self.Position!=other.Position:
				return self.Position>other.Position
			else:
				return self.IsLeft>other.IsLeft
		return NotImplemented
	def __le__(self, other):
		result=self.__gt__(other)
		if result is NotImplemented:
			return result
		return not result
	def __ge__(self, other):
		result=self.__lt__(other)
		if result is NotImplemented:
			return result
		return not result
	def Print(self):
		print("{}\t{}\t{}".format(self.Chr, self.Position, self.IsLeft))
	def ToString(self):
		return ("{}\t{}\t{}".format(self.Chr, self.Position, self.IsLeft))

class SV_t(object):
	def __init__(self, bp1, bp2):
		if bp1<bp2:
			self.BP1=bp1
			self.BP2=bp2
		else:
			self.BP1=bp2
			self.BP2=bp1
	def __eq__(self, other):
		if isinstance(other, SV_t):
			return (self.BP1==other.BP1 and self.BP2==other.BP2)
		return NotImplemented
	def __ne__(self, other):
		result=self.__eq__(other)
		if result is NotImplemented:
			return result
		return not result
	def __lt__(self, other):
		if isinstance(other, SV_t):
			if self.BP1!=other.BP1:
				return self.BP1<other.BP1
			else:
				return self.BP2<other.BP2
		return NotImplemented
	def __gt__(self, other):
		if isinstance(other, SV_t):
			if self.BP1!=other.BP1:
				return self.BP1>other.BP1
			else:
				return self.BP2>other.BP2
		return NotImplemented
	def __le__(self, other):
		result=self.__gt__(other)
		if result is NotImplemented:
			return result
		return not result
	def __ge__(self, other):
		result=self.__lt__(other)
		if result is NotImplemented:
			return result
		return not result
	def Print(self):
		strbp1=self.BP1.ToString()
		strbp2=self.BP2.ToString()
		print(strbp1+"\t|\t"+strbp2)
	def ToString(self):
		strbp1=self.BP1.ToString()
		strbp2=self.BP2.ToString()
		return strbp1+"\t|\t"+strbp2

class Exon_t(object):
	def __init__(self, _GeneID, _Chr, _StartPos, _EndPos):
		self.GeneID=_GeneID
		self.Chr=_Chr
		self.StartPos=_StartPos
		self.EndPos=_EndPos
	def __eq__(self, other):
		if isinstance(other, Exon_t):
			return (self.GeneID==other.GeneID and self.Chr==other.Chr and self.StartPos==other.StartPos, self.EndPos==other.EndPos)
		return NotImplemented
	def __ne__(self, other):
		result=self.__eq__(other)
		if result is NotImplemented:
			return result
		return not result
	def __lt__(self, other):
		if isinstance(other, Exon_t):
			if self.Chr!=other.Chr:
				return self.Chr<other.Chr
			elif self.StartPos!=other.StartPos:
				return self.StartPos<other.StartPos
			else:
				return self.EndPos<other.EndPos
		return NotImplemented
	def __gt__(self, other):
		if isinstance(other, Exon_t):
			if self.Chr!=other.Chr:
				return self.Chr>other.Chr
			elif self.StartPos!=other.StartPos:
				return self.StartPos>other.StartPos
			else:
				return self.EndPos>other.EndPos
		return NotImplemented
	def __le__(self, other):
		result=self.__gt__(other)
		if result is NotImplemented:
			return result
		return not result
	def __ge__(self, other):
		result=self.__lt__(other)
		if result is NotImplemented:
			return result
		return not result

class Gene_t(object):
	def __init__(self, _ID, _Name, _Chr, _StartPos, _EndPos, _Strand):
		self.ID=_ID
		self.Name=_Name
		self.Chr=_Chr
		self.StartPos=_StartPos
		self.EndPos=_EndPos
		self.Exons=[]
		self.Strand=_Strand
	def __eq__(self, other):
		if isinstance(other, Gene_t):
			return (self.Chr==other.Chr and self.StartPos==other.StartPos and self.EndPos==other.EndPos)
		return NotImplemented
	def __ne__(self, other):
		result=self.__eq__(other)
		if result is NotImplemented:
			return result
		return not result
	def __lt__(self, other):
		if isinstance(other, Gene_t):
			if self.Chr!=other.Chr:
				return self.Chr<other.Chr
			elif self.StartPos!=other.StartPos:
				return self.StartPos<other.StartPos
			else:
				return self.EndPos<other.EndPos
		return NotImplemented
	def __gt__(self, other):
		if isinstance(other, Gene_t):
			if self.Chr!=other.Chr:
				return self.Chr>other.Chr
			elif self.StartPos!=other.StartPos:
				return self.StartPos>other.StartPos
			else:
				return self.EndPos>other.EndPos
		return NotImplemented
	def __le__(self, other):
		result=self.__gt__(other)
		if result is NotImplemented:
			return result
		return not result
	def __ge__(self, other):
		result=self.__lt__(other)
		if result is NotImplemented:
			return result
		return not result
	def UniqSortExons(self):
		self.Exons=list(set(self.Exons))
		self.Exons.sort()


def ReadGTF(filename):
	Genes=[]
	UnassignedExons=[]
	fp=open(filename, 'r')
	for line in fp:
		if line[0]=='#':
			continue
		strs=line.strip().split("\t")
		if strs[2]=="gene":
			s=strs[8].index("gene_id")
			t=strs[8].index(";", s+1)
			ID=strs[8][(s+9):(t-1)]
			s=strs[8].index("gene_name")
			t=strs[8].index(";", s+1)
			Name=strs[8][(s+11):(t-1)]
			g=Gene_t(ID, Name, strs[0], int(strs[3])-1, int(strs[4]), strs[6])
			Genes.append(g)
		elif strs[2]=="exon":
			s=strs[8].index("gene_id")
			t=strs[8].index(";", s+1)
			ID=strs[8][(s+9):(t-1)]
			exon=Exon_t(ID, strs[0], int(strs[3])-1, int(strs[4]))
			assigned=False
			for i in range(len(Genes)-2, len(Genes)):
				if exon.GeneID==Genes[i].ID:
					Genes[i].Exons.append(exon)
					assigned=True
					break
			if not assigned:
				UnassignedExons.append(exon)
	fp.close()
	for exon in UnassignedExons:
		for i in range(len(Genes)):
			if exon.GeneID==Genes[i].ID:
				Genes[i].Exons.append(exon)
				break
	Genes.sort()
	for i in range(len(Genes)):
		Genes[i].UniqSortExons()
	return Genes

def ReadTrueSV(filename):
	fp=open(filename, 'r')
	SVs=[]
	for line in fp:
		if line[0]=='#':
			continue
		strs=line.strip().split()
		if strs[8]=="+":
			bp1=BP_t(strs[0], int(strs[2]), False)
		else:
			bp1=BP_t(strs[0], int(strs[1]), True)
		if strs[9]=="+":
			bp2=BP_t(strs[3], int(strs[5]), False)
		else:
			bp2=BP_t(strs[3], int(strs[4]), True)
		sv=SV_t(bp1, bp2)
		SVs.append(sv)
	fp.close()
	return SVs

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
		strs=line.strip().split()
		seg1=strs[1].split(",")
		seg2=strs[4].split(",")
		if (RefName[int(seg1[0])][0]>='0' and RefName[int(seg1[0])][0]<='9') or RefName[int(seg1[0])][0]=='X' or RefName[int(seg1[0])][0]=='Y':
			if (RefName[int(seg2[0])][0]>='0' and RefName[int(seg2[0])][0]<='9') or RefName[int(seg2[0])][0]=='X' or RefName[int(seg2[0])][0]=='Y':
				if strs[2]=='H':
					bp1=BP_t(RefName[int(seg1[0])], int(seg1[1]), True)
				else:
					bp1=BP_t(RefName[int(seg1[0])], int(seg1[1])+int(seg1[2]), False)
				if strs[5]=='H':
					bp2=BP_t(RefName[int(seg2[0])], int(seg2[1]), True)
				else:
					bp2=BP_t(RefName[int(seg2[0])], int(seg2[1])+int(seg2[2]), False)
				sv=SV_t(bp1, bp2)
				SVs.append(sv)
	fp.close()
	return SVs

def ReadFusionCatcher(filename):
	SVs=[]
	fp=open(filename, 'r')
	linecount=0
	for line in fp:
		linecount+=1
		if linecount==1:
			continue
		strs=line.strip().split("\t")
		if ":" in strs[7]:
			bp1=BP_t(strs[7].split(":")[0], int(strs[7].split(":")[1]), (strs[7].split(":")[2]=='-'))
			bp2=BP_t(strs[8].split(":")[0], int(strs[8].split(":")[1]), (strs[8].split(":")[2]=='+'))
		else:
			bp1=BP_t(strs[8].split(":")[0], int(strs[8].split(":")[1]), (strs[8].split(":")[2]=='-'))
			bp2=BP_t(strs[9].split(":")[0], int(strs[9].split(":")[1]), (strs[9].split(":")[2]=='+'))
		sv=SV_t(bp1, bp2)
		SVs.append(sv)
	fp.close()
	return SVs

def ReadJAFFA(filename):
	SVs=[]
	fp=open(filename, 'r')
	linecount=0
	for line in fp:
		linecount+=1
		if linecount==1:
			continue
		strs=line.replace("\"", "").strip().split(",")
		bp1=BP_t(strs[2], int(strs[3]), (strs[4]=='-'))
		bp2=BP_t(strs[5], int(strs[6]), (strs[7]=='+'))
		sv=SV_t(bp1, bp2)
		SVs.append(sv)
	fp.close()
	return SVs

def ReadDefuse(filename):
	SVs=[]
	fp=open(filename, 'r')
	linecount=0
	for line in fp:
		linecount+=1
		if linecount==1:
			continue
		strs=line.strip().split("\t")
		bp1=BP_t(strs[24], int(strs[37]), (strs[43]=='-'))
		bp2=BP_t(strs[25], int(strs[38]), (strs[44]=='-'))
		sv=SV_t(bp1, bp2)
		SVs.append(sv)
	fp.close()
	return SVs

def ReadChimerascan(filename):
	SVs=[]
	fp=open(filename, 'r')
	linecount=0
	for line in fp:
		linecount+=1
		if linecount==1:
			continue
		strs=line.strip().split("\t")
		if strs[8]=="-":
			bp1=BP_t(strs[0], int(strs[1]), True)
		else:
			bp1=BP_t(strs[0], int(strs[2]), False)
		if strs[9]=='+':
			bp2=BP_t(strs[3], int(strs[4]), True)
		else:
			bp2=BP_t(strs[3], int(strs[5]), False)
		sv=SV_t(bp1, bp2)
		SVs.append(sv)
	fp.close()
	return SVs

def ReadIntegrate(filename):
	SVs=[]
	fp=open(filename, 'r')
	linecount=0
	for line in fp:
		linecount+=1
		if linecount==1:
			continue
		strs=line.strip().split("\t")
		if strs[8]=="+":
			bp1=BP_t(strs[0], int(strs[2]), False)
		else:
			bp1=BP_t(strs[0], int(strs[1]), True)
		if strs[9]=="+":
			bp2=BP_t(strs[3], int(strs[4]), True)
		else:
			bp2=BP_t(strs[3], int(strs[5]), False)
		sv=SV_t(bp1, bp2)
		SVs.append(sv)
	fp.close()
	return SVs

def AnnotatebyGenes(SVs, Genes):
	Annot=[0]*len(SVs)
	HitGene1=[False]*len(SVs)
	HitGene2=[False]*len(SVs)
	Orientation1=[]
	Orientation2=[]
	for i in range(len(SVs)):
		Orientation1.append([])
		Orientation2.append([])
	addprefix=False
	delprefix=False
	if "chr" in Genes[0].Chr and "chr" not in SVs[0].BP1.Chr:
		addprefix=True
	if "chr" not in Genes[0].Chr and "chr" in SVs[0].BP1.Chr:
		delprefix=True
	for g in Genes:
		for i in range(len(SVs)):
			sv=SVs[i]
			bp1chr=sv.BP1.Chr
			bp2chr=sv.BP2.Chr
			if addprefix:
				bp1chr="chr"+bp1chr
				bp2chr="chr"+bp2chr
			elif delprefix:
				bp1chr=bp1chr[3:]
				bp2chr=bp2chr[3:]
			if g.Chr==bp1chr and g.StartPos<=sv.BP1.Position+500 and g.EndPos>=sv.BP1.Position-500:
				HitGene1[i]=True
				Orientation1[i].append(g.Strand)
			if g.Chr==bp2chr and g.StartPos<=sv.BP2.Position+500 and g.EndPos>=sv.BP2.Position-500:
				HitGene2[i]=True
				Orientation2[i].append(g.Strand)
	return [HitGene1, HitGene2, Orientation1, Orientation2]

def Annotate(SVs, Genes):
	thresh1=500
	thresh2=20
	Annot=[0]*len(SVs)
	CutExons=[0]*len(SVs)
	HitGene1=[False]*len(SVs)
	HitGene2=[False]*len(SVs)
	CutExon1=[False]*len(SVs)
	CutExon2=[False]*len(SVs)
	addprefix=False
	delprefix=False
	if "chr" in Genes[0].Chr and "chr" not in SVs[0].BP1.Chr:
		addprefix=True
	if "chr" not in Genes[0].Chr and "chr" in SVs[0].BP1.Chr:
		delprefix=True
	for g in Genes:
		for i in range(len(SVs)):
			sv=SVs[i]
			bp1chr=sv.BP1.Chr
			bp2chr=sv.BP2.Chr
			if addprefix:
				bp1chr="chr"+bp1chr
				bp2chr="chr"+bp2chr
			elif delprefix:
				bp1chr=bp1chr[3:]
				bp2chr=bp2chr[3:]
			if g.Chr==bp1chr and g.StartPos<=sv.BP1.Position+500 and g.EndPos>=sv.BP1.Position-500:
				HitGene1[i]=True
				for exon in g.Exons:
					if exon.StartPos<=sv.BP1.Position-thresh2 and exon.EndPos>=sv.BP1.Position+thresh2:
						CutExon1[i]=True
				# if CutExon1[i]==-1:
				# 	for j in range(len(g.Exons)):
				# 		if g.Exons[j].StartPos<=sv.BP1.Position-thresh2 and g.Exons[j].EndPos>=sv.BP1.Position+thresh2 and (j==0 or g.Exons[j-1].EndPos<g.Exons[j].StartPos) and (j+1==len(g.Exons) or g.Exons[j].EndPos<g.Exons[j+1].StartPos):
				# 			CutExon1[i]=1
				# 	if CutExon1[i]!=1:
				# 		CutExon1[i]=0
			if g.Chr==bp2chr and g.StartPos<=sv.BP2.Position+500 and g.EndPos>=sv.BP2.Position-500:
				HitGene2[i]=True
				for exon in g.Exons:
					if exon.StartPos<=sv.BP2.Position-thresh2 and exon.EndPos>=sv.BP2.Position+thresh2:
						CutExon2[i]=True
				# if CutExon2[i]==-1:
				# 	for j in range(len(g.Exons)):
				# 		if g.Exons[j].StartPos<=sv.BP2.Position-thresh2 and g.Exons[j].EndPos>=sv.BP2.Position+thresh2 and (j==0 or g.Exons[j-1].EndPos<g.Exons[j].StartPos) and (j+1==len(g.Exons) or g.Exons[j].EndPos<g.Exons[j+1].StartPos):
				# 			CutExon2[i]=1
				# 	if CutExon2[i]!=1:
				# 		CutExon2[i]=0
	for i in range(len(SVs)):
		if HitGene1[i] and HitGene2[i]:
			Annot[i]=3
		elif HitGene1[i]:
			Annot[i]=1
		elif HitGene2[i]:
			Annot[i]=2
		else:
			Annot[i]=0
	for i in range(len(SVs)):
		if CutExon1[i] and CutExon2[i]:
			CutExons[i]=3
		elif CutExon1[i]:
			CutExons[i]=1
		elif CutExon2[i]:
			CutExons[i]==2
		else:
			CutExons[i]==0
	return [Annot, CutExons]

# def WriteAnnot(SVs, Annot, CutExons, outputfile):
# 	fp=open(outputfile, 'w')
# 	fp.write("# Chr1\tPos1\tIsLeft1\tChr2\tPos2\tIsLeft2\tBP1hitGene\tBP2hitGene\tBP1cutExon\tBP2cutExon\n")
# 	for i in range(len(SVs)):
# 		sv=SVs[i]
# 		bp1Annot=(Annot[i]%2==1)
# 		bp2Annot=(Annot[i]>1)
# 		bp1Cut=(CutExons[i]%2==1)
# 		bp2Cut=(CutExons[i]>1)
# 		fp.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sv.BP1.Chr, sv.BP1.Position, sv.BP1.IsLeft, sv.BP2.Chr, sv.BP2.Position, sv.BP2.IsLeft, bp1Annot, bp2Annot, bp1Cut, bp2Cut))
# 	fp.close()

def WriteAnnot(SVs, HitGene1, HitGene2, Orientation1, Orientation2, outputfile):
	fp=open(outputfile, 'w')
	fp.write("# Chr1\tPos1\tIsLeft1\tChr2\tPos2\tIsLeft2\tBP1hitGene\tBP2hitGene\tBP1Orient\tBP2Orient\tIsFusionGene\n")
	for i in range(len(SVs)):
		sv=SVs[i]
		ori1="".join(Orientation1[i])
		if ori1=="":
			ori1="."
		ori2="".join(Orientation2[i])
		if ori2=="":
			ori2="."
		isfg=False
		if HitGene1[i] and HitGene2[i]:
			if sv.BP1.IsLeft==sv.BP2.IsLeft and (('-' in ori1 and '+' in ori2) or ('+' in ori1 and '-' in ori2)):
				isfg=True
			elif sv.BP1.IsLeft!=sv.BP2.IsLeft and (('-' in ori1 and '-' in ori2) or ('+' in ori1 and '+' in ori2)):
				isfg=True
		fp.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sv.BP1.Chr, sv.BP1.Position, sv.BP1.IsLeft, sv.BP2.Chr, sv.BP2.Position, sv.BP2.IsLeft, HitGene1[i], HitGene2[i], ori1, ori2, isfg))
	fp.close()

if __name__=="__main__":
	if len(sys.argv)<5:
		print("python AnnotTrueSV.py <choice from 1-3> <TrueSV/PredSV> <GTFfile> <Output> (<RNAbam>)")
		print("\t1: True SV")
		print("\t2: SQUID_discordantedges.txt")
		print("\t3: FusionCatcher")
		print("\t4: JAFFA")
		print("\t5: Defuse")
		print("\t6: Chimerascan")
		print("\t7: INTEGRATE")
	else:
		SVFile=sys.argv[2]
		GTFfile=sys.argv[3]
		Output=sys.argv[4]

		if sys.argv[1]=='1':
			SVs=ReadTrueSV(SVFile)
		elif sys.argv[1]=='2':
			RNAbam=sys.argv[5]
			[RefTable, RefLength, RefName]=BuildReference(RNAbam)
			SVs=ReadSQUID(SVFile, RefName)
		elif sys.argv[1]=='3':
			SVs=ReadFusionCatcher(SVFile)
		elif sys.argv[1]=='4':
			SVs=ReadJAFFA(SVFile)
		elif sys.argv[1]=='5':
			SVs=ReadDefuse(SVFile)
		elif sys.argv[1]=='6':
			SVs=ReadChimerascan(SVFile)
		elif sys.argv[1]=='7':
			SVs=ReadIntegrate(SVFile)
		# [Annot, Orientation]=AnnotatebyGTF(SVs, GTFfile)
		Genes=ReadGTF(GTFfile)
		[HitGene1, HitGene2, Orientation1, Orientation2]=AnnotatebyGenes(SVs, Genes)
		WriteAnnot(SVs, HitGene1, HitGene2, Orientation1, Orientation2, Output)
