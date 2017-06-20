#!/bin/python

import sys 

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

def ReadTrueSV(filename):
	fp=open(filename, 'r')
	SVs=[]
	for line in fp:
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

def SQUIDAccuracy(inputfile, outputfile, SVs):
	global thresh
	Hits=[0]*len(SVs)
	addprefix=False
	delprefix=False
	linecount=0
	fpin=open(inputfile, 'r')
	fpout=open(outputfile, 'w')
	hitcount=0
	for line in fpin:
		if line[0]=='#':
			continue
		linecount+=1
		strs=line.strip().split("\t")
		if linecount==1:
			if "chr" in strs[0] and "chr" not in SVs[0].BP1.Chr:
				delprefix=True
			if "chr" not in strs[0] and "chr" in SVs[0].BP1.Chr:
				addprefix=True
		bp1=BP_t(strs[0], int(strs[1]), (strs[8]=='-'))
		bp2=BP_t(strs[3], int(strs[4]), (strs[9]=='-'))
		if not bp1.IsLeft:
			bp1.Position=int(strs[2])
		if not bp2.IsLeft:
			bp2.Position=int(strs[5])
		if addprefix:
			bp1.Chr="chr"+bp1.Chr
			bp2.Chr="chr"+bp2.Chr
		if delprefix:
			bp1.Chr=bp1.Chr[3:]
			bp2.Chr=bp2.Chr[3:]
		sqsv=SV_t(bp1, bp2)

		hit=False
		for i in range(len(SVs)):
			sv=SVs[i]
			if sqsv.BP1.Chr==sv.BP1.Chr and sqsv.BP2.Chr==sv.BP2.Chr and abs(sqsv.BP1.Position-sv.BP1.Position)<thresh and abs(sqsv.BP2.Position-sv.BP2.Position)<thresh and sqsv.BP1.IsLeft==sv.BP1.IsLeft and sqsv.BP2.IsLeft==sv.BP2.IsLeft:
				hit=True
				Hits[i]=1
				break
		if hit:
			hitcount+=1
			fpout.write("TP\t"+line)
		else:
			fpout.write("FP\t"+line)
	fpin.close()
	fpout.close()
	print("RawAccuracy = {}\tRawSensitivity = {}".format(hitcount, sum(Hits)))


def JAFFAaccuracy(inputfile, outputfile, SVs):
	global thresh
	Hits=[0]*len(SVs)
	addprefix=False
	delprefix=False
	linecount=0
	fpin=open(inputfile, 'r')
	fpout=open(outputfile, 'w')
	hitcount=0
	for line in fpin:
		linecount+=1
		if linecount==1:
			continue
		strs=line.replace("\"", "").strip().split(",")
		if linecount==2:
			if "chr" in strs[2] and "chr" not in SVs[0].BP1.Chr:
				delprefix=True
			if "chr" not in strs[2] and "chr" in SVs[0].BP1.Chr:
				addprefix=True
		bp1=BP_t(strs[2], int(strs[3]), (strs[4]=='-'))
		bp2=BP_t(strs[5], int(strs[6]), (strs[7]=='+'))
		if addprefix:
			bp1.Chr="chr"+bp1.Chr
			bp2.Chr="chr"+bp2.Chr
		if delprefix:
			bp1.Chr=bp1.Chr[3:]
			bp2.Chr=bp2.Chr[3:]
		jasv=SV_t(bp1, bp2)

		hit=False
		for i in range(len(SVs)):
			sv=SVs[i]
			if jasv.BP1.Chr==sv.BP1.Chr and jasv.BP2.Chr==sv.BP2.Chr and abs(jasv.BP1.Position-sv.BP1.Position)<thresh and abs(jasv.BP2.Position-sv.BP2.Position)<thresh and jasv.BP1.IsLeft==sv.BP1.IsLeft and jasv.BP2.IsLeft==sv.BP2.IsLeft:
				hit=True
				Hits[i]=1
				break
		if hit:
			hitcount+=1
			fpout.write("TP,"+line)
		else:
			fpout.write("FP,"+line)
	fpin.close()
	fpout.close()
	print("RawAccuracy = {}\tRawSensitivity = {}".format(hitcount, sum(Hits)))


def FCaccuracy(inputfile, outputfile, SVs):
	global thresh
	Hits=[0]*len(SVs)
	addprefix=False
	delprefix=False
	linecount=0
	fpin=open(inputfile, 'r')
	fpout=open(outputfile, 'w')
	hitcount=0
	for line in fpin:
		linecount+=1
		if linecount==1:
			continue
		strs=line.strip().split()
		if linecount==2:
			if "chr" in strs[8] and "chr" not in SVs[0].BP1.Chr:
				delprefix=True
			if "chr" not in strs[8] and "chr" in SVs[0].BP1.Chr:
				addprefix=True
		if ":" in strs[7]:
			bp1=BP_t(strs[7].split(":")[0], int(strs[7].split(":")[1]), (strs[7].split(":")[2]=='-'))
			bp2=BP_t(strs[8].split(":")[0], int(strs[8].split(":")[1]), (strs[8].split(":")[2]=='+'))
		else:
			bp1=BP_t(strs[8].split(":")[0], int(strs[8].split(":")[1]), (strs[8].split(":")[2]=='-'))
			bp2=BP_t(strs[9].split(":")[0], int(strs[9].split(":")[1]), (strs[9].split(":")[2]=='+'))
		if addprefix:
			bp1.Chr="chr"+bp1.Chr
			bp2.Chr="chr"+bp2.Chr
		if delprefix:
			bp1.Chr=bp1.Chr[3:]
			bp2.Chr=bp2.Chr[3:]
		fcsv=SV_t(bp1, bp2)

		hit=False
		for i in range(len(SVs)):
			sv=SVs[i]
			if fcsv.BP1.Chr==sv.BP1.Chr and fcsv.BP2.Chr==sv.BP2.Chr and abs(fcsv.BP1.Position-sv.BP1.Position)<thresh and abs(fcsv.BP2.Position-sv.BP2.Position)<thresh and fcsv.BP1.IsLeft==sv.BP1.IsLeft and fcsv.BP2.IsLeft==sv.BP2.IsLeft:
				hit=True
				Hits[i]=1
				break
		if hit:
			hitcount+=1
			fpout.write("TP\t"+line)
		else:
			fpout.write("FP\t"+line)
	fpin.close()
	fpout.close()
	print("RawAccuracy = {}\tRawSensitivity = {}".format(hitcount, sum(Hits)))


def DefuseAccuracy(inputfile, outputfile, SVs):
	global thresh
	Hits=[0]*len(SVs)
	addprefix=False
	delprefix=False
	linecount=0
	fpin=open(inputfile, 'r')
	fpout=open(outputfile, 'w')
	hitcount=0
	for line in fpin:
		linecount+=1
		if linecount==1:
			continue
		strs=line.strip().split("\t")
		if linecount==2:
			if "chr" in strs[24] and "chr" not in SVs[0].BP1.Chr:
				delprefix=True
			if "chr" not in strs[24] and "chr" in SVs[0].BP1.Chr:
				addprefix=True
		bp1=BP_t(strs[24], int(strs[37]), (strs[43]=='-'))
		bp2=BP_t(strs[25], int(strs[38]), (strs[44]=='-'))
		if addprefix:
			bp1.Chr="chr"+bp1.Chr
			bp2.Chr="chr"+bp2.Chr
		if delprefix:
			bp1.Chr=bp1.Chr[3:]
			bp2.Chr=bp2.Chr[3:]
		dfsv=SV_t(bp1, bp2)

		hit=False
		for i in range(len(SVs)):
			sv=SVs[i]
			if dfsv.BP1.Chr==sv.BP1.Chr and dfsv.BP2.Chr==sv.BP2.Chr and abs(dfsv.BP1.Position-sv.BP1.Position)<thresh and abs(dfsv.BP2.Position-sv.BP2.Position)<thresh and dfsv.BP1.IsLeft==sv.BP1.IsLeft and dfsv.BP2.IsLeft==sv.BP2.IsLeft:
				hit=True
				Hits[i]=1
				break
		if hit:
			hitcount+=1
			fpout.write("TP\t"+line)
		else:
			fpout.write("FP\t"+line)
	fpin.close()
	fpout.close()
	print("RawAccuracy = {}\tRawSensitivity = {}".format(hitcount, sum(Hits)))


def IntegrateAccuracy(inputfile, outputfile, SVs):
	global thresh
	Hits=[0]*len(SVs)
	addprefix=False
	delprefix=False
	linecount=0
	fpin=open(inputfile, 'r')
	fpout=open(outputfile, 'w')
	hitcount=0
	for line in fpin:
		linecount+=1
		strs=line.strip().split("\t")
		if linecount==1:
			if "chr" in strs[0] and "chr" not in SVs[0].BP1.Chr:
				delprefix=True
			if "chr" not in strs[0] and "chr" in SVs[0].BP1.Chr:
				addprefix=True
		if strs[8]=="+":
			bp1=BP_t(strs[0], int(strs[2]), False)
		else:
			bp1=BP_t(strs[0], int(strs[1]), True)
		if strs[9]=="+":
			bp2=BP_t(strs[3], int(strs[4]), True)
		else:
			bp2=BP_t(strs[3], int(strs[5]), False)
		if addprefix:
			bp1.Chr="chr"+bp1.Chr
			bp2.Chr="chr"+bp2.Chr
		if delprefix:
			bp1.Chr=bp1.Chr[3:]
			bp2.Chr=bp2.Chr[3:]
		insv=SV_t(bp1, bp2)

		hit=False
		for i in range(len(SVs)):
			sv=SVs[i]
			if insv.BP1.Chr==sv.BP1.Chr and insv.BP2.Chr==sv.BP2.Chr and abs(insv.BP1.Position-sv.BP1.Position)<thresh and abs(insv.BP2.Position-sv.BP2.Position)<thresh and insv.BP1.IsLeft==sv.BP1.IsLeft and insv.BP2.IsLeft==sv.BP2.IsLeft:
				hit=True
				Hits[i]=1
				break
		if hit:
			hitcount+=1
			fpout.write("TP\t"+line)
		else:
			fpout.write("FP\t"+line)
	fpin.close()
	fpout.close()
	print("RawAccuracy = {}\tRawSensitivity = {}".format(hitcount, sum(Hits)))


def ChimerascanAccuracy(inputfile, outputfile, SVs):
	global thresh
	Hits=[0]*len(SVs)
	addprefix=False
	delprefix=False
	linecount=0
	fpin=open(inputfile, 'r')
	fpout=open(outputfile, 'w')
	hitcount=0
	for line in fpin:
		linecount+=1
		if linecount==1:
			continue
		strs=line.strip().split("\t")
		if linecount==2:
			if "chr" in strs[0] and "chr" not in SVs[0].BP1.Chr:
				delprefix=True
			if "chr" not in strs[0] and "chr" in SVs[0].BP1.Chr:
				addprefix=True
		if strs[8]=="-":
			bp1=BP_t(strs[0], int(strs[1]), True)
		else:
			bp1=BP_t(strs[0], int(strs[2]), False)
		if strs[9]=='+':
			bp2=BP_t(strs[3], int(strs[4]), True)
		else:
			bp2=BP_t(strs[3], int(strs[5]), False)
		if addprefix:
			bp1.Chr="chr"+bp1.Chr
			bp2.Chr="chr"+bp2.Chr
		if delprefix:
			bp1.Chr=bp1.Chr[3:]
			bp2.Chr=bp2.Chr[3:]
		cssv=SV_t(bp1, bp2)

		hit=False
		for i in range(len(SVs)):
			sv=SVs[i]
			if cssv.BP1.Chr==sv.BP1.Chr and cssv.BP2.Chr==sv.BP2.Chr and abs(cssv.BP1.Position-sv.BP1.Position)<thresh and abs(cssv.BP2.Position-sv.BP2.Position)<thresh and cssv.BP1.IsLeft==sv.BP1.IsLeft and cssv.BP2.IsLeft==sv.BP2.IsLeft:
				hit=True
				Hits[i]=1
				break
		if hit:
			hitcount+=1
			fpout.write("TP\t"+line)
		else:
			fpout.write("FP\t"+line)
	fpin.close()
	fpout.close()
	print("RawAccuracy = {}\tRawSensitivity = {}".format(hitcount, sum(Hits)))

if __name__=="__main__":
	thresh=30000
	if len(sys.argv)!=4:
		print("python3  ValidFusionGene.py <option from 1-6>  <TrueSV.bedpe>  <FusionCatcher_pred.txt>")
		print("\t1: squid")
		print("\t2: fusioncatcher")
		print("\t3: jaffa")
		print("\t4: defuse")
		print("\t5: chimerascan")
		print("\t6: integrate")
	else:
		TrueSV=sys.argv[2]
		pred=sys.argv[3]
		start=pred.rfind(".")
		out=pred[:start]+"_hit"+pred[start:]

		SVs=ReadTrueSV(TrueSV)
		if sys.argv[1]=="1":
			SQUIDAccuracy(pred, out, SVs)
		elif sys.argv[1]=="2":
			FCaccuracy(pred, out, SVs)
		elif sys.argv[1]=="3":
			JAFFAaccuracy(pred, out, SVs)
		elif sys.argv[1]=="4":
			DefuseAccuracy(pred, out, SVs)
		elif sys.argv[1]=="5":
			ChimerascanAccuracy(pred, out, SVs)
		elif sys.argv[1]=="6":
			IntegrateAccuracy(pred, out, SVs)
