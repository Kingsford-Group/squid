#!/bin/python

import sys
import numpy as np

class Transcript_t(object):
	def __init__(self, _TransID, _GeneID, _GeneName, _Chr, _Strand, _StartPos, _EndPos):
		self.TransID=_TransID
		self.GeneID=_GeneID
		self.GeneName=_GeneName
		self.Chr = _Chr
		self.Strand = _Strand
		self.StartPos = _StartPos
		self.EndPos = _EndPos
		self.Exons = [] # each exon is a tuple of two integers
	def __eq__(self, other):
		if isinstance(other, Transcript_t):
			return (self.Chr==other.Chr and self.Strand==other.Strand and len(self.Exons)==len(other.Exons) and \
				min([self.Exons[i]==other.Exons[i] for i in range(len(self.Exons))])!=0)
		return NotImplemented
	def __ne__(self, other):
		result=self.__eq__(other)
		if result is NotImplemented:
			return result
		return not result
	def __lt__(self, other):
		if isinstance(other, Transcript_t):
			if self.Chr!=other.Chr:
				return self.Chr<other.Chr
			elif self.StartPos!=other.StartPos:
				return self.StartPos<other.StartPos
			else:
				return self.EndPos<other.EndPos
		return NotImplemented
	def __gt__(self, other):
		if isinstance(other, Transcript_t):
			if self.Chr!=other.Chr:
				return self.Chr>other.Chr
			elif self.StartPos!=other.StartPos:
				return self.StartPos>other.StartPos
			else:
				return self.EndPos<other.EndPos
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


def GetFeature(line, key):
	s=line.index(key)
	t=line.index(";", s+1)
	return line[(s+len(key)+2):(t-1)]


def ReadGTF(gtffile, key_gene_id = "gene_id", key_gene_symbol = "gene_name"):
	Transcripts={}
	strand=""
	fp=open(gtffile, 'r')
	tmptransname=""
	tmptranscript=None
	extraExons = []
	for line in fp:
		if line[0]=='#':
			continue
		strs=line.strip().split("\t")
		if strs[2]=="transcript":
			if tmptransname!="" and not (tmptranscript is None):
				Transcripts[tmptransname]=tmptranscript
			if not "transcript_id" in line:
				print("GTF file attribute column doesn't contain transcript_id: " + line)
				sys.exit()
			if not key_gene_id in line:
				print("GTF file attribute column doesn't contain "+key_gene_id+": " + line)
				sys.exit()
			if not key_gene_symbol in line:
				print("GTF file attribute column doesn't contain "+key_gene_symbol+": " + line)
				sys.exit()
			tmptransname=GetFeature(line, "transcript_id")
			tmpgeneid=GetFeature(line, key_gene_id)
			tmpgenename=GetFeature(line, key_gene_symbol)
			tmptranscript=Transcript_t(tmptransname, tmpgeneid, tmpgenename, strs[0], (strs[6]=="+"), int(strs[3])-1, int(strs[4]))
		elif strs[2]=="exon":
			thistransid=GetFeature(line, "transcript_id")
			if not "transcript_id" in line:
				print("GTF file attribute column doesn't contain transcript_id: " + line)
				sys.exit()
			if thistransid == tmptransname and not (tmptranscript is None):
				tmptranscript.Exons.append((int(strs[3])-1, int(strs[4])))
			else:
				extraExons.append([thistransid, int(strs[3])-1, int(strs[4])])
	if tmptransname!="" and not (tmptranscript is None):
		Transcripts[tmptransname]=tmptranscript
	for e in extraExons:
		assert(e[0] in Transcripts)
		Transcripts[e[0]].Exons.append((e[1],e[2]))
	for t in Transcripts:
		Transcripts[t].Exons.sort(key=lambda x:x[0])
		if not Transcripts[t].Strand:
			Transcripts[t].Exons = Transcripts[t].Exons[::-1]
	fp.close()
	return Transcripts


def Map_Gene_Trans(Transcripts):
	GeneTransMap={}
	TransGeneMap={}
	for v in Transcripts.values():
		TransGeneMap[v.TransID]=v.GeneID
		if v.GeneID in GeneTransMap:
			GeneTransMap[v.GeneID].append(v.TransID)
		else:
			GeneTransMap[v.GeneID]=[v.TransID]
	for g,v in GeneTransMap.items():
		sortedv = sorted(v)
		GeneTransMap[g] = sortedv
	return [GeneTransMap, TransGeneMap]


def GetTransLength(Transcripts):
	TransLength = {t:np.sum([e[1]-e[0] for e in v.Exons]) for t,v in Transcripts.items()}
	return TransLength


def ReadTranscriptFasta(filename, namesplitter = " "):
	TransSequence = {}
	fp = open(filename, 'r')
	name = ""
	seq = ""
	for line in fp:
		if line[0] == '>':
			if len(name) != 0:
				TransSequence[name] = seq
			name = line.strip().split(namesplitter)[0][1:]
			seq = ""
		else:
			seq += line.strip()
	if len(name) != "":
		TransSequence[name] = seq
	fp.close()
	return TransSequence


class GeneLocater(object):
	def __init__(self, Transcripts, GeneTransMap):
		# attributes
		self.GeneRanges = []
		self.GeneExons = []
		self.GeneNames = []
		self.GeneIndex = {}
		# construct these attributes
		for g,v in GeneTransMap.items():
			chrnames = []
			lowerbounds = []
			upperbounds = []
			for t in v:
				chrnames.append(Transcripts[t].Chr)
				lowerbounds.append(Transcripts[t].StartPos)
				upperbounds.append(Transcripts[t].EndPos)
			assert(len(set(chrnames)) == 1)
			lb = np.min(np.array(lowerbounds))
			ub = np.max(np.array(upperbounds))
			self.GeneRanges.append( (chrnames[0], lb, ub) )
			self.GeneNames.append(g)
			tmp_exons = sum([Transcripts[t].Exons for t in v], [])
			tmp_exons.sort()
			union_exons = []
			for e in tmp_exons:
				if len(union_exons) != 0 and union_exons[-1][1] > e[0]:
					union_exons[-1][1] = max(e[1], union_exons[-1][1])
			self.GeneExons.append(tmp_exons)
		# sort by gene locations
		assert(len(self.GeneNames) == len(self.GeneRanges))
		indexes = list(range(len(self.GeneNames)))
		indexes.sort(key = lambda o:self.GeneRanges[o])
		self.GeneRanges = [self.GeneRanges[o] for o in indexes]
		self.GeneNames = [self.GeneNames[o] for o in indexes]
		self.GeneExons = [self.GeneExons[o] for o in indexes]
		self.GeneIndex = {self.GeneNames[o]:o for o in range(len(self.GeneNames))}

	def LocatePosition_generange(self, chr, pos, window = 100000, fuzzy = 50):
		genes = []
		# binary search
		low = 0
		high = len(self.GeneNames)
		while low < high:
			mid = int((low + high) / 2)
			if self.GeneRanges[mid][0] < chr or (self.GeneRanges[mid][0] == chr and self.GeneRanges[mid][2] < pos - fuzzy):
				low = mid + 1
			elif self.GeneRanges[mid][0] == chr and self.GeneRanges[mid][1] <= pos + fuzzy and self.GeneRanges[mid][2] > pos - fuzzy:
				low = mid
				high = mid
			else:
				high = mid - 1
		count_low = 0
		count_high = 0
		if low >= 0 and low != len(self.GeneNames):
			while low >= 0 and ((count_low < 20) or (self.GeneRanges[low][0] == chr and self.GeneRanges[low][2] + window > pos)):
				count_low += 1
				if self.GeneRanges[low][0] == chr and self.GeneRanges[low][1] <= pos + fuzzy and self.GeneRanges[low][2] > pos - fuzzy:
					genes.append(self.GeneNames[low])
				low -= 1
		if high >= 0 and high != len(self.GeneNames):
			while high < len(self.GeneNames) and ((count_high < 20) or (self.GeneRanges[high][0] == chr and self.GeneRanges[high][1] <= pos + fuzzy)):
				count_high += 1
				if self.GeneRanges[high][0] == chr and self.GeneRanges[high][1] <= pos + fuzzy and self.GeneRanges[high][2] > pos - fuzzy:
					genes.append(self.GeneNames[high])
				high += 1
		return list(set(genes))

	def LocatePosition_exonrange(self, chr, pos, window = 100000, fuzzy = 50):
		potential_genes = self.LocatePosition_generange(chr, pos, window, fuzzy)
		final_genes = []
		for g in potential_genes:
			exons = self.GeneExons[self.GeneIndex[g]]
			for e in exons:
				if e[0] <= pos + fuzzy and e[1] > pos - fuzzy:
					final_genes.append(g)
					break
		return final_genes


def Annotate(insquidfile, outputfile, glocater, Transcripts, GeneTransMap):
	fpin = open(insquidfile, 'r')
	fpout = open(outputfile, 'w')
	for line in fpin:
		strs = line.strip().split("\t")
		if line[0] == '#':
			fpout.write("\t".join(strs[:10]) + "\tType\tFusedGenes\n")
		else:
			# extract breakpoint position information
			chr1 = strs[0]
			chr2 = strs[3]
			bp1 = int(strs[1])
			bp2 = int(strs[4])
			bp1strand = (strs[8] == '+')
			bp2strand = (strs[9] == '+')
			if strs[8] == '+':
				bp1 = int(strs[2])
			if strs[9] == '+':
				bp2 = int(strs[5])
			# locate corresponding genes of both breakpoint
			genes1 = glocater.LocatePosition_generange(chr1, bp1)
			genes2 = glocater.LocatePosition_generange(chr2, bp2)
			# find valid fusion genes
			FusionPairs = []
			for g1 in genes1:
				for g2 in genes2:
					strand1 = Transcripts[GeneTransMap[g1][0]].Strand
					strand2 = Transcripts[GeneTransMap[g2][0]].Strand
					# if chr1=="8":
					# 	print([g1, bp1, strand1, bp1strand, g2, bp2, strand2, bp2strand, (strand1 == bp1strand), (strand2 == bp2strand)])
					# check whether strand info is valid for fusion-gene
					# in order to be a fusion-gene, one breakpoint should agree with the gene's strand, and the other should be the opposite
					if (strand1 == bp1strand) != (strand2 == bp2strand):
						# 5' gene to be the first; 3' gene to be the second
						if (strand1 == bp1strand):
							FusionPairs.append( Transcripts[GeneTransMap[g1][0]].GeneName+":"+Transcripts[GeneTransMap[g2][0]].GeneName )
						else:
							FusionPairs.append( Transcripts[GeneTransMap[g2][0]].GeneName+":"+Transcripts[GeneTransMap[g1][0]].GeneName )
			if len(FusionPairs) == 0:
				fpout.write("\t".join(strs[:10]) + "\tnon-fusion-gene\t.\n" )
			else:
				fpout.write("\t".join(strs[:10]) + "\tfusion-gene\t" + ",".join(FusionPairs) + "\n" )
	fpin.close()
	fpout.close()


def ParseArgument(argv):
	key_gene_id = "gene_id"
	key_gene_symbol = "gene_name"
	GTFfile = ""
	SquidPrediction = ""
	OutputFile = ""
	i = 1
	while i < len(argv):
		if argv[i] == "--geneid":
			if i+1 >= len(argv) or argv[i+1][:2] == '--':
				print("GTF gene ID attribute string is empty!")
				sys.exit()
			key_gene_id = argv[i+1]
			i += 2
		elif argv[i] == "--genesymbol":
			if i+1 >= len(argv) or argv[i+1][:2] == '--':
				print("GTF gene symbol attribute string is empty!")
				sys.exit()
			key_gene_symbol = argv[i+1]
			i += 2
		elif argv[i][:2] == "--":
			print("Unknown argument "+argv[i])
			sys.exit()
		else:
			if i+2 >= len(argv):
				print([i, len(argv)])
				print("Missing GTFfile or SquidPrediction or OutputFile")
				sys.exit()
			else:
				GTFfile = argv[i]
				SquidPrediction = argv[i+1]
				OutputFile = argv[i+2]
				break
	return [key_gene_id, key_gene_symbol, GTFfile, SquidPrediction, OutputFile]


if __name__=="__main__":
	if len(sys.argv) == 1:
		print("python AnnotateSQUIDoutput.py [options] <GTFfile> <SquidPrediction> <OutputFile>")
		print("options:")
		print("\t--geneid\tstring\tGTF gene ID attribute string, the attribute name in GTF record that corresponds to the gene ID (default: gene_id)")
		print("\t--genesymbol\tstring\tGTF gene symbol attribute string, the attribute name in GTF record that corresponds to the gene symbol (default: gene_name)")
	else:
		[key_gene_id, key_gene_symbol, GTFfile, SquidPrediction, OutputFile] = ParseArgument(sys.argv)

		Transcripts = ReadGTF(GTFfile, key_gene_id, key_gene_symbol)
		[GeneTransMap, TransGeneMap] = Map_Gene_Trans(Transcripts)
		glocater = GeneLocater(Transcripts, GeneTransMap)
		Annotate(SquidPrediction, OutputFile, glocater, Transcripts, GeneTransMap)
