#!/bin/python

import sys

def ReverseComp(line):
    newline=''
    for i in line:
        if i.upper()=='A':
            newline+='T'
        elif i.upper()=='T':
            newline+='A'
        elif i.upper()=='C':
            newline+='G'
        elif i.upper()=='G':
            newline+='C'
        else:
            newline+=i.upper()
    newline=newline[::-1]
    return newline

def ReadGenome(fafile):
    genome={}
    fp=open(fafile,'r')
    line=fp.readline().strip()
    tmpseq=''
    tmpname=''
    while line!='':
        if line[0]=='>':
            if len(tmpseq)!=0:
                genome[tmpname]=tmpseq
            tmpseq=''
            strs=line.split()
            tmpname=strs[0][1:]
        else:
            tmpseq+=line
        line=fp.readline().strip()
    genome[tmpname]=tmpseq
    fp.close()
    return genome

def WriteChr(outdir, genome):
    for k in genome.keys():
        fp=open(outdir+"/"+str(k)+".fa", 'w')
        fp.write(">"+str(k)+"\n")
        fp.write(genome[k]+"\n")
        fp.close()

if __name__=="__main__":
    InputGenome=sys.argv[1]
    OutputDir=sys.argv[2]
    genome=ReadGenome(InputGenome)
    WriteChr(OutputDir, genome)
