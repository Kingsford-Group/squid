#!/bin/bash


# STAR indexing and aligning
if [[ -e sampledata/genome_rearranged.fa.gz ]] && [[ ! -e sampledata/genome_rearranged.fa ]]; then
	gunzip sampledata/genome_rearranged.fa.gz
fi

if [[ ! -e sampledata/genome_rearranged.fa ]]; then
	echo "No such files: sampledata/genome_rearranged.fa.gz or sampledata/genome_rearranged.fa"
	exit 1
fi

mkdir -p sampledata/STARindex

STAR --runThreadN 4 --runMode genomeGenerate --genomeDir sampledata/STARindex --genomeFastaFiles sampledata/genome_rearranged.fa
mv Log.out sampledata/STARindex/

mkdir -p sampledata/StarAlign
STAR --runThreadN 4 --genomeDir sampledata/STARindex/ --readFilesIn sampledata/RNA1.fastq.gz sampledata/RNA2.fastq.gz --readFilesCommand zcat --outFileNamePrefix sampledata/StarAlign/ --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --chimSegmentMin 20 --outSAMstrandField intronMotif --chimOutType SeparateSAMold
samtools view -Shb sampledata/StarAlign/Chimeric.out.sam -o sampledata/StarAlign/Chimeric.out.bam


# SQUID predicting
./../bin/squid -b sampledata/StarAlign/Aligned.sortedByCoord.out.bam -c sampledata/StarAlign/Chimeric.out.bam -G 1 -CO 1 -o squid
