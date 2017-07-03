#!/bin/bash


# STAR indexing and aligning
if [ $(file sampledata/genome_rearranged.fa.gz) = "gzip" ]; then
	gunzip sampledata/genome_rearranged.fa.gz
fi

mkdir -p sampledata/STARindex

STAR --runThreadN 4 --runMode genomeGenerate --genomeDir sampledata/STARindex --genomeFastaFiles sampledata/genome_rearranged.fa
mv Log.out sampledata/STARindex/

mkdir -p sampledata/StarAlign
STAR --runThreadN 4 --genomeDir sampledata/STARindex/ --readFilesIn sampledata/RNA1.fastq.gz sampledata/RNA2.fastq.gz --outFileNamePrefix sampledata/StarAlign/ --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --chimSegmentMin 20 --outSAMstrandField intronMotif
samtools view -Shb sampledata/StarAlign/Chimeric.out.sam sampledata/StarAlign/Chimeric.out.bam


# SQUID predicting
./../bin/squid -b sampledata/StarAlign/Aligned.sortedByCoord.out.bam -c sampledataStarAlign/Chimeric.out.bam -G 1 -CO 1 -o squid
