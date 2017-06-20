#!/bin/bash

FusionCatcherDataDir= # reference data directory for fusioncatcher
ChimerascanInstallDir= # chimerascan installation directory
IntegrateInstallDir= # INTEGRATE installation directory
OutputDir= # output files for SQUID and fusion gene detection tools on HCC1954/1395 cell lines
ExeDir=$(pwd)

# download genome fasta and annotation gtf
mkdir -p ${OutputDir}/Ensemble75
cd ${OutputDir}/Ensemble75
wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
gunzip Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.75.gtf.gz
GenomeFasta=${OutputDir}/Ensemble75/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa
AnnotationGTF=${OutputDir}/Ensemble75/Homo_sapiens.GRCh37.75.gtf.gz

# preparing STAR index
mkdir -p ${GenomeFasta%/*}/STARIndex
STAR --runThreadN 12 --runMode genomeGenerate --genomeDir ${GenomeFasta%/*}/STARIndex --genomeFastaFiles ${GenomeFasta} --sjdbGTFfile ${AnnotationGTF}
# preparing chimerascan Bowtie index
mkdir -p ${OutputDir}/Ensemble75/CSBowtieIndex
gunzip ${ExeDir}/dataHCC/Homo_sapiens.GRCh37.75_chimerascan.txt.gz
python ${ChimerascanInstallDir}/bin/chimerascan_index.py ${GenomeFasta} ${ExeDir}/dataHCC/Homo_sapiens.GRCh37.75_chimerascan.txt.gz ${OutputDir}/Ensemble75/CSBowtieIndex
# preparing INTEGRATE index
gunzip ${ExeDir}/dataHCC/annot.ensembl.txt.gz
mkdir -p ${OutputDir}/Ensemble75/IntegrateIndex
${IntegrateIntallDir}/bin/Integrate mkbwt -dir ${OutputDir}/Ensemble75/IntegrateIndex ${GenomeFasta}

# download fastq file from SRA using fastq-dump
cd ${ExeDir}/dataHCC
fastq-dump --split-files SRR2532344 # HCC1954
fastq-dump --split-files SRR925710 # HCC1954
fastq-dump --split-files SRR2532336 # HCC1395

cat SRR2532344_1.fastq SRR925710_1.fastq > RNAHCC1954_1.fastq
cat SRR2532344_2.fastq SRR925710_2.fastq > RNAHCC1954_2.fastq

cd $OutputDir
# align reads using STAR
mkdir -p StarAlign_1954
STAR  --runThreadN 12 --genomeDir ${GenomeFasta%/*}/STARIndex --readFilesIn ${ExeDir}/dataHCC/RNAHCC1954_1.fastq ${ExeDir}/dataHCC/RNAHCC1954_2.fastq --outFileNamePrefix $OutputDir/StarAlign_1954/ --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --chimSegmentMin 15 --outReadsUnmapped Fastx
samtools view -Shb StarAlign_1954/Chimeric.out.sam -o StarAlign_1954/Chimeric.out.bam

mkdir -p StarAlign_1395
STAR  --runThreadN 12 --genomeDir ${GenomeFasta%/*}/STARIndex --readFilesIn ${ExeDir}/dataHCC/SRR2532336_1.fastq ${ExeDir}/dataHCC/SRR2532336_2.fastq --outFileNamePrefix $OutputDir/StarAlign_1395/ --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --chimSegmentMin 15 --outReadsUnmapped Fastx
samtools view -Shb StarAlign_1395/Chimeric.out.sam -o StarAlign_1395/Chimeric.out.bam

# SQUID
mkdir -p squid_1954
squid -b StarAlign_1954/Aligned.sortedByCoord.out.bam -c StarAlign_1954/Chimeric.out.bam -o squid_1954/squidtsv
python3 VerifyFusionGene.py

mkdir squid_1395
squid -b StarAlign_1395/Aligned.sortedByCoord.out.bam -c StarAlign_1395/Chimeric.out.bam -o squid_1395/squidtsv

# fusioncatcher
mkdir -p fusioncatcher_1954
mkdir -p fusioncatcher_1954/data
ln -s ${ExeDir}/dataHCC/RNAHCC1954_1.fastq fusioncatcher_1954/data/
ln -s ${ExeDir}/dataHCC/RNAHCC1954_2.fastq fusioncatcher_1954/data/
fusioncatcher -d ${FusionCatcherDataDir} -i fusioncatcher_1954/data/ -o fusioncatcher_1954/

mkdir -p fusioncatcher_1395
mkdir -p fusioncatcher_1395/data
ln -s ${ExeDir}/dataHCC/SRR2532336_1.fastq fusioncatcher_1395/data/
ln -s ${ExeDir}/dataHCC/SRR2532336_2.fastq fusioncatcher_1395/data/
fusioncatcher -d ${FusionCatcherDataDir} -i fusioncatcher_1395/data/ -o fusioncatcher_1395/

# JAFFA
mkdir -p jaffa_1954
mkdir -p jaffa_1954/data
ln -s ${ExeDir}/dataHCC/RNAHCC1954_1.fastq jaffa_1954/data
ln -s ${ExeDir}/dataHCC/RNAHCC1954_2.fastq jaffa_1954/data
cd jaffa_1954
${JaffaDir}/tools/bin/bpipe run ${JaffaDir}/JAFFA_assembly.groovy jaffa_1954/data/RNAHCC1954_*.fastq
cd ..

mkdir -p jaffa_1395
mkdir -p jaffa_1395/data
ln -s ${ExeDir}/dataHCC/SRR2532336_1.fastq jaffa_1395/data/
ln -s ${ExeDir}/dataHCC/SRR2532336_2.fastq jaffa_1395/data/
cd jaffa_1395
${JaffaDir}/tools/bin/bpipe run ${JaffaDir}/JAFFA_assembly.groovy jaffa_1395/data/SRR2532336_*.fastq
cd ..

# deFuse
mkdir -p defuse_1954
${DefuseDir}/scripts/defuse_run.pl -c ${DefuseDir}/scripts/config.txt -d ${DefuseDir}/RefData/ -1 ${ExeDir}/dataHCC/RNAHCC1954_1.fastq -2 ${ExeDir}/dataHCC/RNAHCC1954_2.fastq -o defuse_1954 -p 4

mkdir -p defuse_1395
${DefuseDir}/scripts/defuse_run.pl -c ${DefuseDir}/scripts/config.txt -d ${DefuseDir}/RefData/ -1 ${ExeDir}/dataHCC/SRR2532336_1.fastq -2 ${ExeDir}/dataHCC/SRR2532336_2.fastq -o defuse_1954 -p 4

# Chimerascan
mkdir -p chimerascan_1954
python ${ChimerascanDir}/bin/chimerascan_run.py -v -p 4 ${OutputDir}/Ensemble75/CSBowtieIndex ${ExeDir}/dataHCC/RNAHCC1954_1.fastq ${ExeDir}/dataHCC/RNAHCC1954_2.fastq chimerascan_1954/

mkdir chimerascan_1395
python ${ChimerascanDir}/bin/chimerascan_run.py -v -p 4 ${OutputDir}/Ensemble75/CSBowtieIndex ${ExeDir}/dataHCC/SRR2532336_1.fastq ${ExeDir}/dataHCC/SRR2532336_2.fastq chimerascan_1395/

# INTEGRATE
samtools merge $OutputDir/StarAlign_1954/MergedAlign.bam $OutputDir/StarAlign_1954/Aligned.sortedByCoord.out.bam $OutputDir/StarAlign_1954/Chimeric.out.bam
samtools sort $OutputDir/StarAlign_1954/MergedAlign.bam -o $OutputDir/StarAlign_1954/MergedAlign_sort.bam
samtools index $OutputDir/StarAlign_1954/MergedAlign_sort.bam
java -jar ${}/picard-jar FastqToSam F1=$OutputDir/StarAlign_1954/Unmapped.out.mate1 F2=$OutputDir/StarAlign_1954/Unmapped.out.mate2 O=$OutputDir/StarAlign_1954/Unmapped.bam SM=HCC1954 
mkdir -p integrate_1954
${IntegrateIntallDir}/bin/Integrate fusion ${GenomeFasta} ${ExeDir}/dataHCC/annot.ensembl.txt ${OutputDir}/Ensemble75/IntegrateIndex $OutputDir/StarAlign_1954//MergedAlign_sort.bam $OutputDir/StarAlign_1954/Unmapped.bam

samtools merge $OutputDir/StarAlign_1395/MergedAlign.bam $OutputDir/StarAlign_1395/Aligned.sortedByCoord.out.bam $OutputDir/StarAlign_1395/Chimeric.out.bam
samtools sort $OutputDir/StarAlign_1395/MergedAlign.bam -o $OutputDir/StarAlign_1395/MergedAlign_sort.bam
samtools index $OutputDir/StarAlign_1395/MergedAlign_sort.bam
java -jar ${}/picard-jar FastqToSam F1=$OutputDir/StarAlign_1395/Unmapped.out.mate1 F2=$OutputDir/StarAlign_1395/Unmapped.out.mate2 O=$OutputDir/StarAlign_1395/Unmapped.bam SM=HCC1395
mkdir -p integrate_1395
${IntegrateIntallDir}/bin/Integrate fusion ${GenomeFasta} ${ExeDir}/dataHCC/annot.ensembl.txt ${OutputDir}/Ensemble75/IntegrateIndex $OutputDir/StarAlign_1395//MergedAlign_sort.bam $OutputDir/StarAlign_1395/Unmapped.bam
