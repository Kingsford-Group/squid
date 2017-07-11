#!/bin/bash

DataFolder=./simulationdata
threads=8

cp scripts/SVcalling.sh bin/

for ((i=2; i<9; i+=3)); do
	for ((j=1; j<5; j++)); do
		mkdir -p $DataFolder/SV${i}00_${j}/Alignments

		# Aligning reads with STAR
		if [ ! -d $DataFolder/SV${i}00_${j}/WholeGenome/STAR_genome_rearranged ]; then
			STAR --runThreadN ${threads} --runMode genomeGenerate --genomeDir $DataFolder/SV${i}00_${j}/WholeGenome/STAR_genome_rearranged --genomeFastaFiles $DataFolder/SV${i}00_${j}/WholeGenome/genome_rearranged.fa
		fi
		mkdir -p $DataFolder/SV${i}00_${j}/Alignments/Star_rearranged
		STAR --runThreadN ${threads} --genomeDir $DataFolder/SV${i}00_${j}/WholeGenome/STAR_genome_rearranged/ --readFilesIn $DataFolder/SV${i}00_${j}/Reads/RNA1.fq.gz $DataFolder/SV${i}00_${j}/Reads/RNA2.fq.gz --readFilesCommand gunzip -c --outFileNamePrefix $DataFolder/SV${i}00_${j}/Alignments/Star_rearranged/ --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --chimSegmentMin 20 --outSAMstrandField intronMotif --limitBAMsortRAM 21943468974

		# Aligning reads with SpeedSeq
		mkdir -p $DataFolder/SV${i}00_${j}/Alignments/SpeedSeq
		speedseq align -R "@RG\tID:id\tSM:samplename\tLB:lib" -t ${threads} -o $DataFolder/SV${i}00_${j}/Alignments/SpeedSeq/Aligned $DataFolder/SV${i}00_${j}/WholeGenome/genome_rearranged.fa $DataFolder/SV${i}00_${j}/Reads/RNA1.fq.gz $DataFolder/SV${i}00_${j}/Reads/RNA2.fq.gz

		# Call TSVs with different methods
		./bin/SVcalling.sh -p ${OutDir}/SVRNA${i}00_${j}
	done
done