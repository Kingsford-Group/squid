#!/bin/bash

source ~/cong/Software/PyVirtual/bin/activate

for ((i=2; i<9; i+=3)); do
	for ((j=1; j<5; j++)); do
		if [ ! -d "SVRNA"$i"00_"$j/TransAbyss ]; then
			mkdir "SVRNA"$i"00_"$j/TransAbyss
		fi
		# running transabyss and mummer3
		transabyss -k 25 --pe "SVRNA"$i"00_"$j/Reads/RNA1.fq.gz "SVRNA"$i"00_"$j/Reads/RNA2.fq.gz --outdir "SVRNA"$i"00_"$j/TransAbyss/k25 --name k25 --threads 4 --island 0
		transabyss -k 32 --pe "SVRNA"$i"00_"$j/Reads/RNA1.fq.gz "SVRNA"$i"00_"$j/Reads/RNA2.fq.gz --outdir "SVRNA"$i"00_"$j/TransAbyss/k32 --name k32 --threads 4 --island 0
		transabyss-merge --mink 25 --maxk 32 --prefixes k25. k32. --out "SVRNA"$i"00_"$j/TransAbyss/mergedassembly.fa "SVRNA"$i"00_"$j/TransAbyss/k25/k25-final.fa "SVRNA"$i"00_"$j/TransAbyss/k32/k32-final.fa
		nucmer -p "SVRNA"$i"00_"$j/TransAbyss/nucm "SVRNA"$i"00_"$j/WholeGenome/genome_rearranged.fa "SVRNA"$i"00_"$j/TransAbyss/mergedassembly.fa
		show-coords "SVRNA"$i"00_"$j/TransAbyss/nucm.delta > "SVRNA"$i"00_"$j/TransAbyss/nucm_text.txt
		/home/congm1/cong/Code/Virgorrangement/NucmerSV2 ~/congm/"SVRNA"$i"00_"$j/TransAbyss/nucm_text.txt ~/congm/"SVRNA"$i"00_"$j/TransAbyss/nucm_res.bedpe
		/home/congm1/cong/Code/Virgorrangement/NucmerSV2 ~/congm/"SVRNA"$i"00_"$j/Trinity/nucm_text.txt ~/congm/"SVRNA"$i"00_"$j/Trinity/nucm_res.bedpe

		# running trinity and mummer3
		Trinity --seqType fq --max_memory 50G --left "SVRNA"$i"00_"$j/Reads/RNA1.fq.gz --right "SVRNA"$i"00_"$j/Reads/RNA2.fq.gz --CPU 6 --output "SVRNA"$i"00_"$j/Trinity
		nucmer -p "SVRNA"$i"00_"$j/Trinity/nucm "SVRNA"$i"00_"$j/WholeGenome/genome_rearranged.fa "SVRNA"$i"00_"$j/Trinity/Trinity.fasta
		show-coords "SVRNA"$i"00_"$j/Trinity/nucm.delta > "SVRNA"$i"00_"$j/Trinity/nucm_text.txt
		/home/congm1/cong/Code/Virgorrangement/VerifyNucmerSV ~/congm/"SVRNA"$i"00_"$j/Alignments/SpeedSeq/Aligned.bam ~/congm/"SVRNA"$i"00_"$j/WholeGenome/SV_newpos.txt ~/congm/"SVRNA"$i"00_"$j/TransAbyss/nucm_res.bedpe
		/home/congm1/cong/Code/Virgorrangement/VerifyNucmerSV ~/congm/"SVRNA"$i"00_"$j/Alignments/SpeedSeq/Aligned.bam ~/congm/"SVRNA"$i"00_"$j/WholeGenome/SV_newpos.txt ~/congm/"SVRNA"$i"00_"$j/Trinity/nucm_res.bedpe

		# preparing GMAP
		if [ ! -d "SVRNA"$i"00_"$j/WholeGenome/Chromosomes/ ]; then
			mkdir "SVRNA"$i"00_"$j/WholeGenome/Chromosomes/
		fi
		if [ ! -e "SVRNA"$i"00_"$j/WholeGenome/Chromosomes/chr1.fa ]; then
			python3 ReadGenome.py "SVRNA"$i"00_"$j/WholeGenome/genome_rearranged.fa "SVRNA"$i"00_"$j/WholeGenome/Chromosomes/
		fi
		if [ ! -d "SVRNA"$i"00_"$j/WholeGenome/Chromosomes/"genome"$i"00_"$j ]; then
			gmap_build -D "SVRNA"$i"00_"$j/WholeGenome/Chromosomes/ -d "genome"$i"00_"$j "SVRNA"$i"00_"$j/WholeGenome/Chromosomes/*.fa
		fi
		# running gmap on transabyss and trinity transcripts
		gmap -D "SVRNA"$i"00_"$j/WholeGenome/Chromosomes/ -d "genome"$i"00_"$j "SVRNA"$i"00_"$j/TransAbyss/mergedassembly.fa > "SVRNA"$i"00_"$j/TransAbyss/gmap.out
		gmap -D "SVRNA"$i"00_"$j/WholeGenome/Chromosomes/ -d "genome"$i"00_"$j "SVRNA"$i"00_"$j/Trinity/Trinity.fasta > "SVRNA"$i"00_"$j/Trinity/gmap.out
		/home/congm1/cong/Code/Virgorrangement/GmapSV ~/congm/"SVRNA"$i"00_"$j/TransAbyss/gmap.out ~/congm/"SVRNA"$i"00_"$j/TransAbyss/gmap_res.bedpe
		/home/congm1/cong/Code/Virgorrangement/VerifyNucmerSV ~/congm/"SVRNA"$i"00_"$j/Alignments/SpeedSeq/Aligned.bam ~/congm/"SVRNA"$i"00_"$j/WholeGenome/SV_newpos.txt ~/congm/"SVRNA"$i"00_"$j/TransAbyss/gmap_res.bedpe
		/home/congm1/cong/Code/Virgorrangement/GmapSV ~/congm/"SVRNA"$i"00_"$j/Trinity/gmap.out ~/congm/"SVRNA"$i"00_"$j/Trinity/gmap_res.bedpe
		/home/congm1/cong/Code/Virgorrangement/VerifyNucmerSV ~/congm/"SVRNA"$i"00_"$j/Alignments/SpeedSeq/Aligned.bam ~/congm/"SVRNA"$i"00_"$j/WholeGenome/SV_newpos.txt ~/congm/"SVRNA"$i"00_"$j/Trinity/gmap_res.bedpe

	done
done

deactivate
