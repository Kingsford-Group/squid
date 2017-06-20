#!/bin/bash

args=("$@")
for (( i=0; i<${#args[@]}; )); do
	if [ ${args[$i]} = "-p" ]; then
		if [[ ${args[$(($i+1))]} = /* || ${args[$(($i+1))]} = ~* ]]; then
			ProjectDir=${args[$(($i+1))]}
		else
			ProjectDir=$(pwd)/${args[$(($i+1))]}
		fi
		((i+=2))
	fi
done
echo "Project directory: "$ProjectDir
echo "ILP Executable directory: "$ExecutableDir

# ILP Rearrangement
echo "Running ILP SV calling"
cd $ProjectDir/Alignments/
for (( i=3; i<10; i++)); do
	squid-b Star_rearranged/Aligned.sortedByCoord.out.bam -f ../WholeGenome/genome_rearranged.fa -c Star_rearranged/Chimeric.out.bam -s w$i -w $i
done
for (( i=3; i<10; i++)); do
	squid2 -b SpeedSeq/Aligned.bam -f ../WholeGenome/genome_rearranged.fa -s w$i -w $i
done

# SV metrics (ILP, Delly, Lumpy)
cd $ProjectDir
if [ ! -d $ProjectDir/SVcall ]; then
	echo "making directory "$ProjectDir/SVcall
	mkdir $ProjectDir/SVcall
fi
if [ ! -d $ProjectDir/SVcall/DELLY_BWA ]; then
	echo "making directory "$ProjectDir/SVcall/DELLY_BWA
	mkdir $ProjectDir/SVcall/DELLY_BWA
fi
if [ ! -d $ProjectDir/SVcall/DELLY_STAR ]; then
	echo "making directory "$ProjectDir/SVcall/DELLY_STAR
	mkdir $ProjectDir/SVcall/DELLY_STAR
fi
if [ ! -d $ProjectDir/SVcall/LUMPY ]; then
	echo "making directory "$ProjectDir/SVcall/LUMPY
	mkdir $ProjectDir/SVcall/LUMPY
fi
if [ ! -d $ProjectDir/SVcall/TransAbyss ]; then
	echo "making directory" $ProjectDir/SVcall/TransAbyss
	mkdir $ProjectDir/SVcall/TransAbyss
fi
if [ ! -d $ProjectDir/SVcall/Trinity ]; then
	echo "making directory" $ProjectDir/SVcall/Trinity
	mkdir $ProjectDir/SVcall/Trinity
fi
# Delly_STAR
delly call -t DEL -o SVcall/DELLY_STAR/delly_DEL.bcf -g WholeGenome/genome_rearranged.fa Alignments/Star_rearranged/Merged.bam
delly call -t DUP -o SVcall/DELLY_STAR/delly_DUP.bcf -g WholeGenome/genome_rearranged.fa Alignments/Star_rearranged/Merged.bam
delly call -t TRA -o SVcall/DELLY_STAR/delly_TRA.bcf -g WholeGenome/genome_rearranged.fa Alignments/Star_rearranged/Merged.bam
delly call -t INV -o SVcall/DELLY_STAR/delly_INV.bcf -g WholeGenome/genome_rearranged.fa Alignments/Star_rearranged/Merged.bam
delly call -t INS -o SVcall/DELLY_STAR/delly_INS.bcf -g WholeGenome/genome_rearranged.fa Alignments/Star_rearranged/Merged.bam
bcftools merge --force-samples -m id -O v -o SVcall/DELLY_STAR/delly_SV.vcf SVcall/DELLY_STAR/delly_*.bcf
python3 scripts/VerifySVpred.py 2 WholeGenome/SV_newpos.txt SVcall/DELLY_STAR/delly_SV.vcf SVcall/DELLY_STAR/hit.txt
# Delly_BWA
delly call -t DEL -o SVcall/DELLY_BWA/delly_DEL.bcf -g WholeGenome/genome_rearranged.fa Alignments/SpeedSeq/Aligned.bam
delly call -t DUP -o SVcall/DELLY_BWA/delly_DUP.bcf -g WholeGenome/genome_rearranged.fa Alignments/SpeedSeq/Aligned.bam
delly call -t TRA -o SVcall/DELLY_BWA/delly_TRA.bcf -g WholeGenome/genome_rearranged.fa Alignments/SpeedSeq/Aligned.bam
delly call -t INV -o SVcall/DELLY_BWA/delly_INV.bcf -g WholeGenome/genome_rearranged.fa Alignments/SpeedSeq/Aligned.bam
delly call -t INS -o SVcall/DELLY_BWA/delly_INS.bcf -g WholeGenome/genome_rearranged.fa Alignments/SpeedSeq/Aligned.bam
bcftools merge --force-samples -m id -O v -o SVcall/DELLY_BWA/delly_SV.vcf SVcall/DELLY_BWA/delly_*.bcf
python3 scripts/VerifySVpred.py 2 WholeGenome/SV_newpos.txt SVcall/DELLY_BWA/delly_SV.vcf SVcall/DELLY_BWA/hit.txt

# Lumpy
for (( i=3; i<10; i++)); do
	lumpyexpress -B Alignments/SpeedSeq/Aligned.bam -D Alignments/SpeedSeq/Aligned.discordants.bam -S Alignments/SpeedSeq/Aligned.splitters.bam -o SVcall/LUMPY/lumpy_thresh$i.vcf -m $i
	python3 scripts/VerifySVpred.py 2 WholeGenome/SV_newpos.txt SVcall/LUMPY/lumpy_thresh$i.vcf SVcall/LUMPY/hit_thresh$i.txt
done

# preparing for gmap
if [ ! -d $ProjectDir/WholeGenome/Chromosomes/ ]; then
	mkdir $ProjectDir/WholeGenome/Chromosomes/
fi
if [ ! -e $ProjectDir/WholeGenome/Chromosomes/chr1.fa ]; then
	python3 scripts/ReadGenome.py $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/WholeGenome/Chromosomes/
fi
if [ ! -d $ProjectDir/WholeGenome/Chromosomes/"genome"${ProjectDir##*/} ]; then
	gmap_build -D $ProjectDir/WholeGenome/Chromosomes/ -d "genome"${ProjectDir##*/} $ProjectDir/WholeGenome/Chromosomes/*.fa
fi

# transabyss + mummer3
transabyss -k 25 --pe $ProjectDir/Reads/RNA1.fq.gz $ProjectDir/Reads/RNA2.fq.gz --outdir $ProjectDir/SVcall/TransAbyss/k25 --name k25 --threads 4 --island 0
transabyss -k 32 --pe $ProjectDir/Reads/RNA1.fq.gz $ProjectDir/Reads/RNA2.fq.gz --outdir $ProjectDir/SVcall/TransAbyss/k32 --name k32 --threads 4 --island 0
transabyss-merge --mink 25 --maxk 32 --prefixes k25. k32. --out $ProjectDir/SVcall/TransAbyss/mergedassembly.fa $ProjectDir/SVcall/TransAbyss/k25/k25-final.fa $ProjectDir/SVcall/TransAbyss/k32/k32-final.fa
nucmer -p $ProjectDir/SVcall/TransAbyss/nucm $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/SVcall/TransAbyss/mergedassembly.fa
show-coords $ProjectDir/SVcall/TransAbyss/nucm.delta > $ProjectDir/SVcall/TransAbyss/nucm_text.txt
./bin/NucmerSV2 $ProjectDir/SVcall/TransAbyss/nucm_text.txt $ProjectDir/SVcall/TransAbyss/nucm_res.bedpe
python3 scripts/VerifySVpred.py 5 $ProjectDir/WholeGenome/SV_newpos.txt $ProjectDir/SVcall/TransAbyss/nucm_res.bedpe $ProjectDir/SVcall/TransAbyss/nucm_res_hit.bedpe

# transabyss + gmap
gmap -D $ProjectDir/WholeGenome/Chromosomes/ -d "genome"${ProjectDir##*/} $ProjectDir/SVcall/TransAbyss/mergedassembly.fa > $ProjectDir/SVcall/TransAbyss/gmap.out
./bin/GmapSV $ProjectDir/SVcall/TransAbyss/gmap.out $ProjectDir/SVcall/TransAbyss/gmap_res.bedpe
python3 scripts/VerifySVpred.py 5 $ProjectDir/WholeGenome/SV_newpos.txt $ProjectDir/SVcall/TransAbyss/gmap_res.bedpe $ProjectDir/SVcall/TransAbyss/gmap_res_hit.bedpe

# trinity + mummer3
Trinity --seqType fq --max_memory 50G --left $ProjectDir/Reads/RNA1.fq.gz --right $ProjectDir/Reads/RNA2.fq.gz --CPU 6 --output $ProjectDir/SVcall/Trinity
nucmer -p $ProjectDir/SVcall/Trinity/nucm $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/SVcall/Trinity/Trinity.fasta
show-coords $ProjectDir/SVcall/Trinity/nucm.delta > $ProjectDir/SVcall/Trinity/nucm_text.txt
./bin/NucmerSV2 $ProjectDir/SVcall/Trinity/nucm_text.txt $ProjectDir/SVcall/Trinity/nucm_res.bedpe
python3 scripts/VerifySVpred.py 5 $ProjectDir/WholeGenome/SV_newpos.txt $ProjectDir/SVcall/Trinity/nucm_res.bedpe $ProjectDir/SVcall/Trinity/nucm_res_hit.bedpe

# trinity + gmap
gmap -D $ProjectDir/WholeGenome/Chromosomes/ -d "genome"${ProjectDir##*/} $ProjectDir/SVcall/Trinity/Trinity.fasta > $ProjectDir/SVcall/Trinity/gmap.out
./bin/GmapSV $ProjectDir/SVcall/Trinity/gmap.out $ProjectDir/SVcall/Trinity/gmap_res.bedpe
python3 scripts/VerifySVpred.py 5 $ProjectDir/WholeGenome/SV_newpos.txt $ProjectDir/SVcall/Trinity/gmap_res.bedpe $ProjectDir/SVcall/Trinity/gmap_res_hit.bedpe