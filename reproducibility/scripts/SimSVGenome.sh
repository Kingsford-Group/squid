#!/bin/bash

args=("$@")
for (( i=0; i<${#args[@]}; )); do
    if [ ${args[$i]} = "-g" ]; then
        if [[ ${args[$(($i+1))]} = /* || ${args[$(($i+1))]} = ~* ]]; then
            OldGenomeDir=${args[$(($i+1))]}
        else
            OldGenomeDir=$(pwd)/${args[$(($i+1))]}
        fi
        ((i+=2))
    elif [ ${args[$i]} = "-p" ]; then
        if [[ ${args[$(($i+1))]} = /* || ${args[$(($i+1))]} = ~* ]]; then
            ProjectDir=${args[$(($i+1))]}
        else
            ProjectDir=$(pwd)/${args[$(($i+1))]}
        fi
        ((i+=2))
    elif [ ${args[$i]} = "-a" ]; then
        if [[ ${args[$(($i+1))]} = /* || ${args[$(($i+1))]} = ~* ]]; then
            AnnotationFile=${args[$(($i+1))]}
        else
            AnnotationFile=$(pwd)/${args[$(($i+1))]}
        fi
        ((i+=2))
    elif [ ${args[$i]} = "-seq" ]; then
        NumSEQ=${args[$(($i+1))]}
        ((i+=2))
    elif [ ${args[$i]} = "-inv" ]; then
        NumINV=${args[$(($i+1))]}
        ((i+=2))
    elif [ ${args[$i]} = "-ins" ]; then
        NumINS=${args[$(($i+1))]}
        ((i+=2))
    elif [ ${args[$i]} = "-del" ]; then
        NumDEL=${args[$(($i+1))]}
        ((i+=2))
    elif [ ${args[$i]} = "-dup" ]; then
        NumDUP=${args[$(($i+1))]}
        ((i+=2))
    elif [ ${args[$i]} = "-tra" ]; then
        NumTRA=${args[$(($i+1))]}
        ((i+=2))
    elif [ ${args[$i]} = "-RNA" ]; then
        RNA=true
        ((i++))
    elif [ ${args[$i]} = "-WGS" ]; then
        RNA=false
        ((i++))
    fi
done

if [ ! -d $ProjectDir ]; then
    echo "making directory "$ProjectDir
    mkdir $ProjectDir
fi
# Param log file
for (( i=0; i<${#args[@]}; i++ )); do
    echo -n ${args[$i]}" " >> $ProjectDir/Params.log
done

if [ ! -d $ProjectDir/WholeGenome ]; then
    echo "making directory "$ProjectDir/WholeGenome
    mkdir $ProjectDir/WholeGenome
fi
if [ ! -d $ProjectDir/Reads ]; then
    echo "making directory "$ProjectDir/Reads
    mkdir $ProjectDir/Reads
fi
if [ ! -d $ProjectDir/Alignments ]; then
    echo "making directory "$ProjectDir/Alignments
    mkdir $ProjectDir/Alignments
fi

cd $ProjectDir/WholeGenome
# generate merged old genome
for file in $( ls $OldGenomeDir/ ); do
    cat $OldGenomeDir/$file >> genome_original.fa
done
# RSVSim to add SV to genome
echo "library('RSVSim')" > simSV.R
count=0
for file in $( ls $OldGenomeDir/ ); do
    echo "chr"$count"=readDNAStringSet(\""$OldGenomeDir/$file"\")" >> simSV.R
    ((count++))
done
str="dnalist=list("
for (( i=0;i<$count-1;i++)); do
    str=$str"chr"$i","
done
str=$str"chr"$(($count-1))")"
echo $str >> simSV.R
echo "genome=do.call(c,dnalist)" >> simSV.R
echo >> simSV.R
echo "invSizes=estimateSVSizes(n="$NumINV", minSize=500, maxSize=20000, default=\"inversions\", hist=FALSE)" >> simSV.R
echo "insSizes=estimateSVSizes(n="$NumINS", minSize=500, maxSize=20000, default=\"insertions\", hist=FALSE)" >> simSV.R
echo "delSizes=estimateSVSizes(n="$NumDEL", minSize=500, maxSize=20000, default=\"deletions\", hist=FALSE)" >> simSV.R
echo "dupSizes=estimateSVSizes(n="$NumDUP", minSize=500, maxSize=20000, default=\"tandemDuplications\",hist=FALSE)" >> simSV.R
echo "sim=simulateSV(genome=genome,ins="$NumINS",sizeIns=insSizes,invs="$NumINV",sizeInvs=invSizes,dels="$NumDEL",sizeDels=delSizes,dups="$NumDUP",sizeDups=dupSizes,maxDups=2,trans="$NumTRA")" >> simSV.R
Rscript simSV.R
mv genome_rearranged.fasta genome_rearranged.fa
if $RNA ; then
    mkdir STAR_genome_rearranged
    STAR --runThreadN 8 --runMode genomeGenerate --genomeDir STAR_genome_rearranged --genomeFastaFiles genome_rearranged.fa
    mv Log.out STAR_genome_rearranged/Log.out
fi

# Simulate sequencing Reads
cd $ProjectDir/Reads
if $RNA ; then
    if [ ! -d $ProjectDir/Annotation ]; then
        echo "making directory "$ProjectDir/Annotation
        mkdir $ProjectDir/Annotation
    fi
    ln -s $AnnotationFile $ProjectDir/Annotation/
    echo -e "REF_FILE_NAME\t"$AnnotationFile > fluxsim.par
    echo -e "GEN_DIR\t"$OldGenomeDir >> fluxsim.par
    echo -e "READ_LENGTH\t76\n" >> fluxsim.par
    echo -e "NB_MOLECULES\t40000000\nREAD_NUMBER\t"$NumSEQ"000000\nPAIRED_END\tYES\nFILTERING\tTRUE\n" >> fluxsim.par
    echo -e "# use default 76-bp error model\nERR_FILE\t76\n" >> fluxsim.par
    echo -e "# create a fastq file\nFASTA\tYES\nUNIQUE_IDS\tYES" >> fluxsim.par
    echo -e "# assign tmp directory\nTMP_DIR\ttmp" >> fluxsim.par
    mkdir tmp
    flux-simulator -p fluxsim.par
    HasQual=false;
    if [ $(head -3 fluxsim.fastq | tail -1) = "+" ]; then
        HasQual=true
    fi
    if $HasQual ; then
        awk '{if(NR%8==1) print substr($1,0,length($1)-2); else if(NR%8<5 && NR%8>0) print $0;}' fluxsim.fastq > RNA1.fastq
        awk '{if(NR%8==5) print substr($1,0,length($1)-2); else if(NR%8>4 || NR%8<1) print $0;}' fluxsim.fastq > RNA2.fastq
    else
        awk '{if(NR%4==1) print substr($1,0,length($1)-2); else if(NR%4<3 && NR%4>0) print $0;}' fluxsim.fastq > RNA1.fastq
        awk '{if(NR%4==3) print substr($1,0,length($1)-2); else if(NR%4>2 || NR%4<1) print $0;}' fluxsim.fastq > RNA2.fastq
    fi
    gzip RNA1.fastq
    gzip RNA2.fastq
    rm -fr tmp
    rm -fr fluxsim.fastq
else
    art_illumina -i $ProjectDir/WholeGenome/genome_original.fa -p -l 100 -f 10 -m 200 -s 10 -o WGS
    aln2bed.pl WGS1.bed WGS1.aln
    aln2bed.pl WGS2.bed WGS2.aln
    while read -r a && read -r b <&2; do
        echo -e "$a\n$b" >> WGS.bed
    done <WGS1.bed 2<WGS2.bed
    rm -fr WGS1.aln WGS2.aln WGS1.bed WGS2.bed
fi

# Alignment with STAR
cd $ProjectDir/Alignments
if $RNA ; then
    mkdir Star_rearranged
    STAR --runThreadN 8 --genomeDir $ProjectDir/WholeGenome/STAR_genome_rearranged/ --readFilesIn $ProjectDir/Reads/RNA1.fastq.gz $ProjectDir/Reads/RNA2.fastq.gz --readFilesCommand gunzip -c --outFileNamePrefix Star_rearranged/ --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --chimSegmentMin 20 --outSAMstrandField intronMotif --limitBAMsortRAM 21943468974
    samtools view -Shb Star_rearranged/Chimeric.out.sam > Star_rearranged/Chimeric.out.bam
    # merge star concordant and chimeric alignment files
    ./bin/MergeSTAR/extractSTAR.sh -b Star_rearranged/Aligned.sortedByCoord.out.bam -c Star_rearranged/Chimeric.out.bam
fi

# Alignment with speedseq
mkdir SpeedSeq
if $RNA ; then
    speedseq align -R "@RG\tID:id\tSM:samplename\tLB:lib" -t 8 -o SpeedSeq/Aligned $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/Reads/RNA1.fastq.gz $ProjectDir/Reads/RNA2.fastq.gz
else
    speedseq align -R "@RG\tID:id\tSM:samplename\tLB:lib" -t 8 -o SpeedSeq/Aligned $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/Reads/WGS1.fq $ProjectDir/Reads/WGS2.fq
fi

# Calculating read-covered SV and output in new coordinate
if $RNA ; then
    ./bin/SV2newpos $ProjectDir/WholeGenome/genome_original.fa $ProjectDir/WholeGenome/ $ProjectDir/WholeGenome/SV_newpos.txt $ProjectDir/Reads/fluxsim.bed old $AnnotationFile
fi
