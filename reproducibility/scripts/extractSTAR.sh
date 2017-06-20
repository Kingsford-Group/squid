#!/bin/bash

args=("$@")

BamFile=""
ChimBamFile=""

ID="id"
LB="lib"
SM="sample"
PL="illumina"
PU=10

for (( i=0; i<${#args[@]}; )); do
    if [ ${args[$i]} = "-b" ]; then
        BamFile=${args[$(($i+1))]}
        ((i+=2))
    elif [ ${args[$i]} = "-c" ]; then
        ChimBamFile=${args[$(($i+1))]}
        ((i+=2))
    elif [ ${args[$i]} = "--ID" ]; then
        ID=${args[$(($i+1))]}
        ((i+=2))
    elif [ ${args[$i]} = "--LB" ]; then
        LB=${args[$(($i+1))]}
        ((i+=2))
    elif [ ${args[$i]} = "--SM" ]; then
        SM=${args[$(($i+1))]}
        ((i+=2))
    elif [ ${args[$i]} = "--PL" ]; then
        PL=${args[$(($i+1))]}
        ((i+=2))
    elif [ ${args[$i]} = "--PU" ]; then
        PU=${args[$(($i+1))]}
        ((i+=2))
    else
        echo "Unrecognized argument"
        exit
    fi
done

if [ ${#BamFile} == 0 ] || [ ${#ChimBamFile} == 0 ]; then
    echo "Invalid bamfile or chimeric bamfile"
    exit
fi

BamDir=${BamFile%/*}

echo "BamFile: "$BamFile
echo "ChimBamFile: "$ChimBamFile
echo "BamDir: "$BamDir

java -jar ~/cong/Software/picard-tools-1.141/picard.jar AddOrReplaceReadGroups I=$BamFile O=$BamDir/RG.bam RGID=$ID RGLB=$LB RGPL=$PL RGPU=$PU RGSM=$SM
java -jar ~/cong/Software/picard-tools-1.141/picard.jar AddOrReplaceReadGroups I=$ChimBamFile O=$BamDir/RGchim.bam RGID=$ID RGLB=$LB RGPL=$PL RGPU=$PU RGSM=$SM

~/cong/Code/STARforLUMPY/mergesplit -b $BamDir/RG.bam -c $BamDir/RGchim.bam -p Merged.unsorted
#samtools view -b -F 1294 $BamDir/Merged.unsorted.bam > $BamDir/Merged.unsorted.discordant.bam

samtools sort $BamDir/Merged.unsorted.bam -o $BamDir/Merged.bam
#samtools sort $BamDir/Merged.unsorted.discordant.bam -o $BamDir/Merged.discordant.bam
#samtools sort $BamDir/Merged.unsorted.splitters.bam -o $BamDir/Merged.splitters.bam

rm -fr RG*
rm -fr Merged.unsorted.
