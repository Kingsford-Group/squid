#!/bin/bash

OutDir=./
GenomeDir=./genome
AnnotationFile=./annotation/annot.gtf
BamtoolsDir=./bamtools

mkdir -p OutDir

# build executables
cp scripts/SimSVGenome.sh bin/
cp scripts/SVcalling.sh bin/
cp scripts/extractSTAR.sh bin/
chmod a+x ./bin/SimSVGenome.sh
chmod a+x ./bin/SVcalling.sh
chmod a+x ./bin/extractSTAR.sh
g++ -std=c++11 -o bin/SV2newpos scripts/SV2newpos.cpp scripts/GtfTrans.cpp scripts/SV.cpp scripts/TRA.cpp scripts/SimpleSV.cpp -g
g++ -std=c++11 -o bin/GmapSV scripts/GmapSV.cpp -g
g++ -std=c++11 -o bin/NucmerSV2 scripts/NucmerSV2.cpp -g
g++ -std=c++11 -o bin/mergesplit scripts/mergesplit.cpp -I ${BamtoolsDir}/include/ -L ${BamtoolsDir}/lib/ -lbamtools -lz -g -Wl,-rpath,${BamtoolsDir}/lib

for ((i=200}; i<900; i+=300)); do
	./bin/SimSVGenome.sh -g ${GenomeDir} -p ${OutDir}/SVRNA${i}_1 -a ${AnnotationFile} -seq 21 -inv ${i} -ins ${i} -del ${i} -dup ${i} -tra 2 -RNA
	./bin/SimSVGenome.sh -g ${GenomeDir} -p ${OutDir}/SVRNA${i}_2 -a ${AnnotationFile} -seq 25 -inv ${i} -ins ${i} -del ${i} -dup ${i} -tra 2 -RNA
	./bin/SimSVGenome.sh -g ${GenomeDir} -p ${OutDir}/SVRNA${i}_3 -a ${AnnotationFile} -seq 28 -inv ${i} -ins ${i} -del ${i} -dup ${i} -tra 2 -RNA
	./bin/SimSVGenome.sh -g ${GenomeDir} -p ${OutDir}/SVRNA${i}_4 -a ${AnnotationFile} -seq 30 -inv ${i} -ins ${i} -del ${i} -dup ${i} -tra 2 -RNA
done

for ((i=200; i<900; i+=300)); do
	./bin/SVcalling.sh -p ${OutDir}/SVRNA${i}_1
	./bin/SVcalling.sh -p ${OutDir}/SVRNA${i}_2
	./bin/SVcalling.sh -p ${OutDir}/SVRNA${i}_3
	./bin/SVcalling.sh -p ${OutDir}/SVRNA${i}_4
done
