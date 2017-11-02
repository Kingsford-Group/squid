# OVERVIEW
SQUID is designed to detect transcriptomic structural variations from RNA-seq alignment.

To reproduce the result of applying SQUID on simulation data and previously studied cell lines, follow the instructions from [squidtest](https://github.com/Kingsford-Group/squidtest).

You can find a describtion of our methods in our [preprint](https://www.biorxiv.org/content/early/2017/09/06/162776) manuscript.

# INSTALLING PRE-COMPILED BINARIES
If you want to skip the installation process, you can find the binary release [here](https://github.com/Kingsford-Group/squid/releases)!

# BUILDING FROM SOURCE

You only need to build from source if either the pre-built binaries (see above) don't work on your system or you want to make a change to the SQUID code.

Compiling SQUID requires Boost, GLPK, BamTools. A step by step installation construction can be found [here for linux](doc/Installation_linux.md), and [here for mac](doc/Installation_mac.md).

On Mac, you need to additionly run the following command to dynamicly linking dependent libraries:
```
export DYLD_LIBRARY_PATH=<bamtools_folder>/lib
export DYLD_LIBRARY_PATH=<glpk_folder>/lib
```

# USAGE
```
squid [options] -b <Input_BAM> -o <Output_Prefix>
```
SQUID supports the following options:

 Parameters | Default value | Data type | Description 
 ---        | :---:         | :---:     | ---         
 -c         |               | string    |             
 -f         |               | string    |             
 -pt        |  0            | bool      | Phred type: 0 for Phred33, 1 for Phred64 
 -pl        |  10           | int       | Maximum Length of continuous low Phred score to filter alignment 
 -pm        |  4            | int       | Threshold to count as low Phred score 
 -mq        |  1            | int       | Minimum mapping quality 
 -dp        | 50000         | int       | Maximum paired-end aligning distance to be count as concordant alignment 
 -di        | 20            | int       | Maximum distance of segment indexes to be count as read-through 
 -w         |  5            | int       | Minimum edge weight 
 -G         |  0            | bool      | Whether or not output graph file (0 for not outputing, 1 for outputing) 
 -CO        |  0            | bool      | Whether or not output ordering of connected components (0 for not outputing, 1 for outputing) 
 -TO        |  0            | bool      | Whether or not output ordering of all segments (0 for not outputing, 1 for outputing) 
 -RG        |  0            | bool      | Whether or not output rearranged genome sequence (0 for not outputing, 1 for outputing) 

# OUTPUT SPECIFICATION
+ <Output_Prefix>_sv.txt: a list of predicted TSV in bedpe format. This is the main output of SQUID. All positions in the file are 0-based. Each columns represents:
	- chr1: chromosome name of the first breakpoint.
	- start1: starting position of the segment of the first breakpoint, or the predicted breakpoint position if strand1 is "-".
	- end1: ending position of the segment of the first breakpoint, or the predicted breakpoint position if strand1 is "+".
	- chr2: chromosome name of the second breakpoint.
	- start2: starting position of the segment of the second breakpoint, or the predicted breakpoint position if strand2 is "-".
	- end2: ending position of the segment of the second breakpoint, or the predicted breakpoint position if strand2 is "+".
	- name: TSV is not named yet, this column shows with dot.
	- score: Edge weight, or weighted number of reads supporting this TSV.
	- strand1: strand of the first segment in TSV.
	- strand2: strand of the second segment in TSV.
	example record:
	> 17  38135881  38136308  17  38137195  38137773  .  5480  +  + 
	>
	This means the right end (position 38136308) of segment 38135881-38137195 on chr17 is connected to the right end (position 38137773) of segment 38137195-38137773 also on chr17. And the score of this TSV is 5480.
	> 5  176370330  176370489  8  128043988  128044089  .  328  -  + 
	>
	This means the left end (position 176370330) of segment 176370330-176370489 on chr5 is connected with the right end (position 128044089) of segment 128043988-128044089 on chr8.
+ <Output_Prefix>_graph.txt: genome segment graph, will be output only if -G is set to 1. It has two types of records, nodes (or segment) and edges.
	- node: for each node, the following information are included. ID, start position, end position, label for connected component.
	- edge: for each edge, the following information are included. ID, node id for the first segment, strand of the first segment, node id for the second segment, strand for the second segment, edge weight.
+ <Output_Prefix>_component_pri.txt: ordering of each connected component by ILP, will be output only if -CO is set to 1. In this file, ordering of each connected component will be output into a line. Each segment is represented by its id, with a "-" in the front if the ordering suggests the segment should be using its reverse strand.
+ <Output_Prefix>_component.txt: ordering of the entire genome segments, will be output only is -TO is set to 1. In this file, each newly generated chromosome will be output into one line. Other spefications are the same as <Output_Prefix>_component_pri.txt.
+ <Output_Prefix>_genome.fa: sequence of rearranged genome sequence, where each chromosome corresponds to one line in <Output_Prefix>_component.txt.

# EXAMPLE WORKFLOW
Suppose you have the alignment BAM file, and chimeric BAM file generated by STAR (https://github.com/alexdobin/STAR),  run SQUID with:
```
squid -b alignment.bam -c chimeric.bam -o squidout
```
Or a combined BAM file of both concordant and discordant alignments generated by BWA (http://bio-bwa.sourceforge.net/) or SpeedSeq (https://github.com/hall-lab/speedseq), run SQUID with
```
squid --bwa -b combined_alignment.bam -o squidout
```

An example can be run be downloading the sample data (sampledata.tgz) from (https://cmu.box.com/s/e9u6alp73rfdhfve2a51p6v391vweodq) into example folder, and decompress it with
```
tar -xzvf sampledata.tgz
```
Run SQUID command in example/SQUIDcommand.sh. Or if you want to test the workflow of STAR and SQUID, make sure STAR is in your path, and run bash script example/STARnSQUIDcommand.sh.
```
cd example
./SQUIDcommand.sh
./STARnSQUIDcommand.sh
```
