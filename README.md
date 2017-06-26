# OVERVIEW
SQUID is designed to detect transcriptomic structural variations from RNA-seq alignment.

# INSTALLATION
SQUID requires boost, gurobi, bamtools. After installing these packages, change the INCLUDES and LDFLAGS path in Makefile. And them make.

Gurobi doesn't support gcc5 currently. If you are using gcc5, you need to add -D_GLIBCXX_USE_CXX11_ABI=0 to CXXFLAGS in Makefile

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

