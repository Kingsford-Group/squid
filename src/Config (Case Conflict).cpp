#include "config.h"

using namespace std;

// parameter from BAM
uint16_t ReadLen;

// parameters from user input
	// read BAM file
bool Phred_Type = 0; //0 for phred33, 1 for phred64
uint16_t Max_LowPhred_Len = 10;
uint8_t Min_Phred = 4;
uint16_t Min_MapQual = 1;
	// building genome segment graph
int Concord_Dist_Pos = 50000;
int Concord_Dist_Idx = 20;
int Min_Edge_Weight = 5;
	// input and output
string Input_BAM = "";
string Input_Chim_BAM = "";
string Input_FASTA = "";
string Output_Prefix = "";
bool Print_Graph = false;;
bool Print_Components_Ordering = false;
bool Print_Total_Ordering = false;
bool Print_Rearranged_Genome = false;

int print_help(){
	printf("Usage: squid [options] -b <Input_BAM> -o <Output_Prefix>\n");
	printf("Options:\n");

	printf("\tExtra input ptions:\n");
	printf("\t-c\tstring\tChimeric BAM alignment (If using STAR aligner)\n");
	printf("\t-f\tstring\tGenome FASTA file\n");
	
	printf("\tParsing alignment options:\n");
	printf("\t-pt\tbool\tPhred type: 0 for Phred33, 1 for Phred64 (0)\n");
	printf("\t-pl\tint\tMaximum Length of low Phred score to filter alignment (10)\n");
	printf("\t-pm\tint\tLow Phred score threshold (4)\n");
	printf("\t-mq\tint\tMapping quality to filter alignment (1)\n");

	printf("\tConstructing graph options:\n");
	printf("\t-dp\tint\tMaximum distance of aligning positions for concordant alignment (50000)\n");
	printf("\t-di\tint\tMaximum distance of segment indexes for concordant alignment (20)\n");

	printf("\tOutput options:\n");
	printf("\t-G\tbool\tOutput gragh file (0)\n");
	printf("\t-CO\tbool\tOutput ordering of connected components (0)\n");
	printf("\t-TO\tbool\tOutput ordering of all segments (0)\n");
	printf("\t-RG\tbool\tOutput rearranged genome sequence (0)\n");
};

bool parse_arguments(int argc, char* argv[]){
	bool success=true;
	bool specify_mq=false;
	for(int i=1; i<argc; i++){
		if(string(argv[i])=="-b" && i<argc-1){
			Input_BAM=string(argv[i+1]);
		}
		if(string(argv[i])=="-o" && i<argc-1){
			Output_Prefix=string(argv[i+1]);
		}
		if(string(argv[i])=="-c" && i<argc-1){
			Input_Chim_BAM=string(argv[i+1]);
		}
		if(string(argv[i])=="-f" && i<argc-1){
			Input_FASTA=string(argv[i+1]);
		}
		if(string(argv[i])=="-pt" && i<argc-1){
			if(string(argv[i+1])=="0")
				Phred_Type=0;
			else if(string(argv[i+1])=="1")
				Phred_Type=1;
			else
				success=false;
		}
		if(string(argv[i])=="-pl" && i<argc-1){
			try{
				Max_LowPhred_Len=(uint16_t)atoi(argv[i+1]);
			}
			catch(exception& e){
				success=false;
			}
		}
		if(string(argv[i])=="-pm" && i<argc-1){
			try{
				Min_Phred=(uint8_t)atoi(argv[i+1]);
			}
			catch(exception& e){
				success=false;
			}
		}
		if(string(argv[i])=="-mq" && i<argc-1){
			try{
				Min_MapQual=(uint16_t)atoi(argv[i+1]);
				specify_mq=true;
			}
			catch(exception& e){
				success=false;
			}
		}
		if(string(argv[i])=="-dp" && i<argc-1){
			try{
				Concord_Dist_Pos=atoi(argv[i+1]);
			}
			catch(exception& e){
				success=false;
			}
		}
		if(string(argv[i])=="-di"&& i<argc-1){
			try{
				Concord_Dist_Idx=atoi(argv[i+1]);
			}
			catch(exception& e){
				success=false;
			}
		}
		if(string(argv[i])=="-G" && i<argc-1){
			if(string(argv[i+1])=="0")
				Print_Graph=0;
			else if(string(argv[i+1])=="1")
				Print_Graph=1;
			else
				success=false;
		}
		if(string(argv[i])=="-CO" && i<argc-1){
			if(string(argv[i+1])=="0")
				Print_Components_Ordering=0;
			else if(string(argv[i+1])=="1")
				Print_Components_Ordering=1;
			else
				success=false;
		}
		if(string(argv[i])=="-TO" && i<argc-1){
			if(string(argv[i+1])=="0")
				Print_Total_Ordering=0;
			else if(string(argv[i+1])=="1")
				Print_Total_Ordering=1;
			else
				success=false;
		}
		if(string(argv[i])=="-RG" && i<argc-1){
			if(string(argv[i+1])=="0")
				Print_Rearranged_Genome=0;
			else if(string(argv[i+1])=="1")
				Print_Rearranged_Genome=1;
			else
				success=false;
		}
	}
	if(Input_BAM=="" || Output_Prefix==""){
		print_help();
		success=false;
	}
	if(Input_FASTA=="" && Print_Rearranged_Genome){
		printf("reference FASTA needed to output rearranged genome sequence.\n");
		success=false;
	}
	if(!specify_mq && Input_Chim_BAM!="")
		Min_MapQual=255;
	if(!success)
		printf("Check your argument.\n");
	return success;
};