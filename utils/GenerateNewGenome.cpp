#include "../src/SingleBamRec.h"
#include "../src/ReadRec.h"
#include "../src/BPNode.h"
#include "../src/BPEdge.h"
#include "../src/SegmentGraph.h"
#include "../src/WriteIO.h"
#include "../src/Config.h"

using namespace std;

void print_help_generatenewgenome(){
	printf("GenerateNewGenome takes in squid output and original genome FASTA file, and generate a new FASTA file corresponding to rearrangement.\n");
	printf("Usage: GenerateNewGenome -b <Input_BAM> -f <Input_FASTA> -g <Input_GRAPH> -c <Input_Components> -o <Output_Prefix>\n");
};

int main(int argc, char* argv[]){
	if(argc!=6)
		print_help_generatenewgenome();
	else{
		map<string, int> RefTable;
		vector<string> RefName;
		vector<int> RefLength;
		vector<string> RefSequence;
		string Input_GRAPH;
		string Input_COMP;

		for(int i=1; i<argc; i++){
			if(string(argv[i])=="-b" && i<argc-1)
				Input_BAM=string(argv[i+1]);
			if(string(argv[i])=="-f" && i<argc-1)
				Input_FASTA=string(argv[i+1]);
			if(string(argv[i])=="-g" && i<argc-1)
				Input_GRAPH=string(argv[i+1]);
			if(string(argv[i])=="-c" && i<argc-1)
				Input_COMP=string(argv[i+1]);
			if(string(argv[i])=="-o" && i<argc-1)
				Output_Prefix=string(argv[i+1]);
		}

		BuildRefName(Input_BAM, RefName, RefTable, RefLength);
		bool canbuild=BuildRefSeq(Input_FASTA, RefTable, RefLength, RefSequence);
		if(canbuild){
			SegmentGraph_t SegmentGraph(Input_GRAPH);
			vector< vector<int> > Components=ReadComponents(Input_COMP);
			OutputNewGenome(SegmentGraph, Components, RefSequence, RefName, Output_Prefix+"_genome.fa");
		}
	}
}