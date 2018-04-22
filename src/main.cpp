/*
Part of SQUID transcriptomic structural variation detector
(c) 2017 by  Cong Ma, Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "SingleBamRec.h"
#include "ReadRec.h"
#include "BPNode.h"
#include "BPEdge.h"
#include "SegmentGraph.h"
#include "WriteIO.h"
#include "Config.h"

using namespace std;

int main(int argc, char* argv[]){
	time_t CurrentTime;
	string CurrentTimeStr;
	
	bool success=parse_arguments(argc, argv);

	if(success){
		map<string, int> RefTable;
		vector<string> RefName;
		vector<int> RefLength;

		BuildRefName(Input_BAM, RefName, RefTable, RefLength);
		for(map<string,int>::iterator it=RefTable.begin(); it!=RefTable.end(); it++)
			cout<<"Reference name "<<it->first<<"\t-->\t"<<it->second<<endl;

		SBamrecord_t Chimrecord;
		if(Input_Chim_BAM.size()!=0){
			BuildChimericSBamRecord(Chimrecord, RefName, Input_Chim_BAM);
		}

		SegmentGraph_t SegmentGraph(RefLength, Chimrecord, Input_BAM);

		if(Print_Graph)
			SegmentGraph.OutputGraph(Output_Prefix+"_graph.txt");
		vector< vector<int> > Components=SegmentGraph.Ordering();
		if(Print_Components_Ordering)
			WriteComponents(Output_Prefix+"_component_pri.txt", Components);

		Components=SegmentGraph.SortComponents(Components);
		Components=SegmentGraph.MergeSingleton(Components, RefLength);
		Components=SegmentGraph.SortComponents(Components);
		Components=SegmentGraph.MergeComponents(Components);

		vector< pair<int, int> > Node_NewChr; Node_NewChr.resize(SegmentGraph.vNodes.size());
		for(unsigned int i=0; i<Components.size(); i++)
			for(unsigned int j=0; j<Components[i].size(); j++)
				Node_NewChr[abs(Components[i][j])-1]=make_pair(i, j);

		if(Print_Total_Ordering)
			WriteComponents(Output_Prefix+"_component.txt", Components);
		map<Edge_t, vector< pair<int,int> > > ExactBP;
		SegmentGraph.ExactBreakpoint(Chimrecord, ExactBP);
		map<Edge_t, vector< pair<int,int> > > ExactBP_concord_support;
		SegmentGraph.ExactBPConcordantSupport(Input_BAM, Chimrecord, ExactBP, ExactBP_concord_support);
		SegmentGraph.DeMultiplyDisEdges();
		WriteBEDPE(Output_Prefix+"_sv.txt", SegmentGraph, Components, Node_NewChr, RefName, ExactBP, ExactBP_concord_support);

		if(Print_Rearranged_Genome){
			vector<string> RefSequence;
			bool canbuild=BuildRefSeq(Input_FASTA, RefTable, RefLength, RefSequence);
			if(canbuild)
				OutputNewGenome(SegmentGraph, Components, RefSequence, RefName, Output_Prefix+"_genome.fa");
		}

		time(&CurrentTime);
		CurrentTimeStr=ctime(&CurrentTime);
		cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] Done."<<endl;

	}
}
