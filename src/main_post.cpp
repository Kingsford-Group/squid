#include "SingleBamRec.h"
#include "ReadRec.h"
#include "BPNode.h"
#include "BPEdge.h"
#include "SegmentGraph2.h"
#include "WriteIO.h"
#include "Config.h"

using namespace std;

int main(int argc, char* argv[]){
	bool success=parse_arguments(argc, argv);

	if(success){
		map<string, int> RefTable;
		vector<string> RefName;
		vector<int> RefLength;

		BuildRefName(Input_BAM, RefName, RefTable, RefLength);
		for(map<string,int>::iterator it=RefTable.begin(); it!=RefTable.end(); it++)
			cout<<"Reference name "<<it->first<<"\t-->\t"<<it->second<<endl;

		SegmentGraph_t SegmentGraph(Input_FASTA);

		vector< vector<int> > Components=SegmentGraph.Ordering();

		Components=SegmentGraph.SortComponents(Components);
		Components=SegmentGraph.MergeSingleton(Components, RefLength);
		Components=SegmentGraph.SortComponents(Components);
		Components=SegmentGraph.MergeComponents(Components);

		vector< pair<int, int> > Node_NewChr; Node_NewChr.resize(SegmentGraph.vNodes.size());
		for(int i=0; i<Components.size(); i++)
			for(int j=0; j<Components[i].size(); j++)
				Node_NewChr[abs(Components[i][j])-1]=make_pair(i, j);

		if(Print_Total_Ordering)
			WriteComponents(Output_Prefix+"_component.txt", Components);

		if(Print_Rearranged_Genome){
			vector<string> RefSequence;
			bool canbuild=BuildRefSeq(Input_FASTA, RefTable, RefLength, RefSequence);
			if(canbuild)
				OutputNewGenome(SegmentGraph, Components, RefSequence, RefName, Output_Prefix+"_genome.fa");
		}

		int concordthresh=50000;
		sort(SegmentGraph.vEdges.begin(), SegmentGraph.vEdges.end(),  [](Edge_t a, Edge_t b){return a.Weight>b.Weight;});
		ofstream output3(Output_Prefix+"_discordantedges.txt", ios::out);
		output3<<"# ";
		for(int i=0; i<argc; i++)
			output3<<argv[i]<<" ";
		output3<<endl;
		output3<<"# Ind1\tNode1\tHead1\tInd2\tNode2\tHead2\tWeight\n";
		for(int i=0; i<SegmentGraph.vEdges.size(); i++){
			int ind1=SegmentGraph.vEdges[i].Ind1, ind2=SegmentGraph.vEdges[i].Ind2;
			if(SegmentGraph.vNodes[ind1].Chr!=SegmentGraph.vNodes[ind2].Chr || (SegmentGraph.vNodes[ind2].Position-SegmentGraph.vNodes[ind1].Position-SegmentGraph.vNodes[ind1].Length>Concord_Dist_Pos && ind2-ind1>Concord_Dist_Idx) || SegmentGraph.vEdges[i].Head1!=false || SegmentGraph.vEdges[i].Head2!=true){
				pair<int,int> pos1=Node_NewChr[SegmentGraph.vEdges[i].Ind1], pos2=Node_NewChr[SegmentGraph.vEdges[i].Ind2];
				if(pos1.first==pos2.first && pos1.second<pos2.second && SegmentGraph.vEdges[i].Head1==(Components[pos1.first][pos1.second]<0) && SegmentGraph.vEdges[i].Head2==(Components[pos2.first][pos2.second]>0)){
					output3<<SegmentGraph.vEdges[i].Ind1<<'\t'<<SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Chr<<','<<SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Position<<','<<SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Length<<'\t'<<(SegmentGraph.vEdges[i].Head1?"H\t":"T\t")<<SegmentGraph.vEdges[i].Ind2<<'\t';
					output3<<SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Chr<<','<<SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Position<<','<<SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Length<<'\t'<<(SegmentGraph.vEdges[i].Head2?"H\t":"T\t")<<SegmentGraph.vEdges[i].Weight<<endl;
				}
				else if(pos1.first==pos2.first && pos1.second>pos2.second && SegmentGraph.vEdges[i].Head2==(Components[pos2.first][pos2.second]<0) && SegmentGraph.vEdges[i].Head1==(Components[pos1.first][pos1.second]>0)){
					output3<<SegmentGraph.vEdges[i].Ind1<<'\t'<<SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Chr<<','<<SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Position<<','<<SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Length<<'\t'<<(SegmentGraph.vEdges[i].Head1?"H\t":"T\t")<<SegmentGraph.vEdges[i].Ind2<<'\t';
					output3<<SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Chr<<','<<SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Position<<','<<SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Length<<'\t'<<(SegmentGraph.vEdges[i].Head2?"H\t":"T\t")<<SegmentGraph.vEdges[i].Weight<<endl;
				}
			}
		}
		output3.close();
	}
}
