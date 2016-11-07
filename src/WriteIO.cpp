#include "WriteIO.h"

using namespace std;

void WriteComponents(string outputfile, vector< vector<int> > Components){
	ofstream output(outputfile, ios::out);
	output<<"# component_id\tnodes\n";
	for(int i=0; i<Components.size(); i++){
		output<<i<<'\t';
		for(int j=0; j<Components[i].size()-1; j++)
			output<<Components[i][j]<<",";
		output<<Components[i][Components[i].size()-1]<<endl;
	}
	output.close();
};

void WriteBEDPE(string outputfile, SegmentGraph_t& SegmentGraph, vector< vector<int> >& Components, vector< pair<int, int> >& Node_NewChr, vector<string>& RefName){
	sort(SegmentGraph.vEdges.begin(), SegmentGraph.vEdges.end(),  [](Edge_t a, Edge_t b){return a.Weight>b.Weight;});
	ofstream output(outputfile, ios::out);
	output<<"# chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tstrand1\tstrand2\tweight\n";
	for(int i=0; i<SegmentGraph.vEdges.size(); i++){
		int ind1=SegmentGraph.vEdges[i].Ind1, ind2=SegmentGraph.vEdges[i].Ind2;
		bool flag_chr=(SegmentGraph.vNodes[ind1].Chr==SegmentGraph.vNodes[ind2].Chr);
		bool flag_ori=(SegmentGraph.vEdges[i].Head1==false && SegmentGraph.vEdges[i].Head2==true);
		bool flag_dist=(SegmentGraph.vNodes[ind2].Position-SegmentGraph.vNodes[ind1].Position-SegmentGraph.vNodes[ind1].Length<=concorddis || ind2-ind1<=concordidx);
		if(!flag_chr || !flag_ori || !flag_dist){
			pair<int,int> pos1=Node_NewChr[SegmentGraph.vEdges[i].Ind1];
			pair<int,int> pos2=Node_NewChr[SegmentGraph.vEdges[i].Ind2];
			if(pos1.first==pos2.first && pos1.second<pos2.second && SegmentGraph.vEdges[i].Head1==(Components[pos1.first][pos1.second]<0) && SegmentGraph.vEdges[i].Head2==(Components[pos2.first][pos2.second]>0)){
				output<<RefName[SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Chr]<<'\t';
				output<<SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Position<<'\t';
				output<<(SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Position+SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Length)<<'\t';
				output<<RefName[SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Chr]<<'\t';
				output<<SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Position<<'\t';
				output<<(SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Position+SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Length)<<'\t';
				output<<".\t";
				output<<(SegmentGraph.vEdges[i].Head1?"-\t":"+\t");
				output<<(SegmentGraph.vEdges[i].Head2?"+\t":"-\t");
				output<<SegmentGraph.vEdges[i].Weight<<endl;
			}
			else if(pos1.first==pos2.first && pos1.second>pos2.second && SegmentGraph.vEdges[i].Head2==(Components[pos2.first][pos2.second]<0) && SegmentGraph.vEdges[i].Head1==(Components[pos1.first][pos1.second]>0)){
				output<<RefName[SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Chr]<<'\t';
				output<<SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Position<<'\t';
				output<<(SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Position+SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Length)<<'\t';
				output<<RefName[SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Chr]<<'\t';
				output<<SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Position<<'\t';
				output<<(SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Position+SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Length)<<'\t';
				output<<".\t";
				output<<(SegmentGraph.vEdges[i].Head1?"-\t":"+\t");
				output<<(SegmentGraph.vEdges[i].Head2?"+\t":"-\t");
				output<<SegmentGraph.vEdges[i].Weight<<endl;
			}
		}
	}
	output.close();
};