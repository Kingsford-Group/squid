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
	printf("Usage: GenerateNewGenome [--indirect] -b <Input_BAM> -f <Input_FASTA> -g <Input_GRAPH> -c <Input_Components> -o <Output_Prefix>\n");
};

int main(int argc, char* argv[]){
	if(argc!=11)
		print_help_generatenewgenome();
	else{
		map<string, int> RefTable;
		vector<string> RefName;
		vector<int> RefLength;
		vector<string> RefSequence;
		string Input_GRAPH;
		string Input_COMP;
		bool DirectMode=true;

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
			if(string(argv[i])=="--indirect")
				DirectMode=false;
		}

		BuildRefName(Input_BAM, RefName, RefTable, RefLength);
		bool canbuild=BuildRefSeq(Input_FASTA, RefTable, RefLength, RefSequence);
		if(canbuild){
			SegmentGraph_t SegmentGraph(Input_GRAPH);
			vector< vector<int> > Components=ReadComponents(Input_COMP);
			if(!DirectMode){
				vector< pair<int, int> > Node_NewChr; Node_NewChr.resize(SegmentGraph.vNodes.size());
				for(int i=0; i<Components.size(); i++)
					for(int j=0; j<Components[i].size(); j++)
						Node_NewChr[abs(Components[i][j])-1]=make_pair(i, j);

				vector<Edge_t> UnSatisfied; UnSatisfied.reserve((int)SegmentGraph.vEdges.size()/2);
				int concordthresh=50000;
				for(int i=0; i<SegmentGraph.vEdges.size(); i++){
					int ind1=SegmentGraph.vEdges[i].Ind1, ind2=SegmentGraph.vEdges[i].Ind2;
					bool issatisfied=false;
					pair<int,int> pos1=Node_NewChr[SegmentGraph.vEdges[i].Ind1], pos2=Node_NewChr[SegmentGraph.vEdges[i].Ind2];
					if(pos1.first==pos2.first && pos1.second<pos2.second && SegmentGraph.vEdges[i].Head1==(Components[pos1.first][pos1.second]<0) && SegmentGraph.vEdges[i].Head2==(Components[pos2.first][pos2.second]>0)){
						issatisfied=true;
					}
					else if(pos1.first==pos2.first && pos1.second>pos2.second && SegmentGraph.vEdges[i].Head2==(Components[pos2.first][pos2.second]<0) && SegmentGraph.vEdges[i].Head1==(Components[pos1.first][pos1.second]>0)){
						issatisfied=true;
					}
					if(!issatisfied)
						UnSatisfied.push_back(SegmentGraph.vEdges[i]);
				}
				vector<Edge_t> Diff; Diff.resize(SegmentGraph.vEdges.size());
				vector<Edge_t>::iterator it=set_difference(SegmentGraph.vEdges.begin(), SegmentGraph.vEdges.end(), UnSatisfied.begin(), UnSatisfied.end(), Diff.begin());
				Diff.resize(distance(Diff.begin(), it));
				SegmentGraph.vEdges=Diff;
				SegmentGraph.UpdateNodeLink();
				Diff.clear();

				vector<int> MergeNode(SegmentGraph.vNodes.size(), -1);
				int curNode=0;
				int rightmost=0;
				for(int i=0; i<MergeNode.size(); i++){
					int minDisInd2;
					vector<Edge_t> thisDiscordantEdges, tmpDiscordantEdges;
					if(i!=0 && SegmentGraph.vNodes[i].Chr!=SegmentGraph.vNodes[i-1].Chr && curNode==MergeNode[i-1])
						curNode++;
					for(vector<Edge_t*>::iterator itedge=SegmentGraph.vNodes[i].HeadEdges.begin(); itedge!=SegmentGraph.vNodes[i].HeadEdges.end(); itedge++){
						if(SegmentGraph.IsDiscordant(*itedge))
							thisDiscordantEdges.push_back(**itedge);
						else
							rightmost=(rightmost>max((*itedge)->Ind1, (*itedge)->Ind2))?rightmost:max((*itedge)->Ind1, (*itedge)->Ind2);
					}
					for(vector<Edge_t*>::iterator itedge=SegmentGraph.vNodes[i].TailEdges.begin(); itedge!=SegmentGraph.vNodes[i].TailEdges.end(); itedge++){
						if(SegmentGraph.IsDiscordant(*itedge))
							thisDiscordantEdges.push_back(**itedge);
						else
							rightmost=(rightmost>max((*itedge)->Ind1, (*itedge)->Ind2))?rightmost:max((*itedge)->Ind1, (*itedge)->Ind2);
					}
					if(thisDiscordantEdges.size()!=0){ // remove discordant edges that are in the same group
						minDisInd2 = (thisDiscordantEdges[0].Ind1==i)?thisDiscordantEdges[0].Ind2:(i+20);
						tmpDiscordantEdges.push_back(thisDiscordantEdges[0]);
						for(int k=0; k+1<thisDiscordantEdges.size(); k++){
							const Edge_t& edge1= thisDiscordantEdges[k];
							const Edge_t& edge2= thisDiscordantEdges[k+1];
							bool samegroup = ((edge1.Ind1==i && edge2.Ind1==i) || (edge1.Ind2==i && edge2.Ind2==i));
							if(!(abs(edge1.Ind1-edge2.Ind1)<=Concord_Dist_Idx && abs(edge1.Ind2-edge2.Ind2)<=Concord_Dist_Idx && edge1.Head1==edge2.Head1 && edge1.Head2==edge2.Head2))
								samegroup=false;
							if(!samegroup)
								tmpDiscordantEdges.push_back(thisDiscordantEdges[k+1]);
							int tmpmin = (edge2.Ind1==i)?edge2.Ind2:(i+20);
							if(tmpmin<minDisInd2)
								minDisInd2=tmpmin;
						}
						thisDiscordantEdges=tmpDiscordantEdges;
					}

					if(MergeNode[i]==-1){
						if(thisDiscordantEdges.size()==0 && i<rightmost) // node with only concordant edges, and further connect right nodes
							MergeNode[i]=curNode;
						else if(thisDiscordantEdges.size()==0 && i==rightmost){ // node with only concordant edges, and doesn't connect right nodes
							MergeNode[i]=curNode;
							curNode++;
							rightmost++;
						}
						else{ // first node with discordant in its equivalent group
							if(i!=0 && curNode==MergeNode[i-1])
								curNode++;
							vector<Edge_t> nextDiscordantEdges;
							int j=i+1;
							for(; j<MergeNode.size() && j<i+20 && j<minDisInd2; j++){ // find the first right nodes with discordant edges, possibly in the equivalent group
								for(vector<Edge_t*>::iterator itedge=SegmentGraph.vNodes[j].HeadEdges.begin(); itedge!=SegmentGraph.vNodes[j].HeadEdges.end(); itedge++)
									if(SegmentGraph.IsDiscordant(*itedge))
										nextDiscordantEdges.push_back(**itedge);
								for(vector<Edge_t*>::iterator itedge=SegmentGraph.vNodes[j].TailEdges.begin(); itedge!=SegmentGraph.vNodes[j].TailEdges.end(); itedge++)
									if(SegmentGraph.IsDiscordant(*itedge))
										nextDiscordantEdges.push_back(**itedge);
								if(nextDiscordantEdges.size()!=0)
									break;
							}
							bool equivalent = (nextDiscordantEdges.size()!=0); // check whether indeed in the same equivalent group
							if(nextDiscordantEdges.size()!=0){
								tmpDiscordantEdges.clear();
								tmpDiscordantEdges.push_back(nextDiscordantEdges[0]);
								for(int k=0; k+1<nextDiscordantEdges.size(); k++){
									const Edge_t& edge1=nextDiscordantEdges[k];
									const Edge_t& edge2=nextDiscordantEdges[k+1];
									bool samegroup = ((edge1.Ind1==j && edge2.Ind1==j) || (edge1.Ind2==j && edge2.Ind2==j));
									if(!(abs(edge1.Ind1-edge2.Ind1)<=Concord_Dist_Idx && abs(edge1.Ind2-edge2.Ind2)<=Concord_Dist_Idx && edge1.Head1==edge2.Head1 && edge1.Head2==edge2.Head2))
										samegroup=false;
									if(!samegroup)
										tmpDiscordantEdges.push_back(nextDiscordantEdges[k+1]);
								}
								nextDiscordantEdges=tmpDiscordantEdges;
								tmpDiscordantEdges.clear();
								vector<bool> thisEQ(thisDiscordantEdges.size(), false);
								vector<bool> nextEQ(nextDiscordantEdges.size(), false);
								for(int k=0; k<thisDiscordantEdges.size(); k++)
									for(int l=0; l<nextDiscordantEdges.size(); l++){
										const Edge_t& edge1=thisDiscordantEdges[k];
										const Edge_t& edge2=nextDiscordantEdges[l];
										if(edge1.Ind2>edge2.Ind1 && edge2.Ind2>edge1.Ind1 && abs(edge1.Ind1-edge2.Ind1)<=Concord_Dist_Idx && abs(edge1.Ind2-edge2.Ind2)<=Concord_Dist_Idx && edge1.Head1==edge2.Head1 && edge1.Head2==edge2.Head2){
											thisEQ[k]=true;
											nextEQ[l]=true;
										}
									}
								for(int k=0; k<thisEQ.size(); k++)
									if(!thisEQ[k])
										equivalent=false;
								for(int l=0; l<nextEQ.size(); l++)
									if(!nextEQ[l])
										equivalent=false;
							}
							if(!equivalent){ // equivalend group has only the current item
								MergeNode[i]=curNode;
								curNode++;
							}
							else{ // equivalent group has other items than current item
								for(int k=i; k<=j; k++)
									MergeNode[k]=curNode;
							}
							rightmost=i+1;
						}
					}
					else if(thisDiscordantEdges.size()!=0){ // later nodes with discordant edges in its equivalent group
						vector<Edge_t> nextDiscordantEdges;
						int j=i+1;
						for(; j<MergeNode.size() && j<i+20 && j<minDisInd2; j++){
							for(vector<Edge_t*>::iterator itedge=SegmentGraph.vNodes[j].HeadEdges.begin(); itedge!=SegmentGraph.vNodes[j].HeadEdges.end(); itedge++)
								if(SegmentGraph.IsDiscordant(*itedge))
									nextDiscordantEdges.push_back(**itedge);
							for(vector<Edge_t*>::iterator itedge=SegmentGraph.vNodes[j].TailEdges.begin(); itedge!=SegmentGraph.vNodes[j].TailEdges.end(); itedge++)
								if(SegmentGraph.IsDiscordant(*itedge))
									nextDiscordantEdges.push_back(**itedge);
							if(nextDiscordantEdges.size()!=0)
								break;
						}
						bool equivalent = (nextDiscordantEdges.size()!=0);
						if(nextDiscordantEdges.size()!=0){
							tmpDiscordantEdges.clear();
							tmpDiscordantEdges.push_back(thisDiscordantEdges[0]);
							for(int k=0; k+1<thisDiscordantEdges.size(); k++){
								const Edge_t& edge1= thisDiscordantEdges[k];
								const Edge_t& edge2= thisDiscordantEdges[k+1];
								if(!(abs(edge1.Ind1-edge2.Ind1)<=Concord_Dist_Idx && abs(edge1.Ind2-edge2.Ind2)<=Concord_Dist_Idx && edge1.Head1==edge2.Head1 && edge1.Head2==edge2.Head2))
									tmpDiscordantEdges.push_back(thisDiscordantEdges[k+1]);
							}
							thisDiscordantEdges=tmpDiscordantEdges;
							tmpDiscordantEdges.clear();
							tmpDiscordantEdges.push_back(nextDiscordantEdges[0]);
							for(int k=0; k+1<nextDiscordantEdges.size(); k++){
								const Edge_t& edge1=nextDiscordantEdges[k];
								const Edge_t& edge2=nextDiscordantEdges[k+1];
								if(!(abs(edge1.Ind1-edge2.Ind1)<=Concord_Dist_Idx && abs(edge1.Ind2-edge2.Ind2)<=Concord_Dist_Idx && edge1.Head1==edge2.Head1 && edge1.Head2==edge2.Head2))
									tmpDiscordantEdges.push_back(nextDiscordantEdges[k+1]);
							}
							nextDiscordantEdges=tmpDiscordantEdges;
							tmpDiscordantEdges.clear();
							vector<bool> thisEQ(thisDiscordantEdges.size(), false);
							vector<bool> nextEQ(nextDiscordantEdges.size(), false);
							for(int k=0; k<thisDiscordantEdges.size(); k++)
								for(int l=0; l<nextDiscordantEdges.size(); l++){
									const Edge_t& edge1=thisDiscordantEdges[k];
									const Edge_t& edge2=nextDiscordantEdges[l];
									if(edge1.Ind2>edge2.Ind1 && edge2.Ind2>edge1.Ind1 && abs(edge1.Ind1-edge2.Ind1)<=Concord_Dist_Idx && abs(edge1.Ind2-edge2.Ind2)<=Concord_Dist_Idx && edge1.Head1==edge2.Head1 && edge1.Head2==edge2.Head2){
										thisEQ[k]=true;
										nextEQ[l]=true;
									}
								}
							for(int k=0; k<thisEQ.size(); k++)
								if(!thisEQ[k])
									equivalent=false;
							for(int l=0; l<nextEQ.size(); l++)
								if(!nextEQ[l])
									equivalent=false;
						}
						if(!equivalent) // not able to find the next one with discordant edges in the same equivalent group
							curNode++;
						else{
							for(int k=i; k<=j; k++)
								MergeNode[k]=curNode;
						}
						rightmost=i+1;
					}
				}

				vector<Node_t> newvNodes; newvNodes.reserve(curNode);
				vector<Edge_t> newvEdges; newvEdges.reserve(SegmentGraph.vEdges.size());
				int ind=0;
				while(ind<MergeNode.size()){
					int j=ind;
					for(; j<MergeNode.size() && MergeNode[j]==MergeNode[ind]; j++){}
					Node_t tmpnode(SegmentGraph.vNodes[ind].Chr, SegmentGraph.vNodes[ind].Position, SegmentGraph.vNodes[j-1].Position+SegmentGraph.vNodes[j-1].Length-SegmentGraph.vNodes[ind].Position);
					newvNodes.push_back(tmpnode);
					ind=j;
				}
				for(int i=0; i<SegmentGraph.vEdges.size(); i++){
					const Edge_t& edge=SegmentGraph.vEdges[i];
					if(MergeNode[edge.Ind1]!=MergeNode[edge.Ind2]){
						Edge_t tmpedge(MergeNode[edge.Ind1], edge.Head1, MergeNode[edge.Ind2], edge.Head2, edge.Weight);
						newvEdges.push_back(tmpedge);
					}
				}

				SegmentGraph.vNodes=newvNodes;
				SegmentGraph.vEdges.clear();
				sort(newvEdges.begin(), newvEdges.end());
				for(int i=0; i<newvEdges.size(); i++){
					if(SegmentGraph.vEdges.size()==0 || !(newvEdges[i]==SegmentGraph.vEdges.back()))
						SegmentGraph.vEdges.push_back(newvEdges[i]);
					else
						SegmentGraph.vEdges.back().Weight+=newvEdges[i].Weight;
				}
				SegmentGraph.UpdateNodeLink();
				SegmentGraph.ConnectedComponent();

				Components=SegmentGraph.Ordering();
				Components=SegmentGraph.SortComponents(Components);
				Components=SegmentGraph.MergeSingleton(Components, RefLength);
				Components=SegmentGraph.SortComponents(Components);
				Components=SegmentGraph.MergeComponents(Components);
			}

			Components=SegmentGraph.SortComponents(Components);
			OutputNewGenome(SegmentGraph, Components, RefSequence, RefName, Output_Prefix+"_genome.fa");
		}
	}
}