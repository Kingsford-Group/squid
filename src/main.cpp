#include "SingleBamRec.h"
#include "ReadRec.h"
#include "BPNode.h"
#include "BPEdge.h"
#include "SegmentGraph.h"

using namespace std;

int concorddis=50000, concorddisl=50000;
int concordidx=20, concordidxl=20;
int weightcutoff=5;
string inputbam="";
string inputfasta="";
string inputchimbam="";
string outputpref="";
int LowPhredLenThresh=10;
bool PhredType=0;

int main(int argc, char* argv[]){
    int opt;
    while((opt=getopt(argc,argv,"b:f:c:o:d:i:w:"))!=EOF){
        switch(opt){
            case 'b': inputbam=(string)optarg; break;
            case 'f': inputfasta=(string)optarg; break;
            case 'c': inputchimbam=(string)optarg; break;
            case 'o': outputpref=(string)optarg; break;
            case 'd': concorddis=atoi(optarg); concorddisl=concorddis; break;
            case 'i': concordidx=atoi(optarg); concordidxl=concordidx; break;
            case 'w': weightcutoff=atoi(optarg); break;
        }
    }

    map<string, int> RefTable;
    vector<string> RefName;
    vector<int> RefLength;
    vector<string> RefSequence;

    BuildRefName(inputbam, RefName, RefTable);
    BuildReference(inputfasta, RefTable, RefLength, RefSequence);
    for(map<string,int>::iterator it=RefTable.begin(); it!=RefTable.end(); it++)
        cout<<"Reference name "<<it->first<<"\t-->\t"<<it->second<<endl;

    SBamrecord_t SBamrecord;
    if(inputchimbam.size()!=0){
        SBamrecord=BuildMainSBamRecord(RefTable, inputbam);
        SBamrecord=BuildChimericSBamRecord(SBamrecord, RefTable, inputchimbam);
    }
    else
        SBamrecord=BuildBWASBamRecord(RefTable, inputbam);

    vector< vector<int> > Read_Node;
    SegmentGraph_t SegmentGraph(RefLength, SBamrecord, Read_Node, weightcutoff);
    cout<<"Read_Node size correct? "<<(Read_Node.size()==SBamrecord.size())<<endl;
    for(int i=0; i<SBamrecord.size(); i++){
        if((int)SBamrecord[i].FirstRead.size()+(int)SBamrecord[i].SecondMate.size()!=(int)Read_Node[i].size())
            cout<<"wrong size\n";
        for(int j=0; j<Read_Node[i].size(); j++)
            if(Read_Node[i][j]!=-1 && Read_Node[i][j]>=SegmentGraph.vNodes.size())
                cout<<"out of range\n";
    }

    SegmentGraph.OutputConnectedComponent(outputpref+"/graph.txt");
    vector< vector<int> > Components=SegmentGraph.Ordering();

    ofstream output0(outputpref+"/component_pri.txt", ios::out);
    output0<<"# component_id\tnodes\n";
    for(int i=0; i<Components.size(); i++){
        output0<<i<<'\t';
        for(int j=0; j<Components[i].size()-1; j++)
            output0<<Components[i][j]<<",";
        output0<<Components[i][Components[i].size()-1]<<endl;
    }
    output0.close();

    map<int,int> NewIndex;
    vector<Node_t> NewNodeChr;
    vector< vector<int> > LowSupportNode;
    vector<int> ReferenceNode;
    vector<bool> RelativePosition;
    SegmentGraph.SimplifyComponents(Components, NewIndex, NewNodeChr, LowSupportNode, ReferenceNode, RelativePosition);
    Components=SegmentGraph.SortComponents(Components);
    Components=SegmentGraph.MergeSingleton(Components, RefLength, NewNodeChr);
    SegmentGraph.DesimplifyComponents(Components, NewIndex, LowSupportNode, ReferenceNode, RelativePosition);
    Components=SegmentGraph.SortComponents(Components);
    Components=SegmentGraph.MergeComponents(Components);

    vector< pair<int, int> > Node_NewChr; Node_NewChr.resize(SegmentGraph.vNodes.size());
    for(int i=0; i<Components.size(); i++)
        for(int j=0; j<Components[i].size(); j++)
            Node_NewChr[abs(Components[i][j])-1]=make_pair(i, j);

    ofstream output(outputpref+"/component.txt", ios::out);
    output<<"# component_id\tnodes\n";
    for(int i=0; i<Components.size(); i++){
        output<<i<<'\t';
        for(int j=0; j<Components[i].size()-1; j++)
            output<<Components[i][j]<<",";
        output<<Components[i][Components[i].size()-1]<<endl;
    }
    output.close();

    int concordthresh=50000;
    sort(SegmentGraph.vEdges.begin(), SegmentGraph.vEdges.end(),  [](Edge_t a, Edge_t b){return a.Weight>b.Weight;});
    ofstream output3(outputpref+"/discordantedges.txt", ios::out);
    output3<<"# ";
    for(int i=0; i<argc; i++)
        output3<<argv[i]<<" ";
    output3<<endl;
    output3<<"# Ind1\tNode1\tHead1\tInd2\tNode2\tHead2\tWeight\n";
    for(int i=0; i<SegmentGraph.vEdges.size(); i++){
        int ind1=SegmentGraph.vEdges[i].Ind1, ind2=SegmentGraph.vEdges[i].Ind2;
        if(SegmentGraph.vNodes[ind1].Chr!=SegmentGraph.vNodes[ind2].Chr || (SegmentGraph.vNodes[ind2].Position-SegmentGraph.vNodes[ind1].Position-SegmentGraph.vNodes[ind1].Length>concordthresh && ind2-ind1>20) || SegmentGraph.vEdges[i].Head1!=false || SegmentGraph.vEdges[i].Head2!=true){
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
