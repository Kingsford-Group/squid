#include "WriteIO.h"

using namespace std;

vector< vector<int> > ReadComponents(string file){
	ifstream input(file);
	string line;
	vector< vector<int> > Components;
	while(getline(input,line)){
		if(line[0]=='#')
			continue;
		else{
			size_t start=line.find_first_of('\t');
			line=line.substr(start+1);
			vector<string> strs;
			boost::split(strs, line, boost::is_any_of(","));
			vector<int> tmp;
			for(int i=0; i<strs.size(); i++)
				tmp.push_back(stoi(strs[i]));
			Components.push_back(tmp);
		}
	}
	input.close();
	return Components;
};

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
		bool flag_dist=(SegmentGraph.vNodes[ind2].Position-SegmentGraph.vNodes[ind1].Position-SegmentGraph.vNodes[ind1].Length<=Concord_Dist_Pos || ind2-ind1<=Concord_Dist_Idx);
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
				output<<(SegmentGraph.vEdges[i].Head2?"-\t":"+\t");
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
				output<<(SegmentGraph.vEdges[i].Head2?"-\t":"+\t");
				output<<SegmentGraph.vEdges[i].Weight<<endl;
			}
		}
	}
	output.close();
};

void OutputNewGenome(SegmentGraph_t& SegmentGraph, vector< vector<int> >& Components, const vector<string>& RefSequence, const vector<string>& RefName, string outputfile){
	ofstream output(outputfile, ios::out);
	for(int i=0; i<Components.size(); i++){
		string info="PA:", seq, tmpseq;
		for(int j=0; j<Components[i].size(); j++){
			int k;
			for(k=j+1; k<Components[i].size() && Components[i][k]-Components[i][k-1]==1 && SegmentGraph.vNodes[abs(Components[i][j])-1].Chr==SegmentGraph.vNodes[abs(Components[i][k])-1].Chr; k++){}
			if(Components[i][j]>0){
				int curChr=SegmentGraph.vNodes[abs(Components[i][j])-1].Chr;
				int curStart=SegmentGraph.vNodes[abs(Components[i][j])-1].Position;
				int curLen=SegmentGraph.vNodes[abs(Components[i][k-1])-1].Position+SegmentGraph.vNodes[abs(Components[i][k-1])-1].Length-SegmentGraph.vNodes[abs(Components[i][j])-1].Position;
				tmpseq=RefSequence[curChr].substr(curStart, curLen);
				info+="{"+RefName[curChr]+","+to_string(curStart)+","+to_string(curLen)+"}";
			}
			else{
				int curChr=SegmentGraph.vNodes[abs(Components[i][k-1])-1].Chr;
				int curStart=SegmentGraph.vNodes[abs(Components[i][k-1])-1].Position;
				int curLen=SegmentGraph.vNodes[abs(Components[i][j])-1].Position+SegmentGraph.vNodes[abs(Components[i][j])-1].Length-SegmentGraph.vNodes[abs(Components[i][k-1])-1].Position;
				tmpseq=RefSequence[curChr].substr(curStart, curLen);
				info+="{"+RefName[curChr]+","+to_string(curStart)+","+to_string(curLen)+"}";
			}
			if(Components[i][j]<0)
				ReverseComplement(tmpseq.begin(), tmpseq.end());
			seq+=tmpseq;
			info+=((Components[i][j]<0)?"R-":"F-");
			j=k-1;
		}
		info=info.substr(0, info.size()-1);
		output<<">chr"<<(i+1)<<'\t'<<"LN:"<<seq.size()<<'\t'<<info<<endl;
		output<<seq<<endl;
	}
	output.close();
};