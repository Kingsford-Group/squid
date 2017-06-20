#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include "SV.h"
#include "GtfTrans.h"

using namespace std;

map<string, int> RefTable;
map<int, int> RefLength;
vector<string> RefName;

void BuildRef(string fafile){
	ifstream input(fafile);
	string line, name;
	int count=0, len=0;
	while(getline(input, line)){
		if(line[0]=='>'){
			if(len!=0){
				RefTable[name]=count;
				RefLength[count]=len;
				count++;
			}
			len=0;
			size_t stop=line.find_first_of(' ');
			name=line.substr(1, stop-1);
			RefName.push_back(name);
		}
		else{
			len+=line.size();
		}
	}
	if(len!=0){
		RefTable[name]=count;
		RefLength[count]=len;
	}
};
string PrintRefInfo(){
	string info;
	int i=0;
	for(map<string,int>::iterator it=RefTable.begin(); it!=RefTable.end(); it++, i++){
		info+=it->first+"\tID="+to_string(it->second)+"\tLN:"+to_string(RefLength[it->second])+"\tStored Name:"+RefName[i]+'\n';
	}
	return info;
};

int main(int argc, char* argv[]){
	if(argc>1){
		string fafile(argv[1]);
		BuildRef(fafile);
		cout<<PrintRefInfo();
		string svfilestr[]={argv[2],argv[2],argv[2],argv[2],argv[2]};
		svfilestr[0]+=(svfilestr[0].back()=='/')?"deletions.csv":"/deletions.csv";
		svfilestr[1]+=(svfilestr[1].back()=='/')?"insertions.csv":"/insertions.csv";
		svfilestr[2]+=(svfilestr[2].back()=='/')?"inversions.csv":"/inversions.csv";
		svfilestr[3]+=(svfilestr[3].back()=='/')?"tandemDuplications.csv":"/tandemDuplications.csv";
		svfilestr[4]+=(svfilestr[4].back()=='/')?"translocations.csv":"/translocations.csv";
		SV_t sv(RefTable, 5, svfilestr);
		if(argc>4){
			string str(argv[5]);
			vector<Transcript_t> vTrans;
			vector< vector< pair<AffectExonRelativePos, AffectExonRelativePos> > > AffectedPos;
			vector<string> TransName;
			sv.IsCoveredbyReads(argv[4], RefTable, RefLength, TransName, (str=="old"));
			ReadGTF(argv[6], vTrans, RefTable);
			FilterGTF(TransName, vTrans);
			FindAffected(vTrans, sv, AffectedPos);
			//sv.WriteFilterednewSVPos(RefName, RefLength, argv[3]);
			//sv.WriteFilteredoldSVPos(RefName, RefLength, argv[3]);
			if(!sv.updated)
				sv.Update2NewPos(RefLength);
			WriteBreakpoint(argv[3], sv, RefName, AffectedPos);
		}
		else
			sv.WritenewSVPos(RefName, RefLength, argv[3]);

		/*string sv2filestr[]={argv[3],argv[3],argv[3],argv[3],argv[3]};
		sv2filestr[0]+=(sv2filestr[0].back()=='/')?"deletions.csv":"/deletions.csv";
		sv2filestr[1]+=(sv2filestr[1].back()=='/')?"insertions.csv":"/insertions.csv";
		sv2filestr[2]+=(sv2filestr[2].back()=='/')?"inversions.csv":"/inversions.csv";
		sv2filestr[3]+=(sv2filestr[3].back()=='/')?"tandemDuplications.csv":"/tandemDuplications.csv";
		sv2filestr[4]+=(sv2filestr[4].back()=='/')?"translocations.csv":"/translocations.csv";
		SV_t sv2(RefTable, 5, sv2filestr);
		bool Intersect=sv.IsIntersect(sv2);
		cout<<(Intersect?"intersect":"no intersection")<<endl;
		string outputfile(argv[2]);
		outputfile+=(outputfile.back()=='/')?"simulated_newpos_r1.csv":"/simulated_newpos_r1.csv";
		cout<<outputfile<<endl;
		sv.WritenewSVPos(RefName, RefLength, outputfile);
		vector< tuple<int,int,int,int>> BreakPoint=sv.UpdateNextRound2new(sv2, RefLength);
		sv.WriteNextRoundBPPos(RefName, BreakPoint, argv[4]);*/
	}
	else{
		cout<<"Input: <genome.fasta> <simulated_sv_r1 folder> <simulated_sv_r2 folder> <output_breakpoint.dat>\n";
	}
}
