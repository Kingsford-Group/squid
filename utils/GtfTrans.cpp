#include "GtfTrans.h"

using namespace std;

string Transcript_t::GetFeature(string line, string FeatureName){
	size_t start=line.find(FeatureName);
	size_t stop=line.find(";", start+FeatureName.size()+2);
	return line.substr(start+FeatureName.size()+2, stop-start-FeatureName.size()-3);
};

void ReadGTF(string filename, vector<Transcript_t>& vTrans){
	vTrans.clear();
	vector<Interval_t> AllExons;
	vector<Interval_t> AllCDSs;
	ifstream input(filename);
	string line;
	while(getline(input, line)){
		if(line[0]=='#')
			continue;
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t"));
		if(strs[2]=="transcript"){
			string ID=Transcript_t::GetFeature(line, "transcript_id");
			string Name=Transcript_t::GetFeature(line, "gene_name");
			Transcript_t tmp(ID, Name, strs[0], strs[6][0], stoi(strs[3])-1, stoi(strs[4]));
			vTrans.push_back(tmp);
		}
		else if(strs[2]=="exon"){
			string TID=Transcript_t::GetFeature(line, "transcript_id");
			Interval_t tmp(strs[0], stoi(strs[3])-1, stoi(strs[4]), strs[6][0], TID);
			AllExons.push_back(tmp);
		}
		else if(strs[2]=="CDS"){
			string TID=Transcript_t::GetFeature(line, "transcript_id");
			Interval_t tmp(strs[0], stoi(strs[3])-1, stoi(strs[4]), strs[6][0], TID);
			AllCDSs.push_back(tmp);
		}
	}
	input.close();

	sort(AllExons.begin(), AllExons.end(), Interval_t::TransIDComp);
	sort(AllCDSs.begin(), AllCDSs.end(), Interval_t::TransIDComp);
	sort(vTrans.begin(), vTrans.end(), Transcript_t::TransIDComp);

	vector<Interval_t>::iterator itexon=AllExons.begin();
	vector<Interval_t>::iterator itcds=AllCDSs.begin();
	for(vector<Transcript_t>::iterator it=vTrans.begin(); it!=vTrans.end(); it++){
		while(itexon!=AllExons.end() && itexon->TransID==it->TransID){
			it->vExon.push_back(*itexon);
			itexon++;
		}
		while(itcds!=AllCDSs.end() && itcds->TransID==it->TransID){
			it->vCDS.push_back(*itcds);
			itcds++;
		}
	}
	sort(vTrans.begin(), vTrans.end());
	for(vector<Transcript_t>::iterator it=vTrans.begin(); it!=vTrans.end(); it++){
		sort(it->vExon.begin(), it->vExon.end());
		sort(it->vCDS.begin(), it->vCDS.end());
	}
};