#include "GtfTrans.h"

using namespace std;

string Transcript_t::GetFeature(string line, string FeatureName){
	size_t start=line.find(FeatureName);
	size_t stop=line.find(";", start+FeatureName.size()+2);
	return line.substr(start+FeatureName.size()+2, stop-start-FeatureName.size()-3);
};

bool Transcript_t::IsAffectedbySV(SimpleSV_t& sv, vector< pair<AffectExonRelativePos, AffectExonRelativePos> >& aepos){
	bool flag=false;
	if(sv.Type!=SimpleSVType::INS){
		if(Chr==sv.RefID && TxEnd>=sv.StartPos && TxStart<=sv.EndPos){
			for(int i=0; i<vExon.size(); i++){
				if(vExon[i].StartPos<sv.StartPos && vExon[i].EndPos>sv.StartPos){
					AffectExonRelativePos tmp1(true, true, 0, sv.ID);
					AffectExonRelativePos tmp2(true, false, 0, sv.ID);
					aepos.push_back(make_pair(tmp1, tmp2));
					flag=true;
				}
				if(vExon[i].StartPos<sv.EndPos && vExon[i].EndPos>sv.EndPos){
					AffectExonRelativePos tmp1(false, true, 0, sv.ID);
					AffectExonRelativePos tmp2(false, false, 0, sv.ID);
					aepos.push_back(make_pair(tmp1, tmp2));
					flag=true;
				}
				if(i<vExon.size()-1){
					if(vExon[i].EndPos<=sv.StartPos && vExon[i+1].StartPos>=sv.StartPos && vExon[i+1].EndPos<sv.EndPos){
						AffectExonRelativePos tmp1(true, true, vExon[i].EndPos-sv.StartPos, sv.ID);
						AffectExonRelativePos tmp2(true, false, vExon[i+1].StartPos-sv.StartPos, sv.ID);
						aepos.push_back(make_pair(tmp1, tmp2));
						flag=true;
					}
					if(vExon[i].EndPos<=sv.EndPos && vExon[i+1].StartPos>=sv.EndPos && vExon[i].EndPos>sv.StartPos){
						AffectExonRelativePos tmp1(false, true, vExon[i].EndPos-sv.EndPos, sv.ID);
						AffectExonRelativePos tmp2(false, false, vExon[i+1].StartPos-sv.EndPos, sv.ID);
						aepos.push_back(make_pair(tmp1, tmp2));
						flag=true;
					}
				}
			}
		}
	}
	else{
		if(Chr==sv.RefID && TxEnd>=sv.StartPos && TxStart<=sv.StartPos){
			for(int i=0; i<vExon.size(); i++){
				if(vExon[i].StartPos<sv.StartPos && vExon[i].EndPos>sv.StartPos){
					AffectExonRelativePos tmp1(true, true, 0, sv.ID);
					AffectExonRelativePos tmp2(true, false, 0, sv.ID);
					aepos.push_back(make_pair(tmp1, tmp2));
					flag=true;
				}
				if(i<vExon.size()-1){
					if(vExon[i].EndPos<=sv.StartPos && vExon[i+1].StartPos>sv.StartPos){
						AffectExonRelativePos tmp1(true, true, vExon[i].EndPos-sv.StartPos, sv.ID);
						AffectExonRelativePos tmp2(true, false, vExon[i+1].StartPos-sv.StartPos, sv.ID);
						aepos.push_back(make_pair(tmp1, tmp2));
						flag=true;
					}
				}
			}
		}
	}
	return flag;
};

bool Transcript_t::IsAffectedbySV(TRA_t& sv, vector< pair<AffectExonRelativePos, AffectExonRelativePos> >& aepos){
	bool flag=false;
	if(Chr==sv.Ref1 && TxStart<=sv.Pos1 && TxEnd>=sv.Pos1){
		for(int i=0; i<vExon.size(); i++){
			if(vExon[i].StartPos<sv.Pos1 && vExon[i].EndPos>sv.Pos1){
				AffectExonRelativePos tmp1(true, true, 0, sv.ID);
				AffectExonRelativePos tmp2(true, false, 0, sv.ID);
				aepos.push_back(make_pair(tmp1, tmp2));
				flag=true;
			}
			if(i<vExon.size()-1 && vExon[i].EndPos<=sv.Pos1 && vExon[i+1].StartPos>=sv.Pos1){
				AffectExonRelativePos tmp1(true, true, vExon[i].EndPos-sv.Pos1, sv.ID);
				AffectExonRelativePos tmp2(true, false, vExon[i+1].StartPos-sv.Pos1, sv.ID);
				aepos.push_back(make_pair(tmp1, tmp2));
				flag=true;
			}
		}
	}
	if(Chr==sv.Ref2 && TxStart<=sv.Pos2 && TxEnd>=sv.Pos2){
		for(int i=0; i<vExon.size(); i++){
			if(vExon[i].StartPos<sv.Pos2 && vExon[i].EndPos>sv.Pos2){
				AffectExonRelativePos tmp1(false, true, 0, sv.ID);
				AffectExonRelativePos tmp2(false, false, 0, sv.ID);
				aepos.push_back(make_pair(tmp1, tmp2));
				flag=true;
			}
			if(i<vExon.size()-1 && vExon[i].EndPos<=sv.Pos2 && vExon[i+1].StartPos>=sv.Pos2){
				AffectExonRelativePos tmp1(false, true, vExon[i].EndPos-sv.Pos2, sv.ID);
				AffectExonRelativePos tmp2(false, false, vExon[i+1].StartPos-sv.Pos2, sv.ID);
				aepos.push_back(make_pair(tmp1, tmp2));
				flag=true;
			}
		}
	}
	return flag;
};

void ReadGTF(string filename, vector<Transcript_t>& vTrans, map<string,int>& RefTable){
	vTrans.clear();
	vector<Interval_t> AllExons;
	ifstream input(filename);
	string line;
	while(getline(input, line)){
		if(line[0]=='#')
			continue;
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t"));
		if(RefTable.find(strs[0])==RefTable.end())
			continue;
		if(strs[2]=="exon"){
			string TID=Transcript_t::GetFeature(line, "transcript_id");
			// string Number=Transcript_t::GetFeature(line, "exon_number");
			// Interval_t tmp(RefTable.at(strs[0]), stoi(strs[3])-1, stoi(strs[4]), strs[6][0], TID, stoi(Number));
			Interval_t tmp(RefTable.at(strs[0]), stoi(strs[3])-1, stoi(strs[4]), strs[6][0], TID);
			AllExons.push_back(tmp);
		}
	}
	input.close();

	sort(AllExons.begin(), AllExons.end(), Interval_t::ExonNumberComp);

	Transcript_t tmpTrans;
	for(vector<Interval_t>::iterator it=AllExons.begin(); it!=AllExons.end(); it++){
		if(tmpTrans.vExon.size()==0 || it->TransID==tmpTrans.vExon.back().TransID)
			tmpTrans.vExon.push_back(*it);
		else{
			sort(tmpTrans.vExon.begin(), tmpTrans.vExon.end());
			tmpTrans.TransID=tmpTrans.vExon.front().TransID;
			tmpTrans.Chr=tmpTrans.vExon.front().Chr;
			tmpTrans.TxStart=tmpTrans.vExon.front().StartPos;
			tmpTrans.TxEnd=tmpTrans.vExon.back().EndPos;
			vTrans.push_back(tmpTrans);
			tmpTrans.vExon.clear();
		}
	}
	sort(vTrans.begin(), vTrans.end());
};

void FilterGTF(string fluxpro, vector<Transcript_t>& vTrans){
	ifstream input(fluxpro);
	string line;
	vector<string> Expressed;
	while(getline(input, line)){
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t"));
		if(strs[5]!="0")
			Expressed.push_back(strs[1]);
	}
	input.close();
	sort(Expressed.begin(), Expressed.end());
	vector<Transcript_t> newvTrans;
	newvTrans.reserve(vTrans.size());
	for(vector<Transcript_t>::iterator it=vTrans.begin(); it!=vTrans.end(); it++)
		if(binary_search(Expressed.begin(), Expressed.end(), it->TransID))
			newvTrans.push_back((*it));
	vTrans=newvTrans;
};

void FilterGTF(vector<string>& TransName, vector<Transcript_t>& vTrans){
	vector<Transcript_t> newvTrans;
	newvTrans.reserve(vTrans.size());
	for(vector<Transcript_t>::iterator it=vTrans.begin(); it!=vTrans.end(); it++)
		if(binary_search(TransName.begin(), TransName.end(), it->TransID))
			newvTrans.push_back((*it));
	vTrans=newvTrans;
};

void FindAffected(vector<Transcript_t>& vTrans, SV_t& SVs, vector< vector< pair<AffectExonRelativePos, AffectExonRelativePos> > >& AffectedPos){
	AffectedPos.clear();
	for(int i=0; i<SVs.vSimpleSV.size()+SVs.vTRA.size(); i++){
		vector< pair<AffectExonRelativePos, AffectExonRelativePos> > tmp;
		AffectedPos.push_back(tmp);
	}
	for(vector<Transcript_t>::iterator it=vTrans.begin(); it!=vTrans.end(); it++){
		for(int i=0; i<SVs.vSimpleSV.size(); i++){
			vector< pair<AffectExonRelativePos, AffectExonRelativePos> > aepos;
			bool flag=it->IsAffectedbySV(SVs.vSimpleSV[i], aepos);
			if(flag)
				AffectedPos[i].insert(AffectedPos[i].end(), aepos.begin(), aepos.end());
		}
		for(int i=0; i<SVs.vTRA.size(); i++){
			vector< pair<AffectExonRelativePos, AffectExonRelativePos> > aepos;
			bool flag=it->IsAffectedbySV(SVs.vTRA[i], aepos);
			int k=i+(int)SVs.vSimpleSV.size();
			if(flag)
				AffectedPos[k].insert(AffectedPos[k].end(), aepos.begin(), aepos.end());
		}
	}
	for(int i=0; i<AffectedPos.size(); i++){
		if(AffectedPos[i].size()==0)
			continue;
		sort(AffectedPos[i].begin(), AffectedPos[i].end(), [](pair<AffectExonRelativePos,AffectExonRelativePos> a, pair<AffectExonRelativePos,AffectExonRelativePos> b){if(!(a.first==b.first)) return a.first<b.first; else return a.second<b.second;});
		vector< pair<AffectExonRelativePos,AffectExonRelativePos> >::iterator endit;
		endit=unique(AffectedPos[i].begin(), AffectedPos[i].end(), [](pair<AffectExonRelativePos,AffectExonRelativePos> a, pair<AffectExonRelativePos,AffectExonRelativePos> b){return a.first==b.first && a.second==b.second;});
		AffectedPos[i].resize(distance(AffectedPos[i].begin(), endit));
	}
};

void WriteBreakpoint(string filename, SV_t& SVs, vector<string>& RefName, vector< vector< pair<AffectExonRelativePos, AffectExonRelativePos> > >& AffectedPos){
	ofstream output(filename, ios::out);
	output<<"# chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsvID\tscore\tstrand1\tstrand2\n";
	for(int i=0; i<SVs.vSimpleSV.size(); i++){
		if(AffectedPos[i].size()==0 || !SVs.vSimpleReads[i])
			continue;
		bool flag=true;
		// filter out single duplications and deletions, but keep insertions (transposons / transpositions)
		if(SVs.vSimpleSV[i].Type==SimpleSVType::INS && ((i==0)?true:(SVs.vSimpleSV[i-1].ID!=SVs.vSimpleSV[i].ID)) && ((i==SVs.vSimpleSV.size()-1)?true:(SVs.vSimpleSV[i+1].ID!=SVs.vSimpleSV[i].ID)))
			flag=false;
		else if(SVs.vSimpleSV[i].Type==SimpleSVType::DEL && ((i==0)?true:(SVs.vSimpleSV[i-1].ID!=SVs.vSimpleSV[i].ID)) && ((i==SVs.vSimpleSV.size()-1)?true:(SVs.vSimpleSV[i+1].ID!=SVs.vSimpleSV[i].ID)))
			flag=false;
		if(!flag)
			continue;
		if(SVs.vSimpleSV[i].Type==SimpleSVType::INS){
			const SimpleSV_t& sv=SVs.vSimpleSV[i];
			SimpleSV_t svpair;
			if(i!=0 && SVs.vSimpleSV[i-1].ID==sv.ID)
				svpair=SVs.vSimpleSV[i-1];
			else if(i!=SVs.vSimpleSV.size()-1 && SVs.vSimpleSV[i+1].ID==sv.ID)
				svpair=SVs.vSimpleSV[i+1];
			for(int j=0; j<AffectedPos[i].size(); j++){
				const pair<AffectExonRelativePos, AffectExonRelativePos>& tmpBPs=AffectedPos[i][j];
				assert(tmpBPs.first.svID==sv.ID && tmpBPs.second.svID==sv.ID);
				if(tmpBPs.first.IsStartBP){
					string Chr1=RefName[sv.RefID];
					string Chr2=RefName[svpair.RefID];
					int bp1=sv.StartPos+tmpBPs.first.RelativePos;
					int bp2=svpair.StartPos+tmpBPs.second.RelativePos;
					output<<Chr1<<"\t"<<(bp1-1)<<"\t"<<bp1<<"\t"<<Chr2<<"\t"<<bp2<<"\t"<<(bp2+1)<<"\t"<<sv.ID<<"\t.\t+\t-"<<endl;
				}
				else{
					string Chr1=RefName[svpair.RefID];
					string Chr2=RefName[sv.RefID];
					int bp1=svpair.EndPos+tmpBPs.first.RelativePos;
					int bp2=sv.StartPos+tmpBPs.second.RelativePos;
					output<<Chr1<<"\t"<<(bp1-1)<<"\t"<<bp1<<"\t"<<Chr2<<"\t"<<bp2<<"\t"<<(bp2+1)<<"\t"<<sv.ID<<"\t.\t+\t-"<<endl;
				}
			}
		}
		// else if(SVs.vSimpleSV[i].Type==SimpleSVType::DEL){
		// 	const SimpleSV_t& sv=SVs.vSimpleSV[i];
		// 	SimpleSV_t svpair;
		// 	if(i!=0 && SVs.vSimpleSV[i-1].ID==sv.ID)
		// 		svpair=SVs.vSimpleSV[i-1];
		// 	else if(i!=SVs.vSimpleSV.size()-1 && SVs.vSimpleSV[i+1].ID==sv.ID)
		// 		svpair=SVs.vSimpleSV[i+1];
		// 	for(int j=0; j<AffectedPos[i].size(); j++){
		// 		const pair<AffectExonRelativePos, AffectExonRelativePos>& tmpBPs=AffectedPos[i][j];
		// 		assert(tmpBPs.first.svID==sv.ID && tmpBPs.second.svID==sv.ID);
		// 		if(tmpBPs.first.IsStartBP){
		// 			string Chr=RefName[sv.RefID];
		// 			int bp1=sv.StartPos+tmpBPs.first.RelativePos;
		// 			int bp2=sv.EndPos+tmpBPs.second.RelativePos;
		// 			output<<Chr<<"\t"<<(bp1-1)<<"\t"<<bp1<<"\t"<<Chr<<"\t"<<bp2<<"\t"<<(bp2+1)<<"\t"<<"\t"<<sv.ID<<"\t.\t+\t-"<<endl;
		// 		}
		// 	}
		// }
		else if(SVs.vSimpleSV[i].Type==SimpleSVType::INV){
			const SimpleSV_t& sv=SVs.vSimpleSV[i];
			for(int j=0; j<AffectedPos[i].size(); j++){
				const pair<AffectExonRelativePos, AffectExonRelativePos>& tmpBPs=AffectedPos[i][j];
				assert(tmpBPs.first.svID==sv.ID && tmpBPs.second.svID==sv.ID);
				if(tmpBPs.first.IsStartBP){
					string Chr=RefName[sv.RefID];
					int bp1=sv.StartPos+tmpBPs.first.RelativePos;
					int bp2=sv.EndPos-tmpBPs.second.RelativePos;
					output<<Chr<<"\t"<<(bp1-1)<<"\t"<<bp1<<"\t"<<Chr<<"\t"<<(bp2-1)<<"\t"<<bp2<<"\t"<<sv.ID<<"\t.\t+\t+"<<endl;
				}
				else{
					string Chr=RefName[sv.RefID];
					int bp1=sv.StartPos-tmpBPs.first.RelativePos;
					int bp2=sv.EndPos+tmpBPs.second.RelativePos;
					output<<Chr<<"\t"<<(bp1-1)<<"\t"<<bp1<<"\t"<<Chr<<"\t"<<(bp2-1)<<"\t"<<bp2<<"\t"<<sv.ID<<"\t.\t-\t-"<<endl;
				}
			}
		}
	}
	for(int i=0; i<SVs.vTRA.size(); i++){
		int k=i+(int)SVs.vSimpleSV.size();
		if(AffectedPos[k].size()==0 || !SVs.vTRAReads[i])
			continue;
		const TRA_t& sv=SVs.vTRA[i];
		for(int j=0; j<AffectedPos[k].size(); j++){
			const pair<AffectExonRelativePos, AffectExonRelativePos>& tmpBPs=AffectedPos[k][j];
			assert(tmpBPs.first.svID==sv.ID && tmpBPs.second.svID==sv.ID);
			if(tmpBPs.first.IsStartBP){
				string Chr1, Chr2;
				int bp1, bp2;
				bool IsForward1, IsForward2;
				if(sv.DT1==DirType::right){
					Chr1=RefName[sv.Ref2];
					Chr2=RefName[sv.Ref1];
					bp2=sv.Pos1+tmpBPs.second.RelativePos;
					IsForward2=false;
					if(sv.DT1==sv.DT2){
						bp1=sv.Pos2+tmpBPs.first.RelativePos;
						IsForward1=true;
					}
					else{
						bp1=sv.Pos2-tmpBPs.first.RelativePos;
						IsForward1=false;
					}
				}
				else{
					Chr1=RefName[sv.Ref1];
					Chr2=RefName[sv.Ref2];
					bp1=sv.Pos1+tmpBPs.first.RelativePos;
					IsForward1=true;
					if(sv.DT1==sv.DT2){
						bp2=sv.Pos2+tmpBPs.second.RelativePos;
						IsForward2=false;
					}
					else{
						bp2=sv.Pos2-tmpBPs.second.RelativePos;
						IsForward2=true;
					}
				}
				output<<Chr1<<"\t"<<(IsForward1?(bp1-1):bp1)<<"\t"<<(IsForward1?bp1:(bp1+1))<<"\t"<<Chr2<<"\t"<<(IsForward2?(bp2-1):bp2)<<"\t"<<(IsForward2?bp2:(bp2+1))<<"\t"<<sv.ID<<"\t.\t"<<(IsForward1?"+":"-")<<"\t"<<(IsForward2?"+":"-")<<endl;
			}
			else{
				string Chr1, Chr2;
				int bp1, bp2;
				bool IsForward1, IsForward2;
				if(sv.DT2==DirType::right){
					Chr1=RefName[sv.Ref1];
					Chr2=RefName[sv.Ref2];
					bp2=sv.Pos2+tmpBPs.second.RelativePos;
					IsForward2=false;
					if(sv.DT1==sv.DT2){
						bp1=sv.Pos1+tmpBPs.first.RelativePos;
						IsForward1=true;
					}
					else{
						bp1=sv.Pos2-tmpBPs.first.RelativePos;
						IsForward1=false;
					}
				}
				else{
					Chr1=RefName[sv.Ref2];
					Chr2=RefName[sv.Ref1];
					bp1=sv.Pos2+tmpBPs.first.RelativePos;
					IsForward1=true;
					if(sv.DT1==sv.DT2){
						bp2=sv.Pos1+tmpBPs.second.RelativePos;
						IsForward2=false;
					}
					else{
						bp2=sv.Pos1-tmpBPs.second.RelativePos;
						IsForward2=true;
					}
				}
				output<<Chr1<<"\t"<<(IsForward1?(bp1-1):bp1)<<"\t"<<(IsForward1?bp1:(bp1+1))<<"\t"<<Chr2<<"\t"<<(IsForward2?(bp2-1):bp2)<<"\t"<<(IsForward2?bp2:(bp2+1))<<"\t"<<sv.ID<<"\t.\t"<<(IsForward1?"+":"-")<<"\t"<<(IsForward2?"+":"-")<<endl;
			}
		}
	}
	output.close(); 
};