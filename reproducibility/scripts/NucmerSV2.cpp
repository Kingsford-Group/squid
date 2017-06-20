#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <limits>
#include <ctime>
#include <cmath>
#include <map>
#include <boost/algorithm/string.hpp>

using namespace std;

class Alignment_t{
public:
	string RefChr;
	int ReadStart;
	int ReadEnd;
	int RefStart;
	int RefEnd;
public:
	Alignment_t(){};
	Alignment_t(string RefChr, int ReadStart, int ReadEnd, int RefStart, int RefEnd): RefChr(RefChr), ReadStart(ReadStart), ReadEnd(ReadEnd), RefStart(RefStart), RefEnd(RefEnd){};
	void Inverse(){
		int tmp=ReadStart;
		ReadStart=ReadEnd; ReadEnd=tmp;
		tmp=RefStart;
		RefStart=RefEnd; RefEnd=tmp;
	};
	bool IsReadInverse(){
		return (ReadStart>ReadEnd);
	};
	bool operator < (const Alignment_t& rhs) const{
		return min(ReadStart, ReadEnd)<min(rhs.ReadStart, rhs.ReadEnd);
	};
	static bool IsConcordant(Alignment_t& a, Alignment_t& b){
		if(a.RefChr!=b.RefChr || a.IsReadInverse()!=b.IsReadInverse())
			return false;
		else if(!a.IsReadInverse()){
			if(min(a.ReadStart, a.ReadEnd)<=min(b.ReadStart, b.ReadEnd) && a.RefStart>b.RefStart)
				return false;
			else if(min(a.ReadStart, a.ReadEnd)>=min(b.ReadStart, b.ReadEnd) && a.RefStart<b.RefStart)
				return false;
		}
		else{
			if(min(a.ReadStart, a.ReadEnd)<=min(b.ReadStart, b.ReadEnd) && a.RefStart<b.RefStart)
				return false;
			else if(min(a.ReadStart, a.ReadEnd)>=min(b.ReadStart, b.ReadEnd) && a.RefStart>b.RefStart)
				return false;
		}
		return true;
	};
};

struct BP_t{
	string RefChr;
	int StartPos, EndPos;
	bool IsLeft;

	BP_t(){};
	BP_t(string RefChr, int StartPos, int EndPos, bool IsLeft): RefChr(RefChr), StartPos(StartPos), EndPos(EndPos), IsLeft(IsLeft){};

	bool operator < (const BP_t& rhs) const{
		if(RefChr!=rhs.RefChr)
			return RefChr<rhs.RefChr;
		else if(StartPos!=rhs.StartPos)
			return StartPos<rhs.StartPos;
		else if(EndPos!=rhs.EndPos)
			return EndPos<rhs.EndPos;
		else
			return IsLeft<rhs.IsLeft;
	};
	bool operator == (const BP_t& rhs) const{
		return (RefChr==rhs.RefChr)&&(StartPos==rhs.StartPos)&&(EndPos==rhs.EndPos)&&(IsLeft==rhs.IsLeft);
	};
};

class SV_t{
public:
	BP_t BP1, BP2;
public:
	SV_t(){};
	SV_t(BP_t BP1, BP_t BP2): BP1(BP1), BP2(BP2){};

	bool operator < (const SV_t& rhs) const{
		if(!(BP1==rhs.BP1))
			return BP1<rhs.BP1;
		else
			return BP2<rhs.BP2;
	};
	bool operator == (const SV_t& rhs) const{
		return (BP1==rhs.BP1)&&(BP2==rhs.BP2);
	};
};

bool Overlap(const Alignment_t& a,  const vector<Alignment_t>& v){
	int length=0;
	for(int i=0; i<v.size(); i++){
		int vmin=min(v[i].ReadStart, v[i].ReadEnd), vmax=max(v[i].ReadStart, v[i].ReadEnd);
		int amin=min(a.ReadStart, a.ReadEnd), amax=max(a.ReadStart, a.ReadEnd);
		if(amin<vmin && amax>=vmin)
			length+=min(vmax, amax)-vmin;
		else if(amin>=vmin && amin<=vmax)
			length+=min(vmax, amax)-amin;
	}
	if(length<0.25*abs(a.ReadStart-a.ReadEnd))
		return false;
	else
		return true;
};

bool Overlap(const Alignment_t& a, const Alignment_t& b){
	int amin=min(a.ReadStart, a.ReadEnd), amax=max(a.ReadStart, a.ReadEnd);
	int bmin=min(b.ReadStart, b.ReadEnd), bmax=max(b.ReadStart, b.ReadEnd);
	int length=0;
	if(amin<=bmin && amax>=bmin)
		length+=min(amax, bmax)-bmin;
	else if(amin>=bmin && amin<=bmax)
		length+=min(amax, bmax)-amin;
	if(length<0.25*abs(amax-amin) && length<0.25*abs(bmax-bmin))
		return false;
	else return true;
};

vector<SV_t> FindDiscordant(const vector<Alignment_t>& Alignments){
	vector<Alignment_t> copy=Alignments;
	sort(copy.begin(), copy.end());
	vector<SV_t> SVs;
	int overstart=-1;
	for(int i=0; i+1<copy.size(); i++){
		bool isoverlap=Overlap(copy[i], copy[i+1]);
		if(isoverlap && overstart==-1)
			overstart=i;
		else if(!isoverlap && overstart==-1){
			if(!Alignment_t::IsConcordant(copy[i], copy[i+1])){
				BP_t BP1(copy[i].RefChr, copy[i].RefStart, copy[i].RefEnd, copy[i].IsReadInverse());
				BP_t BP2(copy[i+1].RefChr, copy[i+1].RefStart, copy[i+1].RefEnd, !copy[i+1].IsReadInverse());
				SV_t tmp(BP1, BP2);
				SVs.push_back(tmp);
			}
		}
		else if(!isoverlap){
			bool concordexist=false;
			for(int j=overstart; j<i+1; j++)
				if(Alignment_t::IsConcordant(copy[j], copy[i+1]))
					concordexist=true;
			if(!concordexist){
				for(int j=overstart; j<i+1; j++){
					BP_t BP1(copy[j].RefChr, copy[j].RefStart, copy[j].RefEnd, copy[j].IsReadInverse());
					BP_t BP2(copy[i+1].RefChr, copy[i+1].RefStart, copy[i+1].RefEnd, !copy[i+1].IsReadInverse());
					SV_t tmp(BP1, BP2);
					SVs.push_back(tmp);
				}
			}
			overstart=-1;
		}
	}
	return SVs;
};

vector<SV_t> SVbyNucmer(string filename){
	vector<SV_t> AllSV;

	ifstream input(filename);
	string line;
	bool startreadind=false;
	string curTransName="";
	vector<Alignment_t> Alignments;
	while(getline(input, line)){
		if(line.substr(0,5)=="=====")
			startreadind=true;
		else if(startreadind){
			vector<string> strs;
			boost::trim_if(line, boost::is_any_of("\t "));
			boost::split(strs, line, boost::is_any_of(" \t"), boost::token_compress_on);
			if(curTransName==""){
				Alignments.clear();
				curTransName=strs[12];
			}
			else if(strs[12]!=curTransName){
				vector<SV_t> SVs=FindDiscordant(Alignments);
				if(SVs.size()!=0)
					AllSV.insert(AllSV.end(), SVs.begin(), SVs.end());
				Alignments.clear();
				curTransName=strs[12];
			}
			Alignment_t tmp(strs[11], stoi(strs[3]), stoi(strs[4]), stoi(strs[0]), stoi(strs[1]));
			Alignments.push_back(tmp);
		}
	}
	vector<SV_t> SVs=FindDiscordant(Alignments);
	if(SVs.size()!=0)
		AllSV.insert(AllSV.end(), SVs.begin(), SVs.end());
	input.close();
	return AllSV;
};

vector<SV_t> FilterRepeat(vector<SV_t>& AllSV){
	int thresh=20;
	sort(AllSV.begin(), AllSV.end());
	vector<SV_t> newSV; newSV.reserve(AllSV.size());
	int count=1;
	for(int i=0; i<AllSV.size(); i++){
		if(newSV.size()==0)
			newSV.push_back(AllSV[i]);
		else{
			bool isnear=true;
			if(newSV.back().BP1.RefChr!=AllSV[i].BP1.RefChr || newSV.back().BP2.RefChr!=AllSV[i].BP2.RefChr)
				isnear=false;
			else if(newSV.back().BP1.IsLeft!=AllSV[i].BP1.IsLeft || newSV.back().BP2.IsLeft!=AllSV[i].BP2.IsLeft)
				isnear=false;
			else{
				int insertedbp1=(newSV.back().BP1.IsLeft)?newSV.back().BP1.StartPos:newSV.back().BP1.EndPos;
				int insertedbp2=(newSV.back().BP2.IsLeft)?newSV.back().BP2.StartPos:newSV.back().BP2.EndPos;
				int todobp1=(AllSV[i].BP1.IsLeft)?AllSV[i].BP1.StartPos:AllSV[i].BP1.EndPos;
				int todobp2=(AllSV[i].BP2.IsLeft)?AllSV[i].BP2.StartPos:AllSV[i].BP2.EndPos;
				if(abs(insertedbp1-todobp1)>thresh || abs(insertedbp2-todobp2)>thresh)
					isnear=false;
			}
			if(!isnear){
				newSV.push_back(AllSV[i]);
				count=1;
			}
			else{
				newSV.back().BP1.StartPos=round(newSV.back().BP1.StartPos*count+AllSV[i].BP1.StartPos)/(count+1);
				newSV.back().BP1.EndPos=round(newSV.back().BP1.EndPos*count+AllSV[i].BP1.EndPos)/(count+1);
				newSV.back().BP2.StartPos=round(newSV.back().BP2.StartPos*count+AllSV[i].BP2.StartPos)/(count+1);
				newSV.back().BP2.EndPos=round(newSV.back().BP2.EndPos*count+AllSV[i].BP2.EndPos)/(count+1);
				count++;
			}
		}
	}
	newSV.reserve(newSV.size());
	return newSV;
};

void WriteBedpe(string outputfile, vector<SV_t>& AllSV){
	ofstream output(outputfile, ios::out);
	for(vector<SV_t>::iterator it=AllSV.begin(); it!=AllSV.end(); it++)
		output<<it->BP1.RefChr<<"\t"<<it->BP1.StartPos<<"\t"<<it->BP1.EndPos<<"\t"<<it->BP2.RefChr<<"\t"<<it->BP2.StartPos<<"\t"<<it->BP2.EndPos<<"\t.\t.\t"<<(it->BP1.IsLeft?"-\t":"+\t")<<(it->BP2.IsLeft?"-\n":"+\n");
	output.close();
};

int main(int argc, char* argv[]){
	string InputNucmer(argv[1]);
	string OutputFile(argv[2]);

	vector<SV_t> AllSV=SVbyNucmer(InputNucmer);
	AllSV=FilterRepeat(AllSV);
	WriteBedpe(OutputFile, AllSV);
}