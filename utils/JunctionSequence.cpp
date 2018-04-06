#include "SingleBamRec.h"
#include "ReadRec.h"
#include <limits>

using namespace std;

static std::map<char,char> Nucleotide={{'A','T'},{'C','G'},{'G','C'},{'T','A'},{'R','Y'},{'Y','R'},{'S','W'},{'W','S'},{'K','M'},{'M','K'},{'B','V'},{'V','B'},{'D','H'},{'H','D'}, {'N','N'}, {'.','.'},{'-','-'}};

void ReverseComplement(string::iterator itbegin, string::iterator itend){
	for(string::iterator it=itbegin; it!=itend; it++)
		*it=Nucleotide[toupper(*it)];
	std::reverse(itbegin, itend);
};

class Breakpoint_t{
public:
	int Chr;
	int StartPos;
	int EndPos;
	bool IsLeft;
public:
	Breakpoint_t(){};
	Breakpoint_t(int Chr, int StartPos, int EndPos, bool IsLeft): Chr(Chr), StartPos(StartPos), EndPos(EndPos), IsLeft(IsLeft) {};

	bool operator < (const Breakpoint_t& rhs) const {
		if(Chr!=rhs.Chr)
			return Chr<rhs.Chr;
		else if(StartPos!=rhs.StartPos)
			return StartPos<rhs.StartPos;
		else if(EndPos!=rhs.EndPos)
			return EndPos<rhs.EndPos;
		else
			return IsLeft<rhs.IsLeft;
	};

	bool operator > (const Breakpoint_t& rhs) const {
		if(Chr!=rhs.Chr)
			return Chr>rhs.Chr;
		else if(StartPos!=rhs.StartPos)
			return StartPos>rhs.StartPos;
		else if(EndPos!=rhs.EndPos)
			return EndPos>rhs.EndPos;
		else
			return IsLeft>rhs.IsLeft;
	};

	bool operator == (const Breakpoint_t& rhs) const {
		return (Chr==rhs.Chr && StartPos==rhs.StartPos && EndPos==rhs.EndPos && IsLeft==rhs.IsLeft);
	};

	bool operator != (const Breakpoint_t& rhs) const {
		return !((*this)==rhs);
	};
};

class SV_t{
public:
	Breakpoint_t BP1;
	Breakpoint_t BP2;
public:
	SV_t(){};
	SV_t(Breakpoint_t _BP1, Breakpoint_t _BP2){
		if(_BP1<_BP2){
			BP1=_BP1;
			BP2=_BP2;
		}
		else{
			BP1=_BP2;
			BP2=_BP1;
		}
	};

	bool operator < (const SV_t& rhs) const {
		if(BP1!=rhs.BP1)
			return BP1<rhs.BP1;
		else
			return BP2<rhs.BP2;
	};

	bool operator == (const SV_t& rhs) const {
		return (BP1==rhs.BP1 && BP2==rhs.BP2);
	};
};

void ReadBEDPE(string BEDPEfile, map<string, int>& RefTable, vector<SV_t>& SVs)
{
	SVs.clear();
	ifstream input(BEDPEfile);
	string line;
	while(getline(input, line)){
		if(line[0]=='#')
			continue;
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t"));
		// remove sequence from mitochondra and contigs
		if(strs[0][0]=='M' || strs[0][0]=='G' || strs[0][0]=='K' || strs[3][0]=='M' || strs[3][0]=='G' || strs[3][0]=='K')
			continue;
		Breakpoint_t bp1(RefTable[strs[0]], stoi(strs[1]), stoi(strs[2]), (strs[8]=="-"));
		Breakpoint_t bp2(RefTable[strs[3]], stoi(strs[4]), stoi(strs[5]), (strs[9]=="-"));
		SV_t sv(bp1, bp2);
		SVs.push_back(sv);
	}
	input.close();
};

int SVfromAlignment(ReadRec_t r, vector<SV_t>& tmpSVs)
{
	int flag=-1;
	if(r.FirstRead.size()>0){
		for(int i=0; i<r.FirstRead.size()-1; i++){
			bool is_discordant=false;
			if(r.FirstRead[i].RefID!=r.FirstRead[i+1].RefID || r.FirstRead[i].IsReverse!=r.FirstRead[i+1].IsReverse)
				is_discordant=true;
			else if(!r.FirstRead[i].IsReverse && (r.FirstRead[i].RefPos<r.FirstRead[i+1].RefPos)!=(r.FirstRead[i].ReadPos<r.FirstRead[i+1].ReadPos))
				is_discordant=true;
			else if(r.FirstRead[i].IsReverse && (r.FirstRead[i].RefPos<r.FirstRead[i+1].RefPos)==(r.FirstRead[i].ReadPos<r.FirstRead[i+1].ReadPos))
				is_discordant=true;
			if(is_discordant){
				Breakpoint_t bp1(r.FirstRead[i].RefID, r.FirstRead[i].RefPos, r.FirstRead[i].RefPos+r.FirstRead[i].MatchRef, r.FirstRead[i].IsReverse);
				Breakpoint_t bp2(r.FirstRead[i+1].RefID, r.FirstRead[i+1].RefPos, r.FirstRead[i+1].RefPos+r.FirstRead[i+1].MatchRef, !r.FirstRead[i+1].IsReverse);
				SV_t tmpsv(bp1, bp2);
				tmpSVs.push_back(tmpsv);
				flag=0;
			}
		}
	}
	if(r.SecondMate.size()>0){
		for(int i=0; i<r.SecondMate.size()-1; i++){
			bool is_discordant=false;
			if(r.SecondMate[i].RefID!=r.SecondMate[i+1].RefID || r.SecondMate[i].IsReverse!=r.SecondMate[i+1].IsReverse)
				is_discordant=true;
			else if(!r.SecondMate[i].IsReverse && (r.SecondMate[i].RefPos<r.SecondMate[i+1].RefPos)!=(r.SecondMate[i].ReadPos<r.SecondMate[i+1].ReadPos))
				is_discordant=true;
			else if(r.SecondMate[i].IsReverse && (r.SecondMate[i].RefPos<r.SecondMate[i+1].RefPos)==(r.SecondMate[i].ReadPos<r.SecondMate[i+1].ReadPos))
				is_discordant=true;
			if(is_discordant){
				Breakpoint_t bp1(r.SecondMate[i].RefID, r.SecondMate[i].RefPos, r.SecondMate[i].RefPos+r.SecondMate[i].MatchRef, r.SecondMate[i].IsReverse);
				Breakpoint_t bp2(r.SecondMate[i+1].RefID, r.SecondMate[i+1].RefPos, r.SecondMate[i+1].RefPos+r.SecondMate[i+1].MatchRef, !r.SecondMate[i+1].IsReverse);
				SV_t tmpsv(bp1, bp2);
				tmpSVs.push_back(tmpsv);
				flag=0;
			}
		}
	}
	if(flag<0 && r.FirstRead.size()>0 && r.SecondMate.size()>0){ // only when each paired end read is concordant, consider pair junction
		if(r.IsPairDiscordant(false)){
			bool partial = false;
			if(r.FirstRead.size()!=0 && r.FirstRead.front().ReadPos > 12 && !r.FirstLowPhred)
				partial = true;
			if(r.FirstRead.size()!=0 && r.FirstTotalLen-r.FirstRead.back().ReadPos-r.FirstRead.back().MatchRead > 12 && !r.FirstLowPhred)
				partial = true;
			if(r.SecondMate.size()!=0 && r.SecondMate.front().ReadPos > 12 && !r.SecondLowPhred)
				partial = true;
			if(r.SecondMate.size()!=0 && r.SecondTotalLen-r.SecondMate.back().ReadPos-r.SecondMate.back().MatchRead > 12 && !r.SecondLowPhred)
				partial = true;
			if(partial){
				Breakpoint_t bp1(r.FirstRead.back().RefID, r.FirstRead.back().RefPos, r.FirstRead.back().RefPos+r.FirstRead.back().MatchRef, r.FirstRead.back().IsReverse);
				Breakpoint_t bp2(r.SecondMate.back().RefID, r.SecondMate.back().RefPos, r.SecondMate.back().RefPos+r.SecondMate.back().MatchRef, r.SecondMate.back().IsReverse);
				SV_t tmpsv(bp1, bp2);
				tmpSVs.push_back(tmpsv);
				flag=0;
			}
		}
	}
	return flag;
};

int NearestSV(SV_t newsv, vector<SV_t>& SVs, int thresh1=5, int thresh2=300)
{
	int optarg=-1;
	int optvalue=std::numeric_limits<int>::max();
	for(int i=0; i<SVs.size(); i++){
		if(newsv.BP1.Chr != SVs[i].BP1.Chr || newsv.BP2.Chr != SVs[i].BP2.Chr)
			continue;
		if(newsv.BP1.IsLeft != SVs[i].BP1.IsLeft || newsv.BP2.IsLeft != SVs[i].BP2.IsLeft)
			continue;
		if(newsv.BP1.IsLeft && (newsv.BP1.StartPos<SVs[i].BP1.StartPos-thresh1 || newsv.BP1.StartPos>SVs[i].BP1.StartPos+thresh2))
			continue;
		else if(!newsv.BP1.IsLeft && (newsv.BP1.EndPos<SVs[i].BP1.EndPos-thresh2 || newsv.BP1.EndPos>SVs[i].BP1.EndPos+thresh1))
			continue;
		if(newsv.BP2.IsLeft && (newsv.BP2.StartPos<SVs[i].BP2.StartPos-thresh1 || newsv.BP2.StartPos>SVs[i].BP2.StartPos+thresh2))
			continue;
		else if(!newsv.BP2.IsLeft && (newsv.BP2.EndPos<SVs[i].BP2.EndPos-thresh2 || newsv.BP2.EndPos>SVs[i].BP2.EndPos+thresh1))
			continue;
		int deviation=0;
		if(newsv.BP1.IsLeft)
			deviation+=abs(newsv.BP1.StartPos-SVs[i].BP1.StartPos);
		else
			deviation+=abs(newsv.BP1.EndPos-SVs[i].BP1.EndPos);
		if(newsv.BP2.IsLeft)
			deviation+=abs(newsv.BP2.StartPos-SVs[i].BP2.StartPos);
		else
			deviation+=abs(newsv.BP2.EndPos-SVs[i].BP2.EndPos);
		if(deviation<optvalue){
			optvalue=deviation;
			optarg=i;
		}
	}
	return optarg;
};

void FindReadSupport(SBamrecord_t Chimrecord, vector<SV_t>& SVs, vector< vector<SV_t> >& ReadSVs)
{
	ReadSVs.clear();
	for(int i=0; i<SVs.size(); i++){
		vector<SV_t> tmp;
		ReadSVs.push_back(tmp);
	}
	for(SBamrecord_t::iterator it=Chimrecord.begin(); it!=Chimrecord.end(); it++){
		vector<SV_t> tmp;
		int flag=SVfromAlignment(*it, tmp);
		if(flag!=-1){
			for(int i=0; i<tmp.size(); i++){
				int ind=NearestSV(tmp[i], SVs);
				if(ind!=-1)
					ReadSVs[ind].push_back(tmp[i]);
			}
		}
	}
};

vector<bool> ExactSequence(vector<SV_t>& SVs, vector< vector<SV_t> >& ReadSVs, vector< pair<int,int> >& NumSupports, vector< vector<SV_t> >& AltSVs, int thresh=5) // if an SV has split read support, modify the SV to only contain the exact region
{
	NumSupports.clear();
	AltSVs.clear();
	vector<bool> flags(SVs.size(), false);
	for(int i=0; i<ReadSVs.size(); i++){
		NumSupports.push_back( make_pair(0,0) );
		const vector<SV_t>& readsvs=ReadSVs[i];
		SV_t& sv=SVs[i];
		if(readsvs.size()==0)
			continue;
		// sort each breakpoint of the SV
		vector<Breakpoint_t> BP1s;
		vector<Breakpoint_t> BP2s;
		for(vector<SV_t>::const_iterator it=readsvs.cbegin(); it!=readsvs.cend(); it++){
			BP1s.push_back(it->BP1);
			BP2s.push_back(it->BP2);
		}
		if(sv.BP1.IsLeft)
			sort(BP1s.begin(), BP1s.end());
		else
			sort(BP1s.begin(), BP1s.end(), [](Breakpoint_t a, Breakpoint_t b){if(a.Chr!=b.Chr) return a.Chr<b.Chr; else if(a.EndPos!=b.EndPos) return a.EndPos<b.EndPos; 
				else if(a.StartPos!=b.StartPos) return a.StartPos<b.StartPos; else return a.IsLeft<b.IsLeft;} );
		if(sv.BP2.IsLeft)
			sort(BP2s.begin(), BP2s.end());
		else
			sort(BP2s.begin(), BP2s.end(), [](Breakpoint_t a, Breakpoint_t b){if(a.Chr!=b.Chr) return a.Chr<b.Chr; else if(a.EndPos!=b.EndPos) return a.EndPos<b.EndPos; 
				else if(a.StartPos!=b.StartPos) return a.StartPos<b.StartPos; else return a.IsLeft<b.IsLeft;} );
		// check whether the breakpoints hit predicted breakpoints in thresh window
		bool hitbp1=false, hitbp2=false;
		for(vector<Breakpoint_t>::iterator it=BP1s.begin(); it!=BP1s.end(); it++){
			if((sv.BP1.IsLeft && abs(sv.BP1.StartPos-it->StartPos)<thresh) || (!sv.BP1.IsLeft && abs(sv.BP1.EndPos-it->EndPos)<thresh))
				hitbp1=true;
		}
		for(vector<Breakpoint_t>::iterator it=BP2s.begin(); it!=BP2s.end(); it++){
			if((sv.BP2.IsLeft && abs(sv.BP2.StartPos-it->StartPos)<thresh) || (!sv.BP2.IsLeft && abs(sv.BP2.EndPos-it->EndPos)<thresh))
				hitbp2=true;
		}
		if(!hitbp1 || !hitbp2)
			continue;
		// if both breakpoints has split read hitting them, extend from split read
		// cout<<"old\t"<<(sv.BP1.Chr)<<"\t"<<(sv.BP1.StartPos)<<"\t"<<(sv.BP1.EndPos)<<"\t"<<(sv.BP1.IsLeft)<<"\t";
		// cout<<(sv.BP2.Chr)<<"\t"<<(sv.BP2.StartPos)<<"\t"<<(sv.BP2.EndPos)<<"\t"<<(sv.BP2.IsLeft)<<endl;
		bool flag_bp1 = false;
		bool flag_bp2 = false;
		int support_bp1 = 0;
		int support_bp2 = 0;
		if(sv.BP1.IsLeft){
			vector<Breakpoint_t>::iterator it=BP1s.begin();
			while(abs(it->StartPos-sv.BP1.StartPos)>=thresh && it!=BP1s.end())
				it++;
			assert(it!=BP1s.end());
			int rightmost=it->EndPos;
			for(; it!=BP1s.end(); it++)
				if(it->StartPos < rightmost)
					rightmost=(it->EndPos>rightmost)?(it->EndPos):rightmost;
			if(sv.BP1.StartPos<rightmost){
				sv.BP1.EndPos=(rightmost<sv.BP1.EndPos)?rightmost:(sv.BP1.EndPos);
				flag_bp1 = true;
				support_bp1++;
			}
		}
		else{
			vector<Breakpoint_t>::reverse_iterator it=BP1s.rbegin();
			while(abs(it->EndPos-sv.BP1.EndPos)>=thresh && it!=BP1s.rend())
				it++;
			assert(it!=BP1s.rend());
			int leftmost=it->StartPos;
			for(; it!=BP1s.rend(); it++)
				if(it->EndPos > leftmost)
					leftmost=(it->StartPos<leftmost)?(it->StartPos):leftmost;
			if(leftmost < sv.BP1.EndPos){
				sv.BP1.StartPos=(leftmost>sv.BP1.StartPos)?leftmost:(sv.BP1.StartPos);
				flag_bp1 = true;
				support_bp1++;
			}
		}
		if(sv.BP2.IsLeft){
			vector<Breakpoint_t>::iterator it=BP2s.begin();
			while(abs(it->StartPos-sv.BP2.StartPos)>=thresh && it!=BP2s.end())
				it++;
			assert(it!=BP2s.end());
			int rightmost=it->EndPos;
			for(vector<Breakpoint_t>::iterator it=BP2s.begin(); it!=BP2s.end(); it++)
				if(it->StartPos < rightmost)
					rightmost=(it->EndPos>rightmost)?(it->EndPos):rightmost;
			if(sv.BP2.StartPos < rightmost){
				sv.BP2.EndPos=(rightmost<sv.BP2.EndPos)?rightmost:(sv.BP2.EndPos);
				flag_bp2 = true;
				support_bp2++;
			}
		}
		else{
			vector<Breakpoint_t>::reverse_iterator it=BP2s.rbegin();
			while(abs(it->EndPos-sv.BP2.EndPos)>=thresh && it!=BP2s.rend())
				it++;
			assert(it!=BP2s.rend());
			int leftmost=it->StartPos;
			for(vector<Breakpoint_t>::reverse_iterator it=BP2s.rbegin(); it!=BP2s.rend(); it++)
				if(it->EndPos > leftmost)
					leftmost=(it->StartPos<leftmost)?(it->StartPos):leftmost;
			if(leftmost < sv.BP2.EndPos){
				sv.BP2.StartPos=(leftmost>sv.BP2.StartPos)?leftmost:(sv.BP2.StartPos);
				flag_bp2 = true;
				support_bp2++;
			}
		}
		// record whether have spanning read information
		if(flag_bp1 && flag_bp2){
			flags[i] = true;
			NumSupports.back() = make_pair(support_bp1, support_bp2);
		}
		
		// find alternative junction points
		vector<SV_t> alternativeSVs, tmpalternativeSVs;
		for(vector<SV_t>::const_iterator it=readsvs.cbegin(); it!=readsvs.cend(); it++){
			SV_t altsv(sv.BP1, sv.BP2);
			bool hasalt_bp1 = false;
			bool hasalt_bp2 = false;
			if(sv.BP1.IsLeft==it->BP1.IsLeft){
				if(sv.BP1.IsLeft && abs(sv.BP1.StartPos - it->BP1.StartPos)<thresh){
					altsv.BP1.StartPos = it->BP1.StartPos;
					if(sv.BP1.StartPos != it->BP2.StartPos)
						hasalt_bp1 = true;
				}
				else if(!sv.BP1.IsLeft && abs(sv.BP1.EndPos - it->BP1.EndPos)<thresh){
					altsv.BP1.EndPos = it->BP1.EndPos;
					if(sv.BP1.EndPos != it->BP1.EndPos)
						hasalt_bp1 = true;
				}
			}
			if(sv.BP2.IsLeft == it->BP2.IsLeft){
				if(sv.BP2.IsLeft && abs(sv.BP2.StartPos - it->BP2.StartPos)<thresh){
					altsv.BP2.StartPos = it->BP2.StartPos;
					if(sv.BP2.StartPos != it->BP2.StartPos)
						hasalt_bp2 = true;
				}
				else if(!sv.BP2.IsLeft && abs(sv.BP2.EndPos - it->BP2.EndPos)<thresh){
					altsv.BP2.EndPos = it->BP2.EndPos;
					if(sv.BP2.EndPos != it->BP2.EndPos)
						hasalt_bp2 = true;
				}
			}
			if(hasalt_bp1 && hasalt_bp2)
				tmpalternativeSVs.push_back(altsv);
		}
		// if there is alternative junction point, then the original junction point must be supported by some spanning reads
		if(tmpalternativeSVs.size()!=0)
			assert(flags[i]);
		sort(tmpalternativeSVs.begin(), tmpalternativeSVs.end());
		for(vector<SV_t>::iterator it=tmpalternativeSVs.begin(); it!=tmpalternativeSVs.end(); it++){
			if(alternativeSVs.size()==0 || !(alternativeSVs.back()==(*it)))
				alternativeSVs.push_back(*it);
		}
		AltSVs.push_back(alternativeSVs);
	}
	return flags;
};

void ReadGenome(string FAfile, vector<string>& Genome, map<string,int>& RefTable)
{
	Genome.clear();
	Genome.resize(RefTable.size());
	ifstream input(FAfile);
	string line;
	string prevname="";
	string preseq="";
	while(getline(input, line)){
		if(line[0]=='>'){
			if(prevname!="")
				Genome[RefTable[prevname]]=preseq;
			vector<string> strs;
			boost::split(strs, line, boost::is_any_of(" \t"));
			prevname=strs[0].substr(1);
			preseq="";
		}
		else
			preseq+=line;
	}
	Genome[RefTable[prevname]]=preseq;
	input.close();
};

void WritePreciseJunction(string outfile, vector<SV_t>& SVs, vector<bool>& flags, vector<string>& Genome, vector<string>& RefName)
{
	ofstream output(outfile, ios::out);
	assert(SVs.size()==flags.size());
	for(int i=0; i<SVs.size(); i++){
		if(!flags[i])
			continue;
		const SV_t& sv=SVs[i];
		assert(Genome[sv.BP1.Chr].size()>=sv.BP1.EndPos && Genome[sv.BP2.Chr].size()>=sv.BP2.EndPos);
		string seq1=Genome[sv.BP1.Chr].substr(sv.BP1.StartPos, sv.BP1.EndPos-sv.BP1.StartPos);
		string seq2=Genome[sv.BP2.Chr].substr(sv.BP2.StartPos, sv.BP2.EndPos-sv.BP2.StartPos);
		if(sv.BP1.IsLeft)
			ReverseComplement(seq1.begin(), seq1.end());
		if(!sv.BP2.IsLeft)
			ReverseComplement(seq2.begin(), seq2.end());
		string seq=seq1+seq2;
		output<<">squid_"<<i<<" "<<RefName[sv.BP1.Chr]<<":"<<(sv.BP1.StartPos)<<":"<<(sv.BP1.EndPos)<<":"<<(sv.BP1.IsLeft?"-":"+");
		output<<" "<<RefName[sv.BP2.Chr]<<":"<<(sv.BP2.StartPos)<<":"<<(sv.BP2.EndPos)<<":"<<(sv.BP2.IsLeft?"+":"-")<<endl;
		int count=0;
		while(count<(int)seq.size()){
			output<<(seq.substr(count, min(80,(int)seq.size()-count)))<<endl;
			count+=80;
		}
	}
	output.close();
};

void WriteRelaxedJunction(string outfile, vector<SV_t>& SVs, vector<bool>& flags, vector<string>& Genome, vector<string>& RefName)
{
	ofstream output(outfile, ios::out);
	assert(SVs.size()==flags.size());
	for(int i=0; i<SVs.size(); i++){
		const SV_t& sv=SVs[i];
		SV_t tmp=sv;
		assert(Genome[sv.BP1.Chr].size()>=sv.BP1.EndPos && Genome[sv.BP2.Chr].size()>=sv.BP2.EndPos);
		if(flags[i]){
			if(sv.BP1.IsLeft)
				tmp.BP1.EndPos=min(tmp.BP1.EndPos+1000, (int)Genome[sv.BP1.Chr].size());
			else
				tmp.BP1.StartPos=max(0, tmp.BP1.StartPos-1000);
			if(sv.BP2.IsLeft)
				tmp.BP2.EndPos=min(tmp.BP2.EndPos+1000, (int)Genome[sv.BP2.Chr].size());
			else
				tmp.BP2.StartPos=max(0, tmp.BP2.StartPos-1000);
		}
		string seq1=Genome[tmp.BP1.Chr].substr(tmp.BP1.StartPos, tmp.BP1.EndPos-tmp.BP1.StartPos);
		string seq2=Genome[tmp.BP2.Chr].substr(tmp.BP2.StartPos, tmp.BP2.EndPos-tmp.BP2.StartPos);
		if(tmp.BP1.IsLeft)
			ReverseComplement(seq1.begin(), seq1.end());
		if(!tmp.BP2.IsLeft)
			ReverseComplement(seq2.begin(), seq2.end());
		string seq=seq1+seq2;
		output<<">squid_"<<i<<" "<<RefName[tmp.BP1.Chr]<<":"<<(tmp.BP1.StartPos)<<":"<<(tmp.BP1.EndPos)<<":"<<(tmp.BP1.IsLeft?"-":"+");
		output<<" "<<RefName[tmp.BP2.Chr]<<":"<<(tmp.BP2.StartPos)<<":"<<(tmp.BP2.EndPos)<<":"<<(tmp.BP2.IsLeft?"+":"-")<<endl;
		int count=0;
		while(count<(int)seq.size()){
			output<<(seq.substr(count, min(80,(int)seq.size()-count)))<<endl;
			count+=80;
		}
	}
	output.close();
};

void WriteAlternativeJunction(string outfile, vector< vector<SV_t> >& AltSVs, vector<string>& Genome, vector<string>& RefName)
{
	ofstream output(outfile, ios::out);
	for(int i=0; i<AltSVs.size(); i++){
		for(int j=0; j<AltSVs[i].size(); j++){
			const SV_t& sv=AltSVs[i][j];
			assert(Genome[sv.BP1.Chr].size()>=sv.BP1.EndPos && Genome[sv.BP2.Chr].size()>=sv.BP2.EndPos);
			string seq1=Genome[sv.BP1.Chr].substr(sv.BP1.StartPos, sv.BP1.EndPos-sv.BP1.StartPos);
			string seq2=Genome[sv.BP2.Chr].substr(sv.BP2.StartPos, sv.BP2.EndPos-sv.BP2.StartPos);
			if(sv.BP1.IsLeft)
				ReverseComplement(seq1.begin(), seq1.end());
			if(!sv.BP2.IsLeft)
				ReverseComplement(seq2.begin(), seq2.end());
			string seq=seq1+seq2;

			output<<">squid_"<<i<<"_alt_"<<(j+1)<<" "<<RefName[sv.BP1.Chr]<<":"<<(sv.BP1.StartPos)<<":"<<(sv.BP1.EndPos)<<":"<<(sv.BP1.IsLeft?"-":"+");
			output<<" "<<RefName[sv.BP2.Chr]<<":"<<(sv.BP2.StartPos)<<":"<<(sv.BP2.EndPos)<<":"<<(sv.BP2.IsLeft?"+":"-")<<endl;
			int count=0;
			while(count<(int)seq.size()){
				output<<(seq.substr(count, min(80,(int)seq.size()-count)))<<endl;
				count+=80;
			}
		}
	}
	output.close();
};


uint16_t ReadLen;
bool Phred_Type = 1; //1 for phred33, 0 for phred64
uint16_t Max_LowPhred_Len = 10;
uint8_t Min_Phred = 4;
uint16_t Min_MapQual = 1;

int main(int argc, char* argv[]){
	if(argc==1){
		printf("junctionsequence <BEDPEfile> <Input_Chim_BAM> <FA_genome> <OUTPrefix>\n");
	}
	else{
		string BEDPEfile(argv[1]);
		string Input_Chim_BAM(argv[2]);
		string FAfile(argv[3]);
		string OUTPrefix(argv[4]);

		map<string, int> RefTable;
		vector<string> RefName;
		vector<int> RefLength;
		SBamrecord_t Chimrecord;
		BuildRefName(Input_Chim_BAM, RefName, RefTable, RefLength);
		BuildChimericSBamRecord(Chimrecord, RefName, Input_Chim_BAM);

		vector<SV_t> SVs;
		ReadBEDPE(BEDPEfile, RefTable, SVs);

		vector< vector<SV_t> > ReadSVs;
		FindReadSupport(Chimrecord, SVs, ReadSVs);
		vector< pair<int,int> > NumSupports;
		vector< vector<SV_t> > AltSVs;
		vector<bool> flags=ExactSequence(SVs, ReadSVs, NumSupports, AltSVs);

		vector<string> Genome;
		ReadGenome(FAfile, Genome, RefTable);
		WritePreciseJunction(OUTPrefix+"_junc_precise.fa", SVs, flags, Genome, RefName);
		WriteRelaxedJunction(OUTPrefix+"_junc_relax.fa", SVs, flags, Genome, RefName);
	}
}
