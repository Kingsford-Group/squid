/*
Part of SQUID transcriptomic structural variation detector
(c) 2017 by  Cong Ma, Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "ReadRec.h"
#include "SegmentGraph.h"

ReadRec_t::ReadRec_t(BamAlignment record){
	Qname=record.Name;
	if(Qname.substr(Qname.size()-2)=="/1" || Qname.substr(Qname.size()-2)=="/2")
		Qname=Qname.substr(0, Qname.size()-2);
	MultiFilter=false;
	int32_t ReadPos=0, RefPos=record.Position, TotalLen=0, LowPhredLen=0, tmpLowPhredLen=0;
	for(vector<CigarOp>::const_iterator itcigar=record.CigarData.begin(); itcigar!=record.CigarData.end(); itcigar++)
		if(itcigar->Type=='M' || itcigar->Type=='S' || itcigar->Type=='H' || itcigar->Type=='I' || itcigar->Type=='=' || itcigar->Type=='X')
			TotalLen+=itcigar->Length;
	if(Phred_Type){
		for(unsigned int i=0; i<record.Qualities.size(); i++){
			if(record.Qualities[i]<(char)(33+Min_Phred))
				tmpLowPhredLen++;
			else
				tmpLowPhredLen=0;
			if(LowPhredLen<tmpLowPhredLen)
				LowPhredLen=tmpLowPhredLen;
		}
	}
	else{
		for(unsigned int i=0; i<record.Qualities.size(); i++){
			if(record.Qualities[i]<(char)(64+Min_Phred))
				tmpLowPhredLen++;
			else
				tmpLowPhredLen=0;
			if(LowPhredLen<tmpLowPhredLen)
				LowPhredLen=tmpLowPhredLen;
		}
	}
	if(record.IsFirstMate()){
		FirstTotalLen=TotalLen; SecondTotalLen=0; FirstLowPhred=(LowPhredLen>Max_LowPhred_Len);
	}
	else{
		SecondTotalLen=TotalLen; FirstTotalLen=0; SecondLowPhred=(LowPhredLen>Max_LowPhred_Len);
	}
	for(vector<CigarOp>::const_iterator itcigar=record.CigarData.begin(); itcigar!=record.CigarData.end(); itcigar++){
		if(itcigar->Type=='S' || itcigar->Type=='H'){
			ReadPos+=itcigar->Length;
		}
		else if(itcigar->Type=='M' || itcigar->Type=='='){
			int tmpRead=0, tmpRef=0;
			vector<CigarOp>::const_iterator itcigar2;
			for(itcigar2=itcigar; itcigar2!=record.CigarData.end() && itcigar2->Type!='S' && itcigar2->Type!='H' && itcigar2->Type!='N'; itcigar2++){
				if(itcigar2->Type!='D')
					tmpRead+=itcigar2->Length;
				if(itcigar2->Type!='I')
					tmpRef+=itcigar2->Length;
			}
			// the following line calculate the ratio of poly A and T. If the aligned block is more than 75% A and T, consider it as polyA/T tail, not actual sequence
			int polyAcount=0;
			int polyTcount=0;
			for(int i=ReadPos; i<ReadPos+tmpRead; i++){
				if(record.QueryBases[i]=='a' || record.QueryBases[i]=='A')
					polyAcount++;
				else if(record.QueryBases[i]=='t' || record.QueryBases[i]=='T')
					polyTcount++;
			}
			// add this aligned block only if it is not poly A/T
			if(1.0*polyAcount/tmpRead<0.75 && 1.0*polyTcount/tmpRead<0.75){
				SingleBamRec_t tmp(record.RefID, RefPos, ReadPos, tmpRef, tmpRead, record.MapQuality, record.IsReverseStrand(), record.IsFirstMate());
				if(record.IsReverseStrand())
					tmp.ReadPos=TotalLen-ReadPos-tmpRead;
				if(record.IsFirstMate())
					FirstRead.push_back(tmp);
				else
					SecondMate.push_back(tmp);
			}
			ReadPos+=tmpRead;
			RefPos+=tmpRef;
			itcigar=itcigar2; itcigar--;
		}
		else if(itcigar->Type=='N')
			RefPos+=itcigar->Length;
	}
};

bool ReadRec_t::FrontSmallerThan(const ReadRec_t& lhs, const ReadRec_t& rhs){
	if(lhs.FirstRead.size()!=0 && rhs.FirstRead.size()!=0){
		if(lhs.FirstRead.front().RefID!=rhs.FirstRead.front().RefID)
			return lhs.FirstRead.front().RefID<rhs.FirstRead.front().RefID;
		else
			return lhs.FirstRead.front().RefPos<rhs.FirstRead.front().RefPos;
	}
	else if(lhs.SecondMate.size()!=0 && rhs.SecondMate.size()!=0){
		if(lhs.SecondMate.front().RefID!=rhs.SecondMate.front().RefID)
			return lhs.SecondMate.front().RefID<rhs.SecondMate.front().RefID;
		else
			return lhs.SecondMate.front().RefPos<rhs.SecondMate.front().RefPos;
	}
	else if(lhs.FirstRead.size()!=0 && rhs.SecondMate.size()!=0){
		if(lhs.FirstRead.front().RefID!=rhs.SecondMate.front().RefID)
			return lhs.FirstRead.front().RefID<rhs.SecondMate.front().RefID;
		else
			return lhs.FirstRead.front().RefPos<rhs.SecondMate.front().RefPos;
	}
	else if(lhs.SecondMate.size()!=0 && rhs.FirstRead.size()!=0){
		if(lhs.SecondMate.front().RefID!=rhs.FirstRead.front().RefID)
			return lhs.SecondMate.front().RefID<rhs.FirstRead.front().RefID;
		else
			return lhs.SecondMate.front().RefPos<rhs.FirstRead.front().RefPos;
	}
	else
		return false;
};

bool ReadRec_t::Equal(const ReadRec_t& lhs, const ReadRec_t& rhs){
	bool same1=false;
	bool same2=false;
	if(lhs.FirstRead.size()==rhs.FirstRead.size() && lhs.SecondMate.size()==rhs.SecondMate.size()){
		same1=true;
		for(unsigned int i=0; i<lhs.FirstRead.size(); i++)
			if(lhs.FirstRead[i].RefID!=rhs.FirstRead[i].RefID || lhs.FirstRead[i].RefPos!=rhs.FirstRead[i].RefPos || lhs.FirstRead[i].MatchRef!=rhs.FirstRead[i].MatchRef)
				same1=false;
		for(unsigned int i=0; i<lhs.SecondMate.size(); i++)
			if(lhs.SecondMate[i].RefID!=rhs.SecondMate[i].RefID || lhs.SecondMate[i].RefPos!=rhs.SecondMate[i].RefPos || lhs.SecondMate[i].MatchRef!=rhs.SecondMate[i].MatchRef)
				same1=false;
	}
	if(lhs.FirstRead.size()==rhs.SecondMate.size() && lhs.SecondMate.size()==rhs.FirstRead.size()){
		same2=true;
		for(unsigned int i=0; i<lhs.FirstRead.size(); i++)
			if(lhs.FirstRead[i].RefID!=rhs.SecondMate[i].RefID || lhs.FirstRead[i].RefPos!=rhs.SecondMate[i].RefPos || lhs.FirstRead[i].MatchRef!=rhs.SecondMate[i].MatchRef)
				same2=false;
		for(unsigned int i=0; i<lhs.SecondMate.size(); i++)
			if(lhs.SecondMate[i].RefID!=rhs.FirstRead[i].RefID || lhs.SecondMate[i].RefPos!=rhs.FirstRead[i].RefPos || lhs.SecondMate[i].MatchRef!=rhs.FirstRead[i].MatchRef)
				same2=false;
	}
	return (same1 || same2);
};

void ReadRec_t::SortbyReadPos(){
	sort(FirstRead.begin(), FirstRead.end(), SingleBamRec_t::CompReadPos);
	sort(SecondMate.begin(), SecondMate.end(), SingleBamRec_t::CompReadPos);
};

void ReadRec_t::FilterSplitRecord(){
	for(int i=0; i<(int)FirstRead.size()-1; i++){
		if(FirstRead[i].ReadPos+FirstRead[i].MatchRead-FirstRead[i+1].ReadPos > 10){
			if(FirstRead[i].MapQual > FirstRead[i+1].MapQual){
				FirstRead.erase(FirstRead.begin()+i+1); i--;
			}
			else if(FirstRead[i].MapQual < FirstRead[i+1].MapQual){
				FirstRead.erase(FirstRead.begin()+i); i--;
			}
		}
	}
	for(int i=0; i<(int)SecondMate.size()-1; i++){
		if(SecondMate[i].ReadPos+SecondMate[i].MatchRead-SecondMate[i].ReadPos > 10){
			if(SecondMate[i].MapQual > SecondMate[i+1].MapQual){
				SecondMate.erase(SecondMate.begin()+i+1); i--;
			}
			else if(SecondMate[i].MapQual < SecondMate[i+1].MapQual){
				SecondMate.erase(SecondMate.begin()+i); i--;
			}
		}
	}
};

bool ReadRec_t::IsSingleAnchored() const{
	if((FirstRead.size()==0 || SecondMate.size()==0) && !MultiFilter)
		return true;
	else
		return false;
};

bool ReadRec_t::IsEndDiscordant(bool _isfirst) const{
	if(_isfirst){
		if(FirstRead.size()<=1)
			return false;
		else{
			for(unsigned int i=0; i<FirstRead.size()-1; i++){
				if(FirstRead[i].RefID!=FirstRead[i+1].RefID || FirstRead[i].IsReverse!=FirstRead[i+1].IsReverse)
					return true;
				else if(!FirstRead[i].IsReverse && (FirstRead[i].RefPos<FirstRead[i+1].RefPos)!=(FirstRead[i].ReadPos<FirstRead[i+1].ReadPos))
					return true;
				else if(FirstRead[i].IsReverse && (FirstRead[i].RefPos<FirstRead[i+1].RefPos)==(FirstRead[i].ReadPos<FirstRead[i+1].ReadPos))
					return true;
			}
			return false;
		}
	}
	else{
		if(SecondMate.size()<=1)
			return false;
		else{
			for(unsigned int i=0; i<SecondMate.size()-1; i++){
				if(SecondMate[i].RefID!=SecondMate[i+1].RefID || SecondMate[i].IsReverse!=SecondMate[i+1].IsReverse)
					return true;
				else if(!SecondMate[i].IsReverse && (SecondMate[i].RefPos<SecondMate[i+1].RefPos)!=(SecondMate[i].ReadPos<SecondMate[i+1].ReadPos))
					return true;
				else if(SecondMate[i].IsReverse && (SecondMate[i].RefPos<SecondMate[i+1].RefPos)==(SecondMate[i].ReadPos<SecondMate[i+1].ReadPos))
					return true;
			}
			return false;
		}
	}
};

bool ReadRec_t::IsPairDiscordant(bool needcheck) const{
	if(FirstRead.size()==0 || SecondMate.size()==0)
		return false;
	else{
		if(needcheck){
			if(IsEndDiscordant(true) || IsEndDiscordant(false))
				return true;
		}
		if(FirstRead.front().RefID!=SecondMate.back().RefID || FirstRead.front().IsReverse==SecondMate.back().IsReverse)
			return true;
		else if(!FirstRead.front().IsReverse && FirstRead.front().RefPos-FirstRead.front().ReadPos>SecondMate.back().RefPos-(SecondTotalLen-SecondMate.back().ReadPos-SecondMate.back().MatchRead))
			return true;
		else if(!SecondMate.front().IsReverse && SecondMate.front().RefPos-SecondMate.front().ReadPos>FirstRead.back().RefPos-(FirstTotalLen-FirstRead.back().ReadPos-FirstRead.back().MatchRead))
			return true;
		else
			return false;
	}
};

bool ReadRec_t::IsDiscordant() const{
	return (IsSingleAnchored() || IsEndDiscordant(true) || IsEndDiscordant(false) || IsPairDiscordant(false));
};

int ReadRec_t::ReadCoverageGap() const{
	int prevend=0, firstgap=0, secondgap=0;
	for(vector<SingleBamRec_t>::const_iterator it=FirstRead.cbegin(); it!=FirstRead.cend(); it++){
		if(it->ReadPos>prevend)
			firstgap+=it->ReadPos-prevend;
		prevend=it->ReadPos+it->MatchRead;
	}
	firstgap+=FirstTotalLen-prevend;
	prevend=0;
	for(vector<SingleBamRec_t>::const_iterator it=SecondMate.cbegin(); it!=SecondMate.cend(); it++){
		if(it->ReadPos>prevend)
			secondgap+=it->ReadPos-prevend;
		prevend=it->ReadPos+it->MatchRead;
	}
	secondgap+=SecondTotalLen-prevend;
	return firstgap+secondgap;
};

int ReadRec_t::ReadCoverageGap() {
	return (static_cast<const ReadRec_t*>(this)->ReadCoverageGap());
};

string ReadRec_t::Print(){
	string info=Qname+"\tFirst{";
	for(vector<SingleBamRec_t>::iterator it=FirstRead.begin(); it!=FirstRead.end(); it++)
		info+="("+to_string(it->RefID)+","+to_string(it->RefPos)+","+to_string(it->MatchRef)+","+to_string(it->ReadPos)+","+(it->IsReverse?"-":"+")+") ";
	info+="}\tSecond{";
	for(vector<SingleBamRec_t>::iterator it=SecondMate.begin(); it!=SecondMate.end(); it++)
		info+="("+to_string(it->RefID)+","+to_string(it->RefPos)+","+to_string(it->MatchRef)+","+to_string(it->ReadPos)+","+(it->IsReverse?"-":"+")+") ";
	info+="}";
	return info;
};

void BuildRefName(string bamfile, vector<string>& RefName, std::map<string,int> & RefTable, vector<int>& RefLength){
	RefName.clear();
	RefTable.clear();
	RefLength.clear();
	BamReader bamreader; bamreader.Open(bamfile);
	if(bamreader.IsOpen()){
		int count=0;
		SamHeader header=bamreader.GetHeader();
		for(SamSequenceIterator it=header.Sequences.Begin();it!=header.Sequences.End();it++){
			RefName.push_back(it->Name);
			RefTable[it->Name]=count++;
			RefLength.push_back(stoi(it->Length));
		}
	}
	else
		cout<<"Cannot open bamfile "<<bamfile<<endl;
};

bool BuildRefSeq(string fafile, const std::map<string,int>& RefTable, vector<int>& RefLength, vector<string>& RefSequence){
	ifstream input(fafile);
	RefSequence.resize(RefTable.size());
	for(unsigned int i=0; i<RefLength.size(); i++)
		RefSequence[i].reserve(RefLength[i]);
	string line;
	int ind=-1;
	while(getline(input,line)){
		if(line[0]=='>'){
			vector<string> strs;
			boost::split(strs, line, boost::is_any_of("\t "));
			map<string,int>::const_iterator cit=RefTable.find(strs[0].substr(1));
			if(cit!=RefTable.end())
				ind=cit->second;
			else
				ind=-1;
		}
		else{
			if(ind!=-1)
				RefSequence[ind]+=line;
		}
	}
	for(unsigned int i=0; i<RefSequence.size(); i++)
        if((int)RefSequence[i].size()!=RefLength[i]){
        	cout<<"FASTA file doesn't match BAM file"<<endl;
            RefSequence.clear();
            return false;
        }
    return true;
};

void UpdateReference(const SegmentGraph_t& SegmentGraph, const vector< vector<int> >& Components, vector<int>& RefLength, vector<string>& RefSequence){
	vector<string> newRefSequence(RefSequence.size());
	for(unsigned int i=0; i<Components.size(); i++){
		int length=0;
		for(unsigned int j=0; j<Components[i].size(); j++){
			length+=SegmentGraph.vNodes[abs(Components[i][j])-1].Length;
			newRefSequence[i]+=RefSequence[SegmentGraph.vNodes[abs(Components[i][j])-1].Chr].substr(SegmentGraph.vNodes[abs(Components[i][j])-1].Position, SegmentGraph.vNodes[abs(Components[i][j])-1].Length);
		}
		RefLength[i]=length;
	}
	RefSequence=newRefSequence;
};

void BuildChimericSBamRecord(SBamrecord_t& SBamrecord, const vector<string>& RefName, string bamfile){
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Start reading bam file."<<endl;

	vector<uint16_t> sample_ReadLen; sample_ReadLen.reserve(5);
	SBamrecord.clear();
	SBamrecord.reserve(66536);
	SBamrecord_t newSBamrecord;
	BamReader bamreader; bamreader.Open(bamfile);
	if(bamreader.IsOpen()){
		BamAlignment record;
		while(bamreader.GetNextAlignment(record)){
			if(record.IsMapped() && !record.IsDuplicate()){
				ReadRec_t tmp(record);
				SBamrecord.push_back(tmp);
				if(sample_ReadLen.size()<sample_ReadLen.capacity())
					sample_ReadLen.push_back(max(tmp.FirstTotalLen, tmp.SecondTotalLen));
				if(SBamrecord.capacity()==SBamrecord.size())
					SBamrecord.reserve(SBamrecord.size()*2);
			}
		}
		SBamrecord.reserve(SBamrecord.size());
		sort(SBamrecord.begin(), SBamrecord.end());
		newSBamrecord.reserve(SBamrecord.size());
		for(SBamrecord_t::iterator it=SBamrecord.begin();it!=SBamrecord.end();it++){
			if(newSBamrecord.size()==0 || it->Qname!=newSBamrecord.back().Qname)
				newSBamrecord.push_back((*it));
			else{
				if(newSBamrecord.back().FirstTotalLen==0 && it->FirstTotalLen!=0){
					newSBamrecord.back().FirstTotalLen=it->FirstTotalLen;
					newSBamrecord.back().FirstLowPhred=it->FirstLowPhred;
				}
				if(newSBamrecord.back().SecondTotalLen==0 && it->SecondTotalLen!=0){
					newSBamrecord.back().SecondTotalLen=it->SecondTotalLen;
					newSBamrecord.back().SecondLowPhred=it->SecondLowPhred;
				}
				for(vector<SingleBamRec_t>::iterator itsingle=it->FirstRead.begin(); itsingle!=it->FirstRead.end(); itsingle++)
					newSBamrecord.back().FirstRead.push_back(*itsingle);
				for(vector<SingleBamRec_t>::iterator itsingle=it->SecondMate.begin(); itsingle!=it->SecondMate.end(); itsingle++)
					newSBamrecord.back().SecondMate.push_back(*itsingle);
			}
		}
		newSBamrecord.reserve(newSBamrecord.size());
		for(SBamrecord_t::iterator it=newSBamrecord.begin(); it!=newSBamrecord.end(); it++)
			it->SortbyReadPos();
		// infer read length
		sort(sample_ReadLen.begin(), sample_ReadLen.end());
		ReadLen=sample_ReadLen[(int)sample_ReadLen.size()/2];
		bamreader.Close();
	}
	sort(newSBamrecord.begin(), newSBamrecord.end(), ReadRec_t::FrontSmallerThan);
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Finish sorting Chimeric bam reads."<<endl;

	// remove PCR duplicates
	SBamrecord.clear();
	for(SBamrecord_t::iterator it=newSBamrecord.begin(); it!=newSBamrecord.end(); it++){
		if(SBamrecord.size()==0)
			SBamrecord.push_back(*it);
		else if(it->FirstRead.size()==0 || SBamrecord.back().FirstRead.size()==0)
			SBamrecord.push_back((*it));
		else if(it->FirstRead.front().RefID!=SBamrecord.back().FirstRead.front().RefID || it->FirstRead.front().RefPos!=SBamrecord.back().FirstRead.front().RefPos)
			SBamrecord.push_back(*it);
		else{
			bool isdup=false;
			for(SBamrecord_t::reverse_iterator it2=SBamrecord.rbegin(); it2!=SBamrecord.rend(); it2++){
				if(it2->FirstRead.size()==0 || it->FirstRead.front().RefID!=it2->FirstRead.front().RefID || it->FirstRead.front().RefPos!=it2->FirstRead.front().RefPos)
					break;
				if(ReadRec_t::Equal(*it, *it2)){
					isdup=true;
					break;
				}
			}
			if(!isdup)
				SBamrecord.push_back(*it);
		}
	}
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Finish removing PCR duplicates."<<endl;
};

int AlignmentStat(const SBamrecord_t& SBamrecord){
	int numDiscordant=0;
	for(vector<ReadRec_t>::const_iterator it=SBamrecord.begin(); it!=SBamrecord.end(); it++){
		if(it->IsDiscordant())
			numDiscordant++;
	}
	return numDiscordant;
};

int AlignmentStat(const SBamrecord_t& SBamrecord, string outputfile){
	ofstream output(outputfile, ios::out);
	int numDiscordant=0;
	for(vector<ReadRec_t>::const_iterator it=SBamrecord.begin(); it!=SBamrecord.end(); it++){
		if(it->IsDiscordant()){
			numDiscordant++;
			output<<it->Qname<<endl;
		}
	}
	output.close();
	return numDiscordant;
};
