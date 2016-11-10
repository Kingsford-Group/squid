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
    if(PhredType){
        for(int i=0; i<record.Qualities.size(); i++){
            if(record.Qualities[i]<'$')
                tmpLowPhredLen++;
            else
                tmpLowPhredLen=0;
            if(LowPhredLen<tmpLowPhredLen)
                LowPhredLen=tmpLowPhredLen;
        }
    }
    else{
        for(int i=0; i<record.Qualities.size(); i++){
            if(record.Qualities[i]<'C')
                tmpLowPhredLen++;
            else
                tmpLowPhredLen=0;
            if(LowPhredLen<tmpLowPhredLen)
                LowPhredLen=tmpLowPhredLen;
        }
    }
    if(record.IsFirstMate()){
        FirstTotalLen=TotalLen; SecondTotalLen=0; FirstLowPhred=(LowPhredLen>LowPhredLenThresh);
    }
    else{
        SecondTotalLen=TotalLen; FirstTotalLen=0; SecondLowPhred=(LowPhredLen>LowPhredLenThresh);
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
            SingleBamRec_t tmp(record.RefID, RefPos, ReadPos, tmpRef, tmpRead, record.MapQuality, record.IsReverseStrand(), record.IsFirstMate());
            if(record.IsReverseStrand())
                tmp.ReadPos=TotalLen-ReadPos-tmpRead;
            if(record.IsFirstMate())
                FirstRead.push_back(tmp);
            else
                SecondMate.push_back(tmp);
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
            for(int i=0; i<FirstRead.size()-1; i++){
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
            for(int i=0; i<SecondMate.size()-1; i++){
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

void ReadRec_t::ModifybyGraph(SegmentGraph_t& SegmentGraph, const vector< vector<int> >& Components, vector<int>& singleRead_Node, vector< pair<int, int> >& Node_NewChr){
    vector<int> tmpRead_Node=SegmentGraph.LocateRead(singleRead_Node, *this);
    singleRead_Node=tmpRead_Node;
    bool whetherdelete=false;
    for(int i=0; i<tmpRead_Node.size(); i++)
        if(tmpRead_Node[i]==-1)
            whetherdelete=true;
    if(whetherdelete){
        FirstRead.clear(); SecondMate.clear();
    }
    else{
        for(int i=0; i<FirstRead.size(); i++)
            FirstRead[i].ModifybyGraph_Single(SegmentGraph, Components, tmpRead_Node[i], Node_NewChr);
        for(int i=0; i<SecondMate.size(); i++)
            SecondMate[i].ModifybyGraph_Single(SegmentGraph, Components, tmpRead_Node[(int)FirstRead.size()+i], Node_NewChr);
    }
};

void SingleBamRec_t::ModifybyGraph_Single(SegmentGraph_t& SegmentGraph, const vector< vector<int> >& Components, int nodeidx, vector< pair<int, int> >& Node_NewChr){
    int oldRefPos=RefPos;
    int i=Node_NewChr[nodeidx].first, j=Node_NewChr[nodeidx].second;
    RefID=i;
    int offset=0;
    for(int k=0; k<j; k++)
        offset+=SegmentGraph.vNodes[abs(Components[i][k])-1].Length;
    if(Components[i][j]>0){
        RefPos=offset+oldRefPos-SegmentGraph.vNodes[nodeidx].Position;
    }
    else{
        IsReverse=!IsReverse;
        RefPos=offset+SegmentGraph.vNodes[nodeidx].Position+SegmentGraph.vNodes[nodeidx].Length-oldRefPos-MatchRef;
    }
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

void BuildReference(string fafile, const std::map<string,int>& RefTable, vector<int>& RefLength, vector<string>& RefSequence){
    ifstream input(fafile);
    RefLength.resize(RefTable.size());
    RefSequence.resize(RefTable.size());
    string line;
    int chrlength=0, count=0, curChr=0;
    string tmpseq;
    while(getline(input,line)){
        if(line[0]=='>'){
            vector<string> strs;
            boost::split(strs, line, boost::is_any_of("\t "));
            if(chrlength!=0){
                RefLength[curChr]=chrlength;
                RefSequence[curChr]=tmpseq;
            }
            curChr=RefTable.at(strs[0].substr(1));
            chrlength=0;
            tmpseq="";
        }
        else{
            chrlength+=line.size();
            tmpseq+=line;
        }
    }
    RefLength[curChr]=chrlength;
    RefSequence[curChr]=tmpseq;
    input.close();
};

void UpdateReference(const SegmentGraph_t& SegmentGraph, const vector< vector<int> >& Components, vector<int>& RefLength, vector<string>& RefSequence){
    vector<string> newRefSequence(RefSequence.size());
    for(int i=0; i<Components.size(); i++){
        int length=0;
        for(int j=0; j<Components[i].size(); j++){
            length+=SegmentGraph.vNodes[abs(Components[i][j])-1].Length;
            newRefSequence[i]+=RefSequence[SegmentGraph.vNodes[abs(Components[i][j])-1].Chr].substr(SegmentGraph.vNodes[abs(Components[i][j])-1].Position, SegmentGraph.vNodes[abs(Components[i][j])-1].Length);
        }
        RefLength[i]=length;
    }
    RefSequence=newRefSequence;
};

SBamrecord_t BuildBWASBamRecord(const std::map<string,int> & RefTable, string bamfile){
    time_t CurrentTime;
    string CurrentTimeStr;
    SBamrecord_t SBamrecord; SBamrecord.reserve(66536);
    SBamrecord_t newSBamrecord;
    vector<string> MultiAlignedName; MultiAlignedName.reserve(66536);
    BamReader bamreader; bamreader.Open(bamfile);
    if(bamreader.IsOpen()){
        time(&CurrentTime);
        CurrentTimeStr=ctime(&CurrentTime);
        cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Bam file is opened, start reading bam records.\n";
        BamAlignment record;
        while(bamreader.GetNextAlignment(record)){
            bool XAtag=record.HasTag("XA");
            if(record.IsMapped() && record.IsMateMapped() && !record.IsDuplicate() && record.MapQuality>0 && !XAtag){
                ReadRec_t tmp(record);
                SBamrecord.push_back(tmp);
            }
            else if(XAtag)
                MultiAlignedName.push_back(record.Name);
            if(SBamrecord.capacity()==SBamrecord.size())
                SBamrecord.reserve(SBamrecord.size()*2);
            if(MultiAlignedName.capacity()==MultiAlignedName.size())
                MultiAlignedName.reserve(MultiAlignedName.size()*2);
        }
        SBamrecord.reserve(SBamrecord.size());
        sort(SBamrecord.begin(), SBamrecord.end());
        sort(MultiAlignedName.begin(), MultiAlignedName.end());
        newSBamrecord.reserve(SBamrecord.size());
        for(SBamrecord_t::iterator it=SBamrecord.begin();it!=SBamrecord.end();it++){
            if(newSBamrecord.size()==0 || it->Qname!=newSBamrecord.back().Qname){
                newSBamrecord.push_back((*it));
                if(binary_search(MultiAlignedName.begin(), MultiAlignedName.end(), it->Qname))
                    newSBamrecord.back().MultiFilter=true;
            }
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
        SBamrecord.clear();
        for(SBamrecord_t::iterator it=newSBamrecord.begin(); it!=newSBamrecord.end(); it++){
            int xlen=0; // check whether first read is multi-aligned, but without XA tag
            for(vector<SingleBamRec_t>::iterator itsingle=it->FirstRead.begin(); itsingle!=it->FirstRead.end(); itsingle++)
                xlen+=itsingle->MatchRead;
            if(xlen>=2*it->FirstTotalLen){
                it->FirstRead.clear(); it->MultiFilter=true;
            }
            xlen=0; // check whether second mate is multi-aligned, but without XA tag
            for(vector<SingleBamRec_t>::iterator itsingle=it->SecondMate.begin(); itsingle!=it->SecondMate.end(); itsingle++)
                xlen+=itsingle->MatchRead;
            if(xlen>=2*it->SecondTotalLen){
                it->SecondMate.clear(); it->MultiFilter=true;
            }
            it->SortbyReadPos();
            SBamrecord.push_back((*it));
        }
        time(&CurrentTime);
        CurrentTimeStr=ctime(&CurrentTime);
        cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Finish reading bam records, start sorting\n";

        SBamrecord.reserve(SBamrecord.size());
        sort(SBamrecord.begin(), SBamrecord.end(), ReadRec_t::FrontSmallerThan);
        bamreader.Close();

        time(&CurrentTime);
        CurrentTimeStr=ctime(&CurrentTime);
        cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Finish sorting bam records\n";
    }
    else
        cout<<"Cannot open bamfile!"<<endl;
    return SBamrecord;
};

SBamrecord_t BuildMainSBamRecord(const std::map<string,int> & RefTable, string bamfile){
    SBamrecord_t SBamrecord; SBamrecord.reserve(66536);
    SBamrecord_t newSBamrecord;
    BamReader bamreader; bamreader.Open(bamfile);
    if(bamreader.IsOpen()){
        BamAlignment record;
        while(bamreader.GetNextAlignment(record)){
            if(record.IsMapped() && !record.IsDuplicate() && record.MapQuality==255){
                ReadRec_t tmp(record);
                SBamrecord.push_back(tmp);
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
        for(SBamrecord_t::iterator it=newSBamrecord.begin();it!=newSBamrecord.end();it++)
            it->SortbyReadPos();
        cout<<"Finish building bam record\tsize="<<newSBamrecord.size()<<endl;
    }
    else
        cout<<"Cannot open bamfile!"<<endl;
    bamreader.Close();
    return newSBamrecord;
};

SBamrecord_t BuildChimericSBamRecord(SBamrecord_t& SBamrecord, const std::map<string,int> & RefTable, string bamfile){
    SBamrecord_t tmpBamrecord; tmpBamrecord.reserve(66536);
    SBamrecord_t newSBamrecord;
    BamReader bamreader; bamreader.Open(bamfile);
    if(bamreader.IsOpen()){
        BamAlignment record;
        while(bamreader.GetNextAlignment(record)){
            if(record.IsMapped() && !record.IsDuplicate()){
                ReadRec_t tmp(record);
                tmpBamrecord.push_back(tmp);
                if(tmpBamrecord.capacity()==tmpBamrecord.size())
                    tmpBamrecord.reserve(tmpBamrecord.size()*2);
            }
        }
        tmpBamrecord.reserve(tmpBamrecord.size());
        sort(tmpBamrecord.begin(), tmpBamrecord.end());
        newSBamrecord.reserve(tmpBamrecord.size());
        for(SBamrecord_t::iterator it=tmpBamrecord.begin();it!=tmpBamrecord.end();it++){
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
        tmpBamrecord.clear(); tmpBamrecord.reserve((int)SBamrecord.size()+(int)newSBamrecord.size());
        SBamrecord_t::iterator itold=SBamrecord.begin();
        for(SBamrecord_t::iterator it=newSBamrecord.begin();it!=newSBamrecord.end();it++){
            it->SortbyReadPos();
            for(; itold!=SBamrecord.end() && itold->Qname<it->Qname; itold++)
                tmpBamrecord.push_back((*itold));
            if(itold!=SBamrecord.end() && itold->Qname==it->Qname){
                if(itold->ReadCoverageGap()>it->ReadCoverageGap())
                    tmpBamrecord.push_back(*it);
                else
                    tmpBamrecord.push_back(*itold);
                itold++;
            }
            else
                tmpBamrecord.push_back(*it);
        }
        for(; itold!=SBamrecord.end(); itold++)
            tmpBamrecord.push_back(*itold);
    }
    clock_t starttime=clock();
    sort(tmpBamrecord.begin(), tmpBamrecord.end(), ReadRec_t::FrontSmallerThan);
    cout<<"sorting time="<<(1.0*(clock()-starttime)/CLOCKS_PER_SEC)<<endl;
    return tmpBamrecord;
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
