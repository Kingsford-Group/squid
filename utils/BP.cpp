#include "BP.h"

using namespace std;

string Breakpoint_t::Print(){
    string info=to_string(Chr)+'\t'+to_string(Position)+'\t'+((IsLeft)?"left":"right");
    return info;
};

SVset_t::SVset_t(string filetype, string svfile, const map<string,int>& RefTable,int _FilterLowQual){
    transform(filetype.begin(), filetype.end(), filetype.begin(), ::toupper);
    if(filetype=="VCF")
        VCFRead(svfile, RefTable, _FilterLowQual);
    else if(filetype=="BD")
        BDRead(svfile, RefTable, _FilterLowQual);
    else if(filetype=="BEDPE")
        BedpeRead(svfile, RefTable, _FilterLowQual);
};

void SVset_t::VCFRead(string svfile, const map<string,int>& RefTable, int _FilterLowQual){
    ifstream input(svfile);
    string line;
    Breakpoint_t bp1,bp2;
    while(getline(input,line)){
        if(line[0]=='#')
            continue;
        vector<string> strs;
        boost::split(strs,line,boost::is_any_of("\t"));
        size_t flagpos=strs[4].find_first_of("<>[]");
        if(flagpos== string::npos)
            continue;
        if(_FilterLowQual>0 && strs[6]=="LowQual")
            continue;
        if(strs[4]=="<DEL>"){
            bp1.Chr=RefTable.at(strs[0]); bp1.Position=stoi(strs[1]); bp1.IsLeft=true;
            size_t start=strs[7].find(";END=");
            size_t stop=strs[7].find_first_of(";",start+1);
            bp2.Chr=RefTable.at(strs[0]); bp2.Position=stoi(strs[7].substr(start+5, stop-start-5)); bp2.IsLeft=false;
            vBP.push_back(make_tuple(bp1, bp2, SVCode::DEL, (strs[6]=="LowQual")?0:1, strs[2]));
        }
        else if(strs[4]=="<DUP>"){
            bp1.Chr=RefTable.at(strs[0]); bp1.Position=stoi(strs[1]); bp1.IsLeft=false;
            size_t start=strs[7].find(";END=");
            size_t stop=strs[7].find_first_of(";",start+1);
            bp2.Chr=RefTable.at(strs[0]); bp2.Position=stoi(strs[7].substr(start+5, stop-start-5)); bp2.IsLeft=true;
            vBP.push_back(make_tuple(bp1,bp2,SVCode::DUP, (strs[6]=="LowQual")?0:1, strs[2]));
        }
        else if(strs[4]=="<INV>"){
            bp1.Chr=RefTable.at(strs[0]); bp1.Position=stoi(strs[1]); bp1.IsLeft=false;
            size_t start=strs[7].find(";END=");
            size_t stop=strs[7].find_first_of(";",start+1);
            bp2.Chr=RefTable.at(strs[0]); bp2.Position=stoi(strs[7].substr(start+5, stop-start-5)); bp2.IsLeft=false;
            vBP.push_back(make_tuple(bp1, bp2, SVCode::INV, (strs[6]=="LowQual")?0:1, strs[2]));
        }
        else if(strs[4]=="<TRA>"){
            bp1.Chr=RefTable.at(strs[0]); bp1.Position=stoi(strs[1]);
            size_t start=strs[7].find(";CHR2=");
            size_t stop=strs[7].find_first_of(";",start+1);
            string tmp=strs[7].substr(start+6, stop-start-6);
            bp2.Chr=RefTable.at(tmp);
            start=strs[7].find(";END="); stop=strs[7].find_first_of(";",start+1);
            bp2.Position=stoi(strs[7].substr(start+5, stop-start-5));
            start=strs[7].find(";CT="); stop=strs[7].find_first_of(";",start+1);
            tmp=strs[7].substr(start+4, stop-start-4);
            if(tmp=="3to3"){
                bp1.IsLeft=true; bp2.IsLeft=true;
            }
            else if(tmp=="3to5"){
                bp1.IsLeft=true; bp2.IsLeft=false;
            }
            else if(tmp=="5to3"){
                bp1.IsLeft=false; bp2.IsLeft=true;
            }
            else{
                bp1.IsLeft=false; bp2.IsLeft=false;
            }
            vBP.push_back(make_tuple(bp1, bp2, SVCode::TRA, (strs[6]=="LowQual")?0:1, strs[2]));
        }
        else{
            bp1.Chr=RefTable.at(strs[0]); bp1.Position=stoi(strs[1]);
            size_t start=strs[4].find_first_of("][");
            size_t stop=strs[4].find_first_of("][", start+1);
            if(start==0)
                bp1.IsLeft=false;
            else
                bp1.IsLeft=true;
            vector<string> strs2;
            string ALT=strs[4].substr(start, stop-start);
            boost::split(strs2, ALT, boost::is_any_of(":"));
            bp2.Chr=RefTable.at(strs2[0].substr(1)); bp2.Position=stoi(strs2[1]);
            if(strs[4][start]==']')
                bp2.IsLeft=true;
            else
                bp2.IsLeft=false;
            vBP.push_back(make_tuple(bp1, bp2, SVCode::TRA, (strs[6]=="LowQual")?0:1, strs[2]));
        }
    }
};

void SVset_t::BDRead(string svfile, const map<string,int>& RefTable, int _FilterLowQual){
    ifstream input(svfile);
    string line;
    Breakpoint_t bp1,bp2;
    int count=0;
    while(getline(input,line)){
        if(line[0]=='#')
            continue;
        count++;
        vector<string> strs;
        boost::split(strs,line,boost::is_any_of("\t"));
        if(stoi(strs[9])<_FilterLowQual)
            continue;
        if(strs[6]=="DEL"){
            bp1.Chr=RefTable.at(strs[0]); bp1.Position=stoi(strs[1]); bp1.IsLeft=true;
            bp2.Chr=RefTable.at(strs[3]); bp2.Position=stoi(strs[4]); bp2.IsLeft=false;
            vBP.push_back(make_tuple(bp1, bp2, SVCode::DEL, stoi(strs[9]), to_string(count)));
        }
        else if(strs[6]=="INV"){
            bp1.Chr=RefTable.at(strs[0]); bp1.Position=stoi(strs[1]); bp1.IsLeft=false;
            bp2.Chr=RefTable.at(strs[3]); bp2.Position=stoi(strs[4]); bp2.IsLeft=false;
            vBP.push_back(make_tuple(bp1, bp2, SVCode::INV, stoi(strs[9]), to_string(count)));
        }
        else if(strs[6]=="ITX"){
            bp1.Chr=RefTable.at(strs[0]); bp1.Position=stoi(strs[1]); bp1.IsLeft=false;
            bp2.Chr=RefTable.at(strs[3]); bp2.Position=stoi(strs[4]); bp2.IsLeft=false;
            vBP.push_back(make_tuple(bp1, bp2, SVCode::TRA, stoi(strs[9]), to_string(count)));
        }
        else if(strs[6]=="CTX"){
            bp1.Chr=RefTable.at(strs[0]); bp1.Position=stoi(strs[1]); bp1.IsLeft=false;
            bp2.Chr=RefTable.at(strs[3]); bp2.Position=stoi(strs[4]); bp2.IsLeft=false;
            vBP.push_back(make_tuple(bp1, bp2, SVCode::TRA, stoi(strs[9]), to_string(count)));
        }
    }
};

void SVset_t::BedpeRead(string svfile, const map<string,int>& RefTable, int _FilterLowQual){
    ifstream input(svfile);
    string line;
    Breakpoint_t bp1,bp2;
    while(getline(input,line)){
        if(line[0]=='#')
            continue;
        vector<string> strs;
        boost::split(strs,line,boost::is_any_of("\t"));
        if(stoi(strs[7])<_FilterLowQual)
            continue;
        if(strs[10]=="local_deletion"){
            bp1.Chr=RefTable.at(strs[0]); bp1.Position=(strs[8]=="+")?stoi(strs[2]):stoi(strs[1]); bp1.IsLeft=(strs[8]=="+")?true:false;
            bp2.Chr=RefTable.at(strs[3]); bp2.Position=(strs[9]=="+")?stoi(strs[5]):stoi(strs[4]); bp2.IsLeft=(strs[9]=="+")?true:false;
            vBP.push_back(make_tuple(bp1, bp2, SVCode::DEL, stoi(strs[7]), strs[6]));
        }
        else if(strs[10]=="distant_deletion"){
            bp1.Chr=RefTable.at(strs[0]); bp1.Position=(strs[8]=="+")?stoi(strs[2]):stoi(strs[1]); bp1.IsLeft=(strs[8]=="+")?true:false;
            bp2.Chr=RefTable.at(strs[3]); bp2.Position=(strs[9]=="+")?stoi(strs[5]):stoi(strs[4]); bp2.IsLeft=(strs[9]=="+")?true:false;
            vBP.push_back(make_tuple(bp1, bp2, SVCode::DEL, stoi(strs[7]), strs[6]));
        }
        else if(strs[10]=="local_inversion" || strs[10]=="distant_inversion"){
            bp1.Chr=RefTable.at(strs[0]); bp1.Position=(strs[8]=="+")?stoi(strs[2]):stoi(strs[1]); bp1.IsLeft=(strs[8]=="+")?true:false;
            bp2.Chr=RefTable.at(strs[3]); bp2.Position=(strs[9]=="+")?stoi(strs[5]):stoi(strs[4]); bp2.IsLeft=(strs[9]=="+")?true:false;
            vBP.push_back(make_tuple(bp1, bp2, SVCode::INV, stoi(strs[7]), strs[6]));
        }
        else if(strs[10]=="local_duplication" || strs[10]=="distant_duplication"){
            bp1.Chr=RefTable.at(strs[0]); bp1.Position=(strs[8]=="+")?stoi(strs[2]):stoi(strs[1]); bp1.IsLeft=(strs[8]=="+")?true:false;
            bp2.Chr=RefTable.at(strs[3]); bp2.Position=(strs[9]=="+")?stoi(strs[5]):stoi(strs[4]); bp2.IsLeft=(strs[9]=="+")?true:false;
            vBP.push_back(make_tuple(bp1, bp2, SVCode::DUP, stoi(strs[7]), strs[6]));
        }
        else{
            bp1.Chr=RefTable.at(strs[0]); bp1.Position=(strs[8]=="+")?stoi(strs[2]):stoi(strs[1]); bp1.IsLeft=(strs[8]=="+")?true:false;
            bp2.Chr=RefTable.at(strs[3]); bp2.Position=(strs[9]=="+")?stoi(strs[5]):stoi(strs[4]); bp2.IsLeft=(strs[9]=="+")?true:false;
            vBP.push_back(make_tuple(bp1, bp2, SVCode::OTH, stoi(strs[7]), strs[6]));
        }
    }
};

int SVset_t::CountafterFilter(int _FilterDeletion){
    if(_FilterDeletion>1){
        int count=0;
        for(vector< tuple<Breakpoint_t, Breakpoint_t, SVCode, int, string> >::iterator it=vBP.begin(); it!=vBP.end(); it++)
            if(get<2>(*it)!=SVCode::DEL)
                count++;
        return count;
    }
    else if(_FilterDeletion==1){
        vector<Breakpoint_t> tmpBP, nondelBP;
        vector< tuple<Breakpoint_t, Breakpoint_t, SVCode, int, string> > newvBP; newvBP.reserve(vBP.size());
        for(vector< tuple<Breakpoint_t, Breakpoint_t, SVCode, int, string> >::iterator it=vBP.begin(); it!=vBP.end(); it++){
            if(get<2>(*it)==SVCode::DEL){
                tmpBP.push_back(get<0>(*it));
                tmpBP.push_back(get<1>(*it));
            }
            else{
                nondelBP.push_back(get<0>(*it));
                nondelBP.push_back(get<1>(*it));
            }
        }
        sort(tmpBP.begin(), tmpBP.end());
        sort(nondelBP.begin(), nondelBP.end());
        for(vector< tuple<Breakpoint_t, Breakpoint_t, SVCode, int, string> >::iterator it=vBP.begin(); it!=vBP.end(); it++){
            if(get<2>(*it)!=SVCode::DEL)
                newvBP.push_back(*it);
            else{
                vector<Breakpoint_t>::iterator it1=upper_bound(nondelBP.begin(),nondelBP.end(), get<0>(*it));
                if(it1!=nondelBP.end() && (*it1)<get<1>(*it) && get<0>(*it)<(*it1)) // if DEL overlap with non-DEL SV, keep that DEL; otherwise, don't keep that
                    newvBP.push_back(*it);
            }
        }
        newvBP.reserve(newvBP.size());
        return (int)newvBP.size();
    }
    else
        return (int)vBP.size();
};

vector<int> SVset_t::VerifyBPs(vector<int>& RefLength, vector<INV_t>& INV, vector<INS_t>& INS){
    int thresh=3;
    vector<int> VerifyResult;
    vector< vector<int> > thisChrBP, ChrBP; // thisChrBP records all breakpoints of each chromosome, for building nodes; ChrBP records all true breakpoints, for compairng position
    thisChrBP.resize(RefLength.size()); ChrBP.resize(RefLength.size());
    for(vector< tuple<Breakpoint_t, Breakpoint_t, SVCode, int, string> >::iterator it=vBP.begin(); it!=vBP.end(); it++){
        thisChrBP[get<0>(*it).Chr].push_back(get<0>(*it).Position);
        thisChrBP[get<1>(*it).Chr].push_back(get<1>(*it).Position);
    }
    for(vector<INV_t>::iterator it=INV.begin(); it!=INV.end(); it++){
        ChrBP[get<0>(*it)].push_back(get<1>(*it));
        ChrBP[get<0>(*it)].push_back(get<2>(*it));
    }
    for(vector<INS_t>::iterator it=INS.begin(); it!=INS.end(); it++){
        ChrBP[get<0>(*it)].push_back(get<1>(*it));
        ChrBP[get<2>(*it)].push_back(get<3>(*it)); ChrBP[get<2>(*it)].push_back(get<4>(*it));
    }
    for(int i=0; i<RefLength.size(); i++){
        if(thisChrBP[i].size()!=0)
            sort(thisChrBP[i].begin(), thisChrBP[i].end());
        if(ChrBP[i].size()!=0)
            sort(ChrBP[i].begin(), ChrBP[i].end());
    }
    // for each predicted BP in vBP, compare with INV and INS
    for(vector< tuple<Breakpoint_t, Breakpoint_t, SVCode, int, string> >::iterator it=vBP.begin(); it!=vBP.end(); it++){
        int chr1=get<0>(*it).Chr, chr2=get<1>(*it).Chr, pos1=get<0>(*it).Position, pos2=get<1>(*it).Position;
        // verify linkage
        bool isfind=false;
        int HitID=-1;
        if(chr1==chr2){ // search INV first, require one bp inside true INV, one bp outside
            for(int i=0; i<INV.size(); i++){
                if(chr1==get<0>(INV[i]) && get<1>(INV[i])-thresh<=pos1 && get<2>(INV[i])+thresh>=pos1){
                    if(pos2>pos1 && pos2>get<2>(INV[i])-thresh){
                        vector<int>::iterator low=lower_bound(ChrBP[chr1].begin(), ChrBP[chr1].end(), get<2>(INV[i]));
                        low++;
                        if((low==ChrBP[chr1].end() || pos2<=(*low)) && get<0>(*it).IsLeft==get<1>(*it).IsLeft){
                            isfind=true; HitID=get<3>(INV[i]); break;
                        }
                    }
                    else if(pos2<pos1 && pos2<get<1>(INV[i])+thresh){
                        vector<int>::iterator low=lower_bound(ChrBP[chr1].begin(), ChrBP[chr1].end(), get<1>(INV[i]));
                        if((low==ChrBP[chr1].begin() || *(low-1)<pos2) && get<0>(*it).IsLeft==get<1>(*it).IsLeft){
                            isfind=true; HitID=get<3>(INV[i]); break;
                        }
                    }
                }
                if(chr2==get<0>(INV[i]) && get<1>(INV[i])-thresh<=pos2 && get<2>(INV[i])+thresh>=pos2){
                    if(pos1>pos2 && pos1>get<2>(INV[i])-thresh){
                        vector<int>::iterator low=lower_bound(ChrBP[chr2].begin(), ChrBP[chr2].end(), get<2>(INV[i]));
                        low++;
                        if((low==ChrBP[chr2].end() || pos1<=(*low)) && get<0>(*it).IsLeft==get<1>(*it).IsLeft){
                            isfind=true; HitID=get<3>(INV[i]); break;
                        }
                    }
                    else if(pos1<pos2 && pos1<get<1>(INV[i])+thresh){
                        vector<int>::iterator low=lower_bound(ChrBP[chr2].begin(), ChrBP[chr2].end(), get<1>(INV[i]));
                        if((low==ChrBP[chr2].begin() || *(low-1)<pos1) && get<0>(*it).IsLeft==get<1>(*it).IsLeft){
                            isfind=true; HitID=get<3>(INV[i]); break;
                        }
                    }
                }
            }
        }
        if(!isfind){ // search INS, require one bp within deletion, and one bp within lb and ub of insertion point
            for(int i=0; i<INS.size(); i++){
                if(chr1==get<2>(INS[i]) && get<3>(INS[i])-thresh<=pos1 && get<4>(INS[i])+thresh>=pos1){
                    vector<int>::iterator low=lower_bound(ChrBP[get<0>(INS[i])].begin(), ChrBP[get<0>(INS[i])].end(), get<1>(INS[i]));
                    if(chr2==get<0>(INS[i]) && (low==ChrBP[get<0>(INS[i])].begin() || pos2>=*(low-1)) && (low+1==ChrBP[get<0>(INS[i])].end() || pos2<=*(low+1))){
                        isfind=true; HitID=get<5>(INS[i]); break;
                    }
                }
                if(chr2==get<2>(INS[i]) && get<3>(INS[i])-thresh<=pos2 && get<4>(INS[i])+thresh>=pos2){
                    vector<int>::iterator low=lower_bound(ChrBP[get<0>(INS[i])].begin(), ChrBP[get<0>(INS[i])].end(), get<1>(INS[i]));
                    if(chr1==get<0>(INS[i]) && (low==ChrBP[get<0>(INS[i])].begin() || pos1>=*(low-1)) && (low+1==ChrBP[get<0>(INS[i])].end() || pos1<=*(low+1))){
                        isfind=true; HitID=get<5>(INS[i]); break;
                    }
                }
                /*if(chr1==chr2 && ((pos1<get<3>(INS[i])+thresh && pos2>get<4>(INS[i])-thresh) || (pos2<get<3>(INS[i])+thresh && pos1>get<4>(INS[i])-thresh))){ // deletion part of INS
                    vector<int>::iterator low1=lower_bound(ChrBP[get<2>(INS[i])].begin(), ChrBP[get<2>(INS[i])].end(), get<3>(INS[i]));
                    vector<int>::iterator low2=lower_bound(ChrBP[get<2>(INS[i])].begin(), ChrBP[get<2>(INS[i])].end(), get<4>(INS[i]));
                    if((low1==ChrBP[get<2>(INS[i])].begin() || min(pos1,pos2)>=*(low1-1)) && (low2+1==ChrBP[get<2>(INS[i])].end() || max(pos1,pos2)<=*(low2+1))){
                        isfind=true; HitID=-2; hasInsDelonly=true; break;
                    }
                }*/
            }
        }
        if(isfind)
            VerifyResult.push_back(HitID);
        else
            VerifyResult.push_back(-1);
    }
    return VerifyResult;
};

double SVset_t::CalPrecision(vector<int>& SVHitID, int _FilterDeletion){
    int TrueHit=0;
    vector<int> UniqID;
    for(int i=0; i<SVHitID.size(); i++)
        if(SVHitID[i]!=-1)
            TrueHit++;
    int denominator=CountafterFilter(_FilterDeletion);
    double Precision=1.0*TrueHit/denominator;
    if(Precision>1)
        Precision=1;
    return Precision;
};

double SVset_t::CalSensitivity(vector<int>& SVHitID, vector<INV_t>& INV, vector<INS_t>& INS){
    int TrueHit=0;
    vector<int> UniqID;
    for(int i=0; i<SVHitID.size(); i++)
        if(SVHitID[i]!=-1){
            TrueHit++; UniqID.push_back(SVHitID[i]);
        }
    sort(UniqID.begin(), UniqID.end());
    vector<int>::iterator endit=unique(UniqID.begin(), UniqID.end());
    UniqID.resize(distance(UniqID.begin(), endit));
    double Sensitivity=1.0*(double)UniqID.size()/((double)INV.size()+(double)INS.size());
    return Sensitivity;
};

void SVset_t::Output(string outputfile){
    ofstream output(outputfile, ios::out);
    output<<"## Type\tID\tChr1\tPosition1\tIsLeft1\tChr2\tPosition2\tIsLeft2\n";
    for(vector< tuple<Breakpoint_t, Breakpoint_t, SVCode, int, string> >::iterator it=vBP.begin(); it!=vBP.end(); it++){
        output<<SVString[get<2>(*it)]<<'\t'<<get<3>(*it)<<'\t'<<get<0>(*it).Print()<<'\t'<<get<1>(*it).Print()<<endl;
    }
    output.close();
};

void SVset_t::Output(string outputfile, vector<int>& VerifyResult, double & Precision, double & Sensitivity){
    ofstream output(outputfile, ios::out);
    output<<"## Precision = "<<Precision<<"\tSensitivity = "<<Sensitivity<<endl;
    output<<"## Type\tID\tChr1\tPosition1\tIsLeft1\tChr2\tPosition2\tIsLeft2\tQuality\tHitorNot\tHitID\n";
    int i=0;
    for(vector< tuple<Breakpoint_t, Breakpoint_t, SVCode, int, string> >::iterator it=vBP.begin(); it!=vBP.end(); it++){
        output<<SVString[get<2>(*it)]<<'\t'<<get<4>(*it)<<'\t'<<get<0>(*it).Print()<<'\t'<<get<1>(*it).Print()<<'\t'<<get<3>(*it)<<'\t'<<(VerifyResult[i]>=0)<<'\t'<<VerifyResult[i]<<endl;
        i++;
    }
    output.close();
};
