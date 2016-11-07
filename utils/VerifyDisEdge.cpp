#include "BPNode.h"
#include "BPEdge.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace BamTools;

typedef tuple<int,int,int,int> INV_t; // Chr, StartPos, EndPos, ID
typedef tuple<int,int,int,int,int,int> INS_t; // INSChr, INSPos, DELChr, DELStartPos, DELEndPos, ID
typedef tuple<int,int,bool,int,int,bool,int> TRA_t; //Chr1, Pos1, IsLeft1, Chr2, Pos2, IsLeft2, ID

void BuildRefName(string bamfile, vector<string>& RefName, std::map<string,int> & RefTable){
    BamReader bamreader; bamreader.Open(bamfile);
    if(bamreader.IsOpen()){
        int count=0;
        SamHeader header=bamreader.GetHeader();
        for(SamSequenceIterator it=header.Sequences.Begin();it!=header.Sequences.End();it++){
            RefName.push_back(it->Name);
            RefTable[it->Name]=count++;
        }
    }
    else
        cout<<"Cannot open bamfile "<<bamfile<<endl;
};

void ReadSVfile(string file, const map<string,int>& RefTable, vector< vector<int> >& ChrBreakpoint, vector<INV_t>& INV, vector<INS_t>& INS){
    ChrBreakpoint.resize(RefTable.size());
    ifstream input(file);
    string line;
    // pass first few records that are insertions without reference
    while(getline(input, line)){
        if(line.substr(0,3)!="INS" && line[0]!='#')
            break;
    }
    // insertions
    do{
        if(line.substr(0,3)=="INV")
            break;
        vector<string> strs;
        boost::split(strs, line, boost::is_any_of("\t"));
        int chr=RefTable.at(strs[1]);
        ChrBreakpoint[chr].push_back(stoi(strs[2]));
        if(strs[0]=="DEL"){
            ChrBreakpoint[chr].push_back(stoi(strs[5]));
            INS.push_back(make_tuple(0,0,chr, stoi(strs[2]), stoi(strs[5]), stoi(strs[7])));
        }
        else if(strs[0]=="INS"){
            get<0>(INS.back())=chr;
            get<1>(INS.back())=stoi(strs[2]);
        }
    }while(getline(input, line));
    // inversions
    do{
        if(line.substr(0,3)!="INV")
            break;
        vector<string> strs;
        boost::split(strs, line, boost::is_any_of("\t"));
        int chr=RefTable.at(strs[1]);
        ChrBreakpoint[chr].push_back(stoi(strs[2]));
        ChrBreakpoint[chr].push_back(stoi(strs[5]));
        INV.push_back(make_tuple(chr, stoi(strs[2]), stoi(strs[5]), stoi(strs[7])));
    }while(getline(input, line));
    // skip the next deletions
    do{
        if(line.substr(0,3)!="DEL")
            break;
    }while(getline(input, line));
    // the last is TRA
    do{
        vector<string> strs;
        boost::split(strs, line, boost::is_any_of("\t"));
        ChrBreakpoint[RefTable.at(strs[1])].push_back(stoi(strs[2]));
        ChrBreakpoint[RefTable.at(strs[4])].push_back(stoi(strs[5]));
    }while(getline(input, line));

    for(int i=0; i<ChrBreakpoint.size(); i++)
        if(ChrBreakpoint[i].size()!=0)
            sort(ChrBreakpoint[i].begin(), ChrBreakpoint[i].end());
};

vector< pair<int,int> > VerifyEdges(string file, vector< vector<int> >& ChrBreakpoint, vector<INV_t>& INV, vector<INS_t>& INS, vector<bool>& SVHit, vector<int>& TPDistances){
    int thresh=3; // threshold for 3 bp, if end of region is very near breakpoint
    SVHit.clear(); SVHit.resize((int)INV.size()+(int)INS.size(), false);
    TPDistances.clear(); TPDistances.reserve(4*(int)INV.size()+4*(int)INS.size());
    vector< pair<int,int> > VerifyResult;
    ifstream input(file);
    string line;
    while(getline(input, line)){
        if(line[0]=='#')
            continue;
        vector<string> strs1, strs2;
        boost::split(strs1, line, boost::is_any_of("\t"));
        boost::split(strs2, strs1[1], boost::is_any_of(","));
        int chr1=stoi(strs2[0]), pos1=(strs1[2]=="H")?stoi(strs2[1]):(stoi(strs2[1])+stoi(strs2[2]));
        boost::split(strs2, strs1[4], boost::is_any_of(","));
        int chr2=stoi(strs2[0]), pos2=(strs1[5]=="H")?stoi(strs2[1]):(stoi(strs2[1])+stoi(strs2[2]));
        bool isfind=false;
        int HitID=-1;
        int tppos1=-1, tppos2=-1;
        if(chr1==chr2){ // search INV first
            for(int i=0; i<INV.size(); i++){
                if(chr1==get<0>(INV[i]) && get<1>(INV[i])-thresh<=pos1 && get<2>(INV[i])+thresh>=pos1){
                    if(pos2>pos1 && pos2>get<2>(INV[i])-thresh){
                        vector<int>::iterator low=lower_bound(ChrBreakpoint[chr1].begin(), ChrBreakpoint[chr1].end(), get<2>(INV[i]));
                        low++;
                        if((low==ChrBreakpoint[chr1].end() || pos2<=(*low)) && strs1[2]=="H" && strs1[5]=="H"){
                            isfind=true; SVHit[i]=true; HitID=get<3>(INV[i]); tppos1=get<1>(INV[i]); tppos2=get<2>(INV[i]); break;
                        }
                    }
                    else if(pos2<pos1 && pos2<get<1>(INV[i])+thresh){
                        vector<int>::iterator low=lower_bound(ChrBreakpoint[chr1].begin(), ChrBreakpoint[chr1].end(), get<1>(INV[i]));
                        if((low==ChrBreakpoint[chr1].begin() || *(low-1)<pos2) && strs1[2]=="T" && strs1[5]=="T"){
                            isfind=true; SVHit[i]=true; HitID=get<3>(INV[i]); tppos1=get<2>(INV[i]); tppos2=get<1>(INV[i]); break;
                        }
                    }
                }
                if(chr2==get<0>(INV[i]) && get<1>(INV[i])-thresh<=pos2 && get<2>(INV[i])+thresh>=pos2){
                    if(pos1>pos2 && pos1>get<2>(INV[i])-thresh){
                        vector<int>::iterator low=lower_bound(ChrBreakpoint[chr2].begin(), ChrBreakpoint[chr2].end(), get<2>(INV[i]));
                        low++;
                        if((low==ChrBreakpoint[chr2].end() || pos1<=(*low)) && strs1[2]=="H" && strs1[5]=="H"){
                            isfind=true; SVHit[i]=true; HitID=get<3>(INV[i]); tppos1=get<2>(INV[i]); tppos2=get<1>(INV[i]); break;
                        }
                    }
                    else if(pos1<pos2 && pos1<get<1>(INV[i])+thresh){
                        vector<int>::iterator low=lower_bound(ChrBreakpoint[chr2].begin(), ChrBreakpoint[chr2].end(), get<1>(INV[i]));
                        if((low==ChrBreakpoint[chr2].begin() || *(low-1)<pos1) && strs1[2]=="T" && strs1[5]=="T"){
                            isfind=true; SVHit[i]=true; HitID=get<3>(INV[i]); tppos1=get<1>(INV[i]); tppos2=get<2>(INV[i]); break;
                        }
                    }
                }
            }
        }
        if(!isfind){
            for(int i=0; i<INS.size(); i++){
                if(chr1==get<2>(INS[i]) && get<3>(INS[i])-thresh<=pos1 && get<4>(INS[i])+thresh>=pos1){
                    vector<int>::iterator low=lower_bound(ChrBreakpoint[get<0>(INS[i])].begin(), ChrBreakpoint[get<0>(INS[i])].end(), get<1>(INS[i]));
                    if(chr2==get<0>(INS[i]) && (low==ChrBreakpoint[get<0>(INS[i])].begin() || pos2>=*(low-1)) && (low+1==ChrBreakpoint[get<0>(INS[i])].end() || pos2<=*(low+1))){
                        isfind=true; SVHit[(int)INV.size()+i]=true; HitID=get<5>(INS[i]); tppos1=(strs1[2]=="H")?get<3>(INS[i]):get<4>(INS[i]); tppos2=get<1>(INS[i]); break;
                    }
                }
                if(chr2==get<2>(INS[i]) && get<3>(INS[i])-thresh<=pos2 && get<4>(INS[i])+thresh>=pos2){
                    vector<int>::iterator low=lower_bound(ChrBreakpoint[get<0>(INS[i])].begin(), ChrBreakpoint[get<0>(INS[i])].end(), get<1>(INS[i]));
                    if(chr1==get<0>(INS[i]) && (low==ChrBreakpoint[get<0>(INS[i])].begin() || pos1>=*(low-1)) && (low+1==ChrBreakpoint[get<0>(INS[i])].end() || pos1<=*(low+1))){
                        isfind=true; SVHit[(int)INV.size()+i]=true; HitID=get<5>(INS[i]); tppos1=get<1>(INS[i]); tppos2=(strs1[5]=="H")?get<3>(INS[i]):get<4>(INS[i]); break;
                    }
                }
            }
        }
        if(isfind){
            VerifyResult.push_back(make_pair(stoi(strs1[6]), HitID));
            if(tppos1!=-1 && tppos2!=-1){
                TPDistances.push_back(abs(tppos1-pos1));
                TPDistances.push_back(abs(tppos2-pos2));
            }
            else
                cout<<"Error: true position being -1 but still report isfind=true"<<endl;
        }
        else
            VerifyResult.push_back(make_pair(stoi(strs1[6]), -1));
    }
    return VerifyResult;
};

void WriteResult(string infile, string outfile, vector< pair<int,int> >& VerifyResult){
    ifstream input(infile);
    ofstream output(outfile, ios::out);
    string line;
    int count=0;
    while(getline(input, line)){
        if(line[0]=='#')
            output<<line<<'\t'<<"HitorNot\tHittedSVID\n";
        else{
            output<<line<<'\t'<<(VerifyResult[count].second>=0)<<'\t'<<VerifyResult[count].second<<endl; count++;
        }
    }
    input.close(); output.close();
};

void WriteMetrics(string outfile, vector<INV_t>& INV, vector<INS_t>& INS, vector< pair<int,int> >& VerifyResult){
    int TotalSV=(int)INV.size()+(int)INS.size();
    ofstream output(outfile, ios::out);
    output<<"## Weight\tTruePos\tFalsepos\tPrecision\tNumHit\tSensitivity\n";
    vector<int> uniqweight;
    for(int i=0; i<VerifyResult.size(); i++)
        uniqweight.push_back(VerifyResult[i].first);
    sort(uniqweight.begin(), uniqweight.end());
    vector<int>::iterator endit=unique(uniqweight.begin(), uniqweight.end());
    uniqweight.resize(distance(uniqweight.begin(), endit));
    for(int i=0; i<uniqweight.size(); i++){
        int truepos=0, falsepos=0, numhit=0;
        vector<int> HitID;
        for(int j=0; j<VerifyResult.size(); j++){
            if(VerifyResult[j].first>=uniqweight[i] && VerifyResult[j].second==-1)
                falsepos++;
            else if(VerifyResult[j].first>=uniqweight[i]){
                truepos++;
                HitID.push_back(VerifyResult[j].second);
            }
        }
        sort(HitID.begin(), HitID.end());
        endit=unique(HitID.begin(), HitID.end());
        HitID.resize(distance(HitID.begin(), endit));
        output<<uniqweight[i]<<'\t'<<truepos<<'\t'<<falsepos<<'\t'<<(1.0*truepos/(truepos+falsepos))<<'\t'<<(int)HitID.size()<<'\t'<<(1.0*(double)HitID.size()/TotalSV)<<endl;
    }
    output.close();
};

void WriteDistances(string outfile, vector<int>& TPDistances){
    FILE *fp;
    fp=fopen(outfile.c_str(), "wb");
    for(int i=0; i<TPDistances.size(); i++)
        fwrite(&TPDistances[i], sizeof(int), 1, fp);
    fclose(fp);
};

int main(int argc, char* argv[]){
    map<string, int> RefTable;
    vector<string> RefName;
    BuildRefName(argv[1], RefName, RefTable);

    vector<INV_t> INV;
    vector<INS_t> INS;
    vector< vector<int> > ChrBreakpoint;
    ReadSVfile(argv[2], RefTable, ChrBreakpoint, INV, INS);

    vector<bool> SVHit;
    vector<int> TPDistances;
    vector< pair<int,int> > VerifyResult=VerifyEdges(argv[3], ChrBreakpoint, INV, INS, SVHit, TPDistances);
    string hitoutput(argv[3]);
    string metricoutput(argv[3]);
    string distoutput(argv[3]);
    string suffix;
    if(argc>4){
        suffix=argv[4];
        suffix="_"+suffix;
    }
    else
        suffix="";
    size_t start=hitoutput.find_last_of(".");
    hitoutput=hitoutput.substr(0,start)+"_hit.txt";
    start=metricoutput.find_last_of("/");
    metricoutput=metricoutput.substr(0, start)+"/MetricTable"+suffix+".txt";
    distoutput=distoutput.substr(0,start)+"/TPDist"+suffix+".dat";
    WriteResult(argv[3], hitoutput, VerifyResult);
    WriteMetrics(metricoutput, INV, INS, VerifyResult);
    WriteDistances(distoutput, TPDistances);
    int hitcount=0;
    for(int i=0; i<SVHit.size(); i++)
        if(SVHit[i])
            hitcount++;
    cout<<"Sensitivity = "<<hitcount<<"/"<<SVHit.size()<<" = "<<(1.0*hitcount/(int)SVHit.size())<<endl;
}
