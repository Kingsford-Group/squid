#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include <boost/algorithm/string.hpp>
#include "BP_VCF.h"

using namespace std;
using namespace BamTools;

void BuildRefName(string bamfile, vector<string>& RefName, map<string,int> & RefTable, vector<int>& RefLength){
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

void ReadTrueSVfile(string file, const map<string,int>& RefTable, vector<INV_t>& INV, vector<INS_t>& INS){
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
        if(strs[0]=="DEL"){
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
    }while(getline(input, line));
};

int main(int argc, char* argv[]){
    string filetype="VCF";
    string inputbam;
    string inputRefSV;
    string inputPredSV;
    string outputFile;
    int _FilterDeletion=1;
    int _FilterLowQual=0;

    int opt;
    while((opt=getopt(argc,argv,"b:r:p:o:d:l:t:"))!=EOF){
        switch(opt){
            case 'b': inputbam=(string)optarg; break;
            case 'r': inputRefSV=(string)optarg; break;
            case 'p': inputPredSV=(string)optarg; break;
            case 'o': outputFile=(string)optarg; break;
            case 'd': _FilterDeletion=atoi(optarg); break;
            case 'l': _FilterLowQual=atoi(optarg); break;
            case 't': filetype=(string)optarg; break;
        }
    }
    map<string, int> RefTable;
    vector<string> RefName;
    vector<int> RefLength;
    BuildRefName(inputbam, RefName, RefTable, RefLength);

    vector<INV_t> INV;
    vector<INS_t> INS;
    ReadTrueSVfile(inputRefSV, RefTable, INV, INS);

    SVset_t SV(filetype, inputPredSV, RefTable, _FilterLowQual);
    vector<int> SVHitID=SV.VerifyBPs(RefLength, INV, INS);
    double Precision=SV.CalPrecision(SVHitID, _FilterDeletion);
    double Sensitivity=SV.CalSensitivity(SVHitID, INV, INS);
    cout<<"Precision = "<<Precision<<"\tSensitivity = "<<Sensitivity<<endl;
    SV.Output(outputFile, SVHitID, Precision, Sensitivity);
}
