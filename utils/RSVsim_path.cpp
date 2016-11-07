#include "BPNode.h"

using namespace std;

static map<char,char> Nucleotide={{'A','T'},{'C','G'},{'G','C'},{'T','A'},{'R','Y'},{'Y','R'},{'S','W'},{'W','S'},{'K','M'},{'M','K'},{'B','V'},{'V','B'},{'D','H'},{'H','D'}, {'N','N'}, {'.','.'},{'-','-'}};
void ReverseComplement(string::iterator itbegin, string::iterator itend){
    for(string::iterator it=itbegin; it!=itend; it++)
        *it=Nucleotide[toupper(*it)];
    std::reverse(itbegin, itend);
};

void BuildReference(string fafile, map<string,int>& RefTable, vector<string>& RefName, vector<int>& RefLength, vector<string>& RefSequence){
    ifstream input(fafile);
    string line, tmpseq;
    int curlength=0, count=0;
    while(getline(input, line)){
        if(line[0]=='>'){
            if(curlength!=0){
                RefLength.push_back(curlength);
                RefSequence.push_back(tmpseq);
            }
            curlength=0;
            tmpseq.clear();
            size_t start=line.find_first_of(" \t");
            string chrname=line.substr(1, start-1);
            RefTable[chrname]=count++;
            RefName.push_back(chrname);
        }
        else{
            curlength+=line.size();
            for(int i=0; i<line.size(); i++)
                line[i]=toupper(line[i]);
            tmpseq+=line;
        }
    }
    RefLength.push_back(curlength);
    RefSequence.push_back(tmpseq);
    input.close();
};

void ReadGenome(string fafile, const map<string,int>& RefTable, vector<string>& RefSequence){
    ifstream input(fafile);
    string line, tmpseq;
    int prevchr=0;
    RefSequence.resize(RefTable.size());
    while(getline(input, line)){
        if(line[0]=='>'){
            if(tmpseq.size()!=0){
                RefSequence[prevchr]=tmpseq;
            }
            tmpseq.clear();
            size_t start=line.find_first_of(" \t");
            string chrname=line.substr(1, start-1);
            prevchr=RefTable.at(chrname);
        }
        else{
            for(int i=0; i<line.size(); i++)
                line[i]=toupper(line[i]);
            tmpseq+=line;
        }
    }
    RefSequence[prevchr]=tmpseq;
    input.close();
};

void BuildNodes(vector<string> rsvsimfiles, vector<Node_t>& vNodes, const map<string, int>& RefTable, const vector<int>& RefLength){
    vNodes.clear();
    vector< vector<int> > breakpoints; breakpoints.resize(RefTable.size());
    for(int i=0; i<rsvsimfiles.size(); i++){
        ifstream input(rsvsimfiles[i]);
        string line;
        getline(input, line);
        vector<string> strs;
        boost::split(strs,line,boost::is_any_of("\t"));
        if(strs.size()<=7){ //Deletions, Inversion, Duplication
            while(getline(input, line)){
                boost::split(strs,line,boost::is_any_of("\t"));
                if(strs[5]=="1")
                    continue;
                size_t start=strs[1].find_first_of(" ");
                string chrname=strs[1].substr(0, start);
                int chr=RefTable.at(chrname);
                breakpoints[chr].push_back(stoi(strs[2])-1);
                breakpoints[chr].push_back(stoi(strs[3]));
            }
        }
        else if(strs.size()==12 && strs[8]=="Copied"){ //Insertions, which cuts or copies from chrA and paste to chr B, aka transposon. =deletion+insertion
            while(getline(input, line)){
                boost::split(strs,line,boost::is_any_of("\t"));
                size_t start=strs[1].find_first_of(" ");
                string chrname=strs[1].substr(0, start);
                int chr=RefTable.at(chrname);
                breakpoints[chr].push_back(stoi(strs[2])-1);
                breakpoints[chr].push_back(stoi(strs[3]));
                start=strs[4].find_first_of(" ");
                chrname=strs[4].substr(0, start);
                chr=RefTable.at(chrname);
                breakpoints[chr].push_back(stoi(strs[5])-1);
            }
        }
        else if(strs.size()==12 && strs[9]=="Balanced"){ //Translocations
            while(getline(input, line)){
                boost::split(strs,line,boost::is_any_of("\t"));
                size_t start=strs[1].find_first_of(" ");
                string chrname=strs[1].substr(0, start);
                int chr=RefTable.at(chrname);
                breakpoints[chr].push_back((strs[2]=="1")?stoi(strs[3]):stoi(strs[2])-1);
                start=strs[5].find_first_of(" ");
                chrname=strs[5].substr(0, start);
                chr=RefTable.at(chrname);
                breakpoints[chr].push_back((strs[6]=="1")?stoi(strs[7]):stoi(strs[6])-1);
            }
        }
        input.close();
    }
    for(int i=0; i<breakpoints.size(); i++)
        sort(breakpoints[i].begin(), breakpoints[i].end());
    for(int i=0; i<breakpoints.size(); i++){
        for(int j=0; j<breakpoints[i].size(); j++){
            if(breakpoints[i][j]==0)
                continue;
            else{
                Node_t tmp(i, breakpoints[i][j-1], breakpoints[i][j]-breakpoints[i][j-1]);
                vNodes.push_back(tmp);
            }
        }
        if(vNodes.back().Position+vNodes.back().Length<RefLength[i]){
            Node_t tmp(i, vNodes.back().Position+vNodes.back().Length, RefLength[i]-vNodes.back().Position-vNodes.back().Length);
            vNodes.push_back(tmp);
        }
    }
};

void BuildPaths(vector<string> rsvsimfiles, vector<Node_t>& vNodes, vector< vector<int> >& Paths, const map<string, int>& RefTable){
    Paths.clear();
    vector<int> tmp;
    for(int i=0; i<vNodes.size(); i++){
        if(tmp.size()==0 || vNodes[i].Chr==vNodes[tmp.back()-1].Chr){
            tmp.push_back(i+1);
        }
        else{
            Paths.push_back(tmp);
            tmp.clear();
            tmp.push_back(i+1);
        }
    }
    Paths.push_back(tmp);
    vector<int>::iterator itpath;
    for(int i=0; i<rsvsimfiles.size(); i++){
        ifstream input(rsvsimfiles[i]);
        string line;
        getline(input, line);
        vector<string> strs;
        boost::split(strs,line,boost::is_any_of("\t"));
        if(strs.size()==6){ //Deletions
            while(getline(input, line)){
                boost::split(strs,line,boost::is_any_of("\t"));
                size_t start=strs[1].find_first_of(" ");
                string chrname=strs[1].substr(0, start);
                int chr=RefTable.at(chrname), position=stoi(strs[2])-1;
                for(itpath=Paths[chr].begin(); itpath!=Paths[chr].end(); itpath++)
                    if(vNodes[abs(*itpath)-1].Position==position){
                        Paths[chr].erase(itpath); break;
                    }
            }
        }
        else if(strs.size()==7 && strs[5]=="Duplications"){ //Duplications, the first copy is normal not duplication, start from the second copy
            while(getline(input, line)){
                boost::split(strs,line,boost::is_any_of("\t"));
                size_t start=strs[1].find_first_of(" ");
                string chrname=strs[1].substr(0, start);
                int chr=RefTable.at(chrname), position=stoi(strs[2])-1, duptimes=stoi(strs[5])-1;
                if(duptimes!=0){
                    for(itpath=Paths[chr].begin(); itpath!=Paths[chr].end(); itpath++)
                        if(vNodes[abs(*itpath)-1].Position==position){
                            int index=(*itpath);
                            for(int j=0; j<duptimes; j++)
                                Paths[chr].insert(itpath, index);
                            break;
                        }
                }
            }
        }
        else if(strs.size()==7){ //Inversions
            while(getline(input,line)){
                boost::split(strs,line,boost::is_any_of("\t"));
                size_t start=strs[1].find_first_of(" ");
                string chrname=strs[1].substr(0, start);
                int chr=RefTable.at(chrname), position=stoi(strs[2])-1;
                for(itpath=Paths[chr].begin(); itpath!=Paths[chr].end(); itpath++)
                    if(vNodes[abs(*itpath)-1].Position==position){
                        (*itpath)=-abs(*itpath); break;
                    }
            }
        }
        else if(strs.size()==12 && strs[8]=="Copied"){ //Insertions, which cuts or copies from chrA and paste to chr B
            while(getline(input, line)){
                boost::split(strs,line,boost::is_any_of("\t"));
                size_t start=strs[1].find_first_of(" ");
                string chrname=strs[1].substr(0, start);
                int chr=RefTable.at(chrname), position=stoi(strs[2])-1, index;
                for(itpath=Paths[chr].begin(); itpath!=Paths[chr].end(); itpath++)
                    if(vNodes[abs(*itpath)-1].Position==position){
                        index=abs(*itpath);
                        if(strs[8]=="FALSE")
                            Paths[chr].erase(itpath);
                        break;
                    }
                start=strs[4].find_first_of(" ");
                chrname=strs[4].substr(0, start);
                chr=RefTable.at(chrname); position=stoi(strs[5])-1;
                for(itpath=Paths[chr].begin(); itpath!=Paths[chr].end(); itpath++)
                    if(vNodes[abs(*itpath)-1].Position==position){
                        Paths[chr].insert(itpath, index);
                        break;
                    }
            }
        }
        else if(strs.size()==12 && strs[9]=="Balanced"){
            while(getline(input, line)){
                boost::split(strs,line,boost::is_any_of("\t"));
                size_t start=strs[1].find_first_of(" ");
                string chrname=strs[1].substr(0, start);
                int chr=RefTable.at(chrname), position=stoi(strs[2])-1, index;
                vector<int>::iterator itpath2;
                for(itpath=Paths[chr].begin(); itpath!=Paths[chr].end(); itpath++)
                    if(vNodes[abs(*itpath)-1].Position==position)
                        break;
                start=strs[5].find_first_of(" ");
                chrname=strs[5].substr(0, start);
                chr=RefTable.at(chrname); position=stoi(strs[6])-1;
                for(itpath2=Paths[chr].begin(); itpath2!=Paths[chr].end(); itpath2++)
                    if(vNodes[abs(*itpath2)-1].Position==position)
                        break;
                if((strs[2]=="1")==(strs[6]=="1")){
                    index=(*itpath);
                    (*itpath)=(*itpath2); (*itpath2)=index;
                }
                else{
                    index=(*itpath);
                    (*itpath)=-(*itpath2); (*itpath2)=-index;
                }
            }
        }
        input.close();
    }
};

bool CompareLength(vector<string>& RefSequence0, vector<string>& RefSequence1, vector<Node_t>& vNodes, vector< vector<int> >& Paths){
    if(RefSequence1.size()!=Paths.size()){
        cout<<"number of chromosomes is inconsistant\n";
        return false;
    }
    else{
        bool correct=true;
        for(int i=0; i<RefSequence1.size(); i++){
            int curlength=0;
            for(int j=0; j<Paths[i].size(); j++)
                curlength+=vNodes[abs(Paths[i][j])-1].Length;
            if(curlength!=RefSequence1[i].size()){
                correct=false;
                cout<<i<<'\t'<<curlength<<'\t'<<"False\n";
            }
            else
                cout<<i<<'\t'<<curlength<<'\t'<<"True\n";
        }
        return correct;
    }
};

bool CompareSequence(vector<string>& RefSequence0, vector<string>& RefSequence1, vector<Node_t>& vNodes, vector< vector<int> >& Paths){
    if(RefSequence1.size()!=Paths.size()){
        cout<<"number of chromosomes is inconsistant\n";
        return false;
    }
    else{
        bool correct=true;
        for(int i=0; i<Paths.size(); i++){
            int curposition=0;
            for(int j=0; j<Paths[i].size(); j++){
                if(Paths[i][j]>0){
                    string oldseq=RefSequence0[vNodes[Paths[i][j]-1].Chr].substr(vNodes[Paths[i][j]-1].Position, vNodes[Paths[i][j]-1].Length);
                    string newseq=RefSequence1[i].substr(curposition, vNodes[Paths[i][j]-1].Length);
                    if(oldseq!=newseq){
                        correct=false;
                        cout<<i<<'\t'<<j<<'\t'<<"False\n";
                    }
                }
                else{
                    string oldseq=RefSequence0[vNodes[-Paths[i][j]-1].Chr].substr(vNodes[-Paths[i][j]-1].Position, vNodes[-Paths[i][j]-1].Length);
                    string newseq=RefSequence1[i].substr(curposition, vNodes[-Paths[i][j]-1].Length);
                    ReverseComplement(oldseq.begin(), oldseq.end());
                    if(oldseq!=newseq){
                        correct=false;
                        cout<<i<<'\t'<<j<<'\t'<<"False\n";
                    }
                }
                curposition+=vNodes[abs(Paths[i][j])-1].Length;
            }
        }
        return correct;
    }
};

void OutputPaths(string outputfile, const vector<string>&RefName, vector<Node_t>& vNodes, vector< vector<int> >& Paths){
    ofstream output(outputfile, ios::out);
    for(int i=0; i<Paths.size(); i++){
        output<<i<<'\t';
        int j;
        for(j=0; j<Paths[i].size()-1; j++){
            output<<"{"<<RefName[vNodes[abs(Paths[i][j])-1].Chr]<<","<<vNodes[abs(Paths[i][j])-1].Position<<","<<vNodes[abs(Paths[i][j])-1].Length<<"}";
            if(Paths[i][j]>0)
                output<<"F-";
            else
                output<<"R-";
        }
        output<<"{"<<RefName[vNodes[abs(Paths[i][j])-1].Chr]<<","<<vNodes[abs(Paths[i][j])-1].Position<<","<<vNodes[abs(Paths[i][j])-1].Length<<"}";
        if(Paths[i][j]>0)
            output<<"F\n";
        else
            output<<"R\n";
    }
    output.close();
};

int main(int argc, char* argv[]){
    vector<string> rsvsimfiles;
    for(int i=3; i<argc-1; i++){
        string tmp(argv[i]);
        rsvsimfiles.push_back(tmp);
    }
    map<string, int> RefTable;
    vector<string> RefName;
    vector<int> RefLength;
    vector<string> RefSequence0;
    BuildReference(argv[1], RefTable, RefName, RefLength, RefSequence0);

    vector<string> RefSequence1;
    ReadGenome(argv[2], RefTable, RefSequence1);

    vector<Node_t> vNodes;
    vector< vector<int> > Paths;
    BuildNodes(rsvsimfiles, vNodes, RefTable, RefLength);
    BuildPaths(rsvsimfiles, vNodes, Paths, RefTable);

    CompareLength(RefSequence0, RefSequence1, vNodes, Paths);
    CompareSequence(RefSequence0, RefSequence1, vNodes, Paths);
    OutputPaths(argv[argc-1], RefName, vNodes, Paths);
}
