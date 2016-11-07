#ifndef __TRA_H__
#define __TRA_H__
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/algorithm/string.hpp>
//#include "SimpleSV.h"

using namespace std;

enum DirType {left=0, right=1}; //which part is kept in the same chr; DT1 is infered from supported reads, DT2 is the opposite of supported reads

class SimpleSV_t;
class TRA_t{
public:
    int Ref1, Ref2;
    int Pos1, Pos2;
    DirType DT1, DT2;
    int ID;
public:
    TRA_t(){}
    TRA_t(int Ref1, int Pos1, DirType DT1, int Ref2, int Pos2, DirType DT2, int ID):Ref1(Ref1), Pos1(Pos1), DT1(DT1), Ref2(Ref2), Pos2(Pos2), DT2(DT2), ID(ID){}
    void DecideSwap();
    bool operator < (const TRA_t& rhs) const{
        if(Ref1!=rhs.Ref1)
            return Ref1<rhs.Ref1;
        else if(Pos1!=rhs.Pos1)
            return Pos1<rhs.Pos1;
        else if(Ref2!=rhs.Ref2)
            return Ref2<rhs.Ref2;
        else
            return Pos2<rhs.Pos2;
    };
    string Print(){
        string info="(Ref1="+to_string(Ref1)+", Pos1="+to_string(Pos1)+", DirType="+to_string(DT1)+", Ref2="+to_string(Ref2)+", Pos1="+to_string(Pos2)+", DirType="+to_string(DT2)+", ID="+to_string(ID)+")\n";
        return info;
    };
    pair<int,int> UpdatePoint(map<int, int>& RefLength, pair<int,int> BP, DirType DT);
    void UpdateSimpleSV(map<int, int>& RefLength, vector<SimpleSV_t>::iterator itssv);
    void UpdateTRA(map<int, int>& RefLength, vector<TRA_t>::iterator ittra);
    void EditnReverse(map<int, int>& RefLength);
};

#endif
