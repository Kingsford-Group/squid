#ifndef __SIMPLESV_H__
#define __SIMPLESV_H__
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include "TRA.h"

using namespace std;

enum SimpleSVType {INS=0, INV=1, DEL=2};
static const char* SimpleSVTypeString[]={"INS", "INV", "DEL"};

//class TRA_t;
class SimpleSV_t{
public:
    int RefID;
    int StartPos;
    int EndPos; //for INS, it is the length of insertion
    SimpleSVType Type;
    int ID;
public:
    SimpleSV_t(){}
    SimpleSV_t(int RefID, int StartPos, int EndPos, SimpleSVType Type, int ID): RefID(RefID), StartPos(StartPos), EndPos(EndPos), Type(Type), ID(ID){};
    void DecideSwap();
    bool operator < (const SimpleSV_t& rhs) const {
        if(RefID!=rhs.RefID)
            return RefID<rhs.RefID;
        else if(StartPos!=rhs.StartPos)
            return StartPos<rhs.StartPos;
        else
            return EndPos<rhs.EndPos;
    };
    string Print(){
        string info="(RefID="+to_string(RefID)+", Start="+to_string(StartPos)+", End="+to_string(EndPos)+", "+SimpleSVTypeString[Type]+", ID="+to_string(ID)+")\n";
        return info;
    };
    pair<int,int> UpdatePoint(pair<int,int> BP);
    void UpdateSimpleSV(vector<SimpleSV_t>::iterator itssv);
    void UpdateTRA(vector<TRA_t>::iterator ittra);
    void EditnReverse(map<int, int>& RefLength);
};

#endif
