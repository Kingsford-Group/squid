#ifndef __BP_H__
#define __BP_H__

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
#include "../src/BPNode.h"

using namespace std;

enum SVCode {INV, DUP, DEL, TRA, OTH};
const string SVString[4]={"INV", "DUP", "DEL", "TRA", "OTH"};

typedef tuple<int,int,int,int> INV_t; // Chr, StartPos, EndPos, ID
typedef tuple<int,int,int,int,int,int> INS_t; // INSChr, INSPos, DELChr, DELStartPos, DELEndPos, ID
typedef tuple<int,int,bool,int,int,bool,int> TRA_t; //Chr1, Pos1, IsLeft1, Chr2, Pos2, IsLeft2, ID

class Breakpoint_t{
public:
    int Chr;
    int Position;
    bool IsLeft;
public:
    Breakpoint_t(){}
    Breakpoint_t(int Chr, int Position): Chr(Chr), Position(Position){IsLeft=false;}
    Breakpoint_t(int Chr, int Position, bool IsLeft): Chr(Chr), Position(Position), IsLeft(IsLeft){};
    bool operator < (const Breakpoint_t& rhs) const{ //why const is needed but not the rest?
        if(Chr!=rhs.Chr)
            return Chr<rhs.Chr;
        else if(Position!=rhs.Position)
            return Position<rhs.Position;
        else
            return IsLeft<rhs.IsLeft;
    };
    bool operator == (const Breakpoint_t& rhs) const{
        return (Chr==rhs.Chr && Position==rhs.Position);
    };
    string Print();
};

class SVset_t{
public:
    vector< tuple<Breakpoint_t, Breakpoint_t, SVCode, int, string> > vBP; // 2 breakpoints, SV type, quality score, ID in vcf
public:
    SVset_t(){};
    SVset_t(string filetype, string svfile, const map<string,int>& RefTable, int _FilterLowQual);
    void VCFRead(string svfile, const map<string,int>& RefTable, int _FilterLowQual);
    void BDRead(string svfile, const map<string,int>& RefTable, int _FilterLowQual);
    void BedpeRead(string svfile, const map<string,int>& RefTable, int _FilterLowQual);
    int CountafterFilter(int _FilterDeletion);
    vector<int> VerifyBPs(vector<int>& RefLength, vector<INV_t>& INV, vector<INS_t>& INS);
    double CalPrecision(vector<int>& SVHitID, int _FilterDeletion);
    double CalSensitivity(vector<int>& SVHitID, vector<INV_t>& INV, vector<INS_t>& INS);
    void Output(string outputfile);
    void Output(string outputfile, vector<int>& VerifyResult, double & Precision, double & Sensitivity);
};

#endif
