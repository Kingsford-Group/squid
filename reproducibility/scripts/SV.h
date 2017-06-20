#ifndef __SV_H__
#define __SV_H__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include "SimpleSV.h"
#include "TRA.h"

using namespace std;

struct SV_t{
public:
    vector<SimpleSV_t> vSimpleSV;
    vector<TRA_t> vTRA;
    vector<bool> vSimpleReads;
    vector<bool> vTRAReads;
    bool updated;
    int num;
public:
    SV_t(){};
    SV_t(map<string,int>& RefTable, int numfile, string files[]);
    void Update2NewPos(map<int,int>& RefLength);
    static bool tupleCompare(const std::tuple<int,int,int,int>& lhs, const std::tuple<int,int,int,int>& rhs) {
        if(get<2>(lhs)!=get<2>(rhs))
            return get<2>(lhs)<get<2>(rhs);
        else
            return get<3>(lhs)<get<3>(rhs);
    }
    static bool tupleEqual(const std::tuple<int,int,int,int>& lhs, const std::tuple<int,int,int,int>& rhs) {
        if(get<0>(lhs)!=get<0>(rhs) || get<1>(lhs)!=get<1>(rhs))
            return false;
        else
            return true;
    }
    vector< tuple<int,int,int,int> > UpdateNextRound2new(SV_t svnext, map<int,int>& RefLength);
    vector< tuple<int,int,int,int> > UpdateNextRound2old(SV_t svnext, map<int,int>& RefLength);
    void WriteNextRoundBPPos(vector< tuple<int,int,int,int> >& BreakPoint, char* outputfile);
    void WriteNextRoundBPPos(vector<string>& RefName, vector< tuple<int,int,int,int> >& BreakPoint, char* outputfile);
    void WritenewSVPos(vector<string>& RefName, map<int,int>& RefLength, string outputfile);
    void WriteoldSVPos(vector<string>& RefName, map<int,int>& RefLength, string outputfile);
    void WriteFilterednewSVPos(vector<string>& RefName, map<int,int>& RefLength, string outputfile);
    void WriteFilteredoldSVPos(vector<string>& RefName, map<int,int>& RefLength, string outputfile);
    bool IsIntersect(SV_t svrhs);
    void IsCoveredbyReads(string bedfile, map<string,int>& RefTable, map<int,int>& RefLength, vector<string>& TransName, bool is_old_coordicate, int threshold=2);
};

#endif
