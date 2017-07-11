#ifndef __SINGLEBAMREC_H__
#define __SINGLEBAMREC_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include <limits>
#include <ctime>
#include <cmath>
#include <iomanip>

using namespace std;

class SegmentGraph_t;

class SingleBamRec_t{
public:
    int32_t RefID;
    int32_t RefPos;
    int32_t ReadPos;
    int32_t MatchRead;
    int32_t MatchRef;
    uint8_t MapQual;
    bool IsReverse;
    bool IsFirstRead;
public:
    SingleBamRec_t(){};
    SingleBamRec_t(int32_t RefID, int32_t RefPos, int32_t ReadPos, int32_t MatchRef, int32_t MatchRead, uint16_t MapQual, bool IsReverse, bool IsFirstRead):
        RefID(RefID), RefPos(RefPos), ReadPos(ReadPos), MatchRef(MatchRef), MatchRead(MatchRead), MapQual(MapQual), IsReverse(IsReverse), IsFirstRead(IsFirstRead){};
    bool operator < (const SingleBamRec_t& rhs) const{
        if(RefID!=rhs.RefID)
            return RefID<rhs.RefID;
        else
            return RefPos<rhs.RefPos;
    };
    bool operator > (const SingleBamRec_t& rhs) const{
        if(RefID!=rhs.RefID)
            return RefID>rhs.RefID;
        else
            return RefPos>rhs.RefPos;
    };
    bool operator == (const SingleBamRec_t& rhs) const{
        return (RefID==rhs.RefID && RefPos==rhs.RefPos);
    };
    bool Same(const SingleBamRec_t& rhs) const{
        return (RefID==rhs.RefID && RefPos==rhs.RefPos && ReadPos==rhs.ReadPos && MatchRead==rhs.MatchRead && MatchRef==rhs.MatchRef && IsReverse==rhs.IsReverse && IsFirstRead==rhs.IsFirstRead);
    };
    static bool CompReadPos (const SingleBamRec_t& lhs, const SingleBamRec_t& rhs){
        return lhs.ReadPos<rhs.ReadPos;
    };
    void ModifybyGraph_Single(SegmentGraph_t& SegmentGraph, const vector< vector<int> >& Components, int nodeidx, vector< pair<int, int> >& Node_NewChr); //see ReadRec.cpp
};

#endif
