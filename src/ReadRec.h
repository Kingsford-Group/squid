#ifndef __ReadREC_H__
#define __ReadREC_H__

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
#include "SingleBamRec.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include "api/BamWriter.h"
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace BamTools;

const uint8_t BAM_XS_Forward='+';
const uint8_t BAM_XS_Reverse='-';

extern bool PhredType; // 0 for phred33 (default), 1 for phred64
extern int LowPhredLenThresh;

class SegmentGraph_t;

class ReadRec_t{
public:
    string Qname;
    vector<SingleBamRec_t> FirstRead, SecondMate;
    int FirstTotalLen, SecondTotalLen;
    bool FirstLowPhred, SecondLowPhred;
    bool MultiFilter;
public:
    ReadRec_t(){};
    ReadRec_t(BamAlignment record);
    bool operator < (const ReadRec_t& rhs) const {return Qname<rhs.Qname;};
    static bool FrontSmallerThan(const ReadRec_t& lhs, const ReadRec_t& rhs);
    void SortbyReadPos();
    void FilterSplitRecord();
    bool IsSingleAnchored() const;
    bool IsDiscordant() const;
    bool IsEndDiscordant(bool _isfirst) const;
    bool IsPairDiscordant(bool needcheck=true) const;
    int ReadCoverageGap() const;
    int ReadCoverageGap();
    void ModifybyGraph(SegmentGraph_t& SegmentGraph, const vector< vector<int> >& Components, vector<int>& singleRead_Node, vector< pair<int, int> >& Node_NewChr);
    string Print();
};

typedef vector<ReadRec_t> SBamrecord_t;

void BuildRefName(string bamfile, vector<string>& RefName, std::map<string,int>& RefTable);

void BuildReference(string fafile, const std::map<string,int>& RefTable, vector<int>& RefLength, vector<string>& RefSequence);

void UpdateReference(const SegmentGraph_t& SegmentGraph, const vector< vector<int> >& Components, vector<int>& RefLength, vector<string>& RefSequence);

SBamrecord_t BuildBWASBamRecord(const std::map<string,int> & RefTable, string bamfile);

SBamrecord_t BuildMainSBamRecord(const std::map<string,int> & RefTable, string bamfile);

SBamrecord_t BuildChimericSBamRecord(SBamrecord_t& SBamrecord, const std::map<string,int> & RefTable, string bamfile);

int AlignmentStat(const SBamrecord_t& SBamrecord);
int AlignmentStat(const SBamrecord_t& SBamrecord, string outputfile);

#endif
