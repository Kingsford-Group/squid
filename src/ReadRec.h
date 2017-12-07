/*
Part of SQUID transcriptomic structural variation detector
(c) 2017 by  Cong Ma, Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

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
#include "boost/algorithm/string.hpp"

using namespace std;
using namespace BamTools;

extern uint16_t ReadLen;
extern bool Phred_Type;
extern uint16_t Max_LowPhred_Len;
extern uint8_t Min_Phred;
extern uint16_t Min_MapQual;

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
    static bool Equal(const ReadRec_t& lhs, const ReadRec_t& rhs);
    void SortbyReadPos();
    void FilterSplitRecord();
    bool IsSingleAnchored() const;
    bool IsDiscordant() const;
    bool IsEndDiscordant(bool _isfirst) const;
    bool IsPairDiscordant(bool needcheck=true) const;
    int ReadCoverageGap() const;
    int ReadCoverageGap();
    string Print();
};

typedef vector<ReadRec_t> SBamrecord_t;

void BuildRefName(string bamfile, vector<string>& RefName, std::map<string,int>& RefTable, vector<int>& RefLength);

bool BuildRefSeq(string fafile, const std::map<string,int>& RefTable, vector<int>& RefLength, vector<string>& RefSequence);

void BuildChimericSBamRecord(SBamrecord_t& SBamrecord, const vector<string>& RefName, string bamfile);

int AlignmentStat(const SBamrecord_t& SBamrecord);
int AlignmentStat(const SBamrecord_t& SBamrecord, string outputfile);

#endif
