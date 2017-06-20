#ifndef __GTFTRANS_H__
#define __GTFTRANS_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include <boost/algorithm/string.hpp>
#include "GtfInterval.h"
#include "SV.h"

using namespace std;

struct AffectExonRelativePos{
	bool IsStartBP; // whether nearest exon start/end position is relative to start or end breakpoint
	bool IsLeft; // whether this exon is to the left or right of corresponding breakpoint
	int RelativePos;
	int svID;
	AffectExonRelativePos(){};
	AffectExonRelativePos(bool IsStartBP, bool IsLeft, int RelativePos,int svID): IsStartBP(IsStartBP), IsLeft(IsLeft), RelativePos(RelativePos), svID(svID){};

	bool operator < (const AffectExonRelativePos& rhs) const{
		if(IsStartBP!=rhs.IsStartBP)
			return IsStartBP<rhs.IsStartBP;
		else if(IsLeft!=rhs.IsLeft)
			return IsLeft<rhs.IsLeft;
		else
			return RelativePos<rhs.RelativePos; 
	};
	bool operator == (const AffectExonRelativePos& rhs) const{
		return (IsStartBP==rhs.IsStartBP && IsLeft==rhs.IsLeft && RelativePos==rhs.RelativePos);
	}
};

class Transcript_t{
public:
	string TransID;
	string GeneName;
	int Chr;
	char Strand;
	int TxStart;
	int TxEnd;
	vector<Interval_t> vExon;
public:
	Transcript_t(){};
	Transcript_t(int Chr, int TxStart, int TxEnd, char Strand, string TransID, string GeneName): Chr(Chr), TxStart(TxStart), TxEnd(TxEnd), Strand(Strand), TransID(TransID), GeneName(GeneName) {};
	Transcript_t(int Chr, int TxStart, int TxEnd, char Strand, string TransID, string GeneName, vector<Interval_t> vExon): 
		Chr(Chr), TxStart(TxStart), TxEnd(TxEnd), Strand(Strand), TransID(TransID), GeneName(GeneName), vExon(vExon) {};
	bool operator < (const Transcript_t& rhs) const{
		if(Chr!=rhs.Chr)
			return Chr<rhs.Chr;
		if(TxStart!=rhs.TxStart)
			return TxStart<rhs.TxStart;
		if(TxEnd!=rhs.TxEnd)
			return TxEnd<rhs.TxEnd;
	};
	static bool TransIDComp(const Transcript_t& lhs, const Transcript_t& rhs) {
		return lhs.TransID<rhs.TransID;
	};
	static string GetFeature(string line, string FeatureName);

	bool IsAffectedbySV(SimpleSV_t& sv, vector< pair<AffectExonRelativePos, AffectExonRelativePos> >& aepos);
	bool IsAffectedbySV(TRA_t& sv, vector< pair<AffectExonRelativePos, AffectExonRelativePos> >& aepos);
};

void ReadGTF(string filename, vector<Transcript_t>& vTrans, map<string,int>& RefTable);
void FilterGTF(string fluxpro, vector<Transcript_t>& vTrans);
void FilterGTF(vector<string>& TransName, vector<Transcript_t>& vTrans);
void FindAffected(vector<Transcript_t>& vTrans, SV_t& SVs, vector< vector< pair<AffectExonRelativePos, AffectExonRelativePos> > >& AffectedPos);
void WriteBreakpoint(string filename, SV_t& SVs, vector<string>& RefName, vector< vector< pair<AffectExonRelativePos, AffectExonRelativePos> > >& AffectedPos); // SVs is with updated coordinate, AffectedExonRelative is derived by non-updated ones

#endif