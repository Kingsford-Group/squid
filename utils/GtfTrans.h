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

using namespace std;

class Transcript_t{
public:
	string TransID;
	string GeneName;
	string Chr;
	char Strand;
	int TxStart;
	int TxEnd;
	vector<Interval_t> vExon;
	vector<Interval_t> vCDS;
public:
	Transcript_t(){};
	Transcript_t(string TransID, string GeneName, string Chr, char Strand, int TxStart, int TxEnd): TransID(TransID), GeneName(GeneName), Chr(Chr), Strand(Strand), TxStart(TxStart), TxEnd(TxEnd) {};
	Transcript_t(string TransID, string GeneName, string Chr, char Strand, int TxStart, int TxEnd, vector<Interval_t> vExon, vector<Interval_t> vCDS): 
		TransID(TransID), GeneName(GeneName), Chr(Chr), Strand(Strand), TxStart(TxStart), TxEnd(TxEnd), vExon(vExon), vCDS(vCDS) {};
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
};

void ReadGTF(string filename, vector<Transcript_t>& vTrans);

#endif