#ifndef __GTFINTERVAL_H__
#define __GTFINTERVAL_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>

using namespace std;

class Interval_t{
public:
	int Chr;
	int StartPos;
	int EndPos;
	string TransID;
	int Number;
	char Strand;
public:
	Interval_t(){};
	Interval_t(int Chr, int StartPos, int EndPos, char Strand, string TransID):Chr(Chr), StartPos(StartPos), EndPos(EndPos), Strand(Strand), TransID(TransID) {};
	Interval_t(int Chr, int StartPos, int EndPos, char Strand, string TransID, int Number):Chr(Chr), StartPos(StartPos), EndPos(EndPos), Strand(Strand), TransID(TransID), Number(Number) {};
	bool operator < (const Interval_t& rhs) const{
		if(Chr!=rhs.Chr)
			return Chr<rhs.Chr;
		if(StartPos!=rhs.StartPos)
			return StartPos<rhs.StartPos;
		if(EndPos!=rhs.EndPos)
			return EndPos<rhs.EndPos;
	};
	static bool TransIDComp(const Interval_t& lhs, const Interval_t& rhs) {
		return lhs.TransID<rhs.TransID;
	};
	static bool ExonNumberComp(const Interval_t& lhs, const Interval_t& rhs) {
		if(lhs.TransID!=rhs.TransID)
			return lhs.TransID<rhs.TransID;
		else
			return lhs.Number<rhs.Number;
	};
};

#endif