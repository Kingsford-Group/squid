#ifndef __BP_H__
#define __BP_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>

using namespace std;

class BP_t{
public:
	string Chr;
	int StartPos;
	int EndPos;
	bool IsLeft;
public:
	BP_t(){}
	BP_t(string Chr, int StartPos, int EndPos, bool IsLeft): Chr(Chr), StartPos(StartPos), EndPos(EndPos), IsLeft(IsLeft){};
	bool operator < (const BP_t& rhs) const{
		if(Chr!=rhs.Chr)
			return Chr<rhs.Chr;
		else if(StartPos!=rhs.StartPos)
			return StartPos<rhs.StartPos;
		else if(EndPos!=rhs.EndPos)
			return EndPos<rhs.EndPos;
		else
			return IsLeft<rhs.IsLeft;
	};
	bool operator == (const BP_t& rhs) const{
		return (Chr==rhs.Chr && StartPos==rhs.StartPos && EndPos==rhs.EndPos && IsLeft==rhs.IsLeft);
	};
    bool operator != (const BP_t& rhs) const{
        return !(*this==rhs);
    };
	string Print(){
		return Chr+"\t"+to_string(StartPos)+"\t"+to_string(EndPos)+"\t"+to_string(IsLeft);
	};
};

#endif