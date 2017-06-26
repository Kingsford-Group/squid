#ifndef __SV_H__
#define __SV_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include "BP.h"

using namespace std;

class SV_t{
public:
	BP_t BP1, BP2;
public:
	SV_t(){};
	SV_t(BP_t _BP1, BP_t _BP2){
		if(_BP1<_BP2){
			BP1=_BP1;
			BP2=_BP2;
		}
		else{
			BP1=_BP2;
			BP2=_BP1;
		}
	};

	bool operator < (const SV_t& rhs) const{
		if(BP1!=rhs.BP1)
			return BP1<rhs.BP1;
		else
			return BP2<rhs.BP2;
	};
	bool operator == (const SV_t& rhs) const{
		return (BP1==rhs.BP1 && BP2==rhs.BP2);
	};
};

#endif