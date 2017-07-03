#ifndef __BPEDGE_H__
#define __BPEDGE_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <limits>
#include <ctime>
#include <cmath>
#include <map>
#include "boost/algorithm/string.hpp"
#include "BPNode.h"

using namespace std;

class Edge_t{
public:
    int Ind1, Ind2;
    bool Head1, Head2;
    int Weight, GroupWeight;
public:
    Edge_t(){};
    Edge_t(int _Ind1, bool _Head1, int _Ind2, bool _Head2){
        Weight=1; GroupWeight=0;
        if(_Ind1>_Ind2){
            Ind1=_Ind2; Head1=_Head2;
            Ind2=_Ind1; Head2=_Head1;
        }
        else{
            Ind1=_Ind1; Head1=_Head1;
            Ind2=_Ind2; Head2=_Head2;
        }
    };
    Edge_t(int _Ind1, bool _Head1, int _Ind2, bool _Head2, int Weight): Weight(Weight) {
        GroupWeight=0;
        if(_Ind1>_Ind2){
            Ind1=_Ind2; Head1=_Head2;
            Ind2=_Ind1; Head2=_Head1;
        }
        else{
            Ind1=_Ind1; Head1=_Head1;
            Ind2=_Ind2; Head2=_Head2;
        }
    };
    static bool WeakEqual(const Edge_t& lhs, const Edge_t& rhs){
        if((lhs.Ind1==rhs.Ind1 && lhs.Ind2==rhs.Ind2) || (lhs.Ind1==rhs.Ind2 && lhs.Ind2==rhs.Ind1))
            return true;
        else
            return false;
    };
    bool operator < (const Edge_t& rhs) const{
        if(Ind1!=rhs.Ind1)
            return Ind1<rhs.Ind1;
        else if(Ind2!=rhs.Ind2)
            return Ind2<rhs.Ind2;
        else if(Head1!=rhs.Head1)
            return (int)Head1<(int)rhs.Head1;
        else if(Head2!=rhs.Head2)
            return (int)Head2<(int)rhs.Head2;
        else
            return false;
    };
    bool operator == (const Edge_t& rhs) const{
        if(Ind1==rhs.Ind1 && Ind2==rhs.Ind2 && Head1==rhs.Head1 && Head2==rhs.Head2)
            return true;
        else
            return false;
    }
};

#endif
