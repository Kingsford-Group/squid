#ifndef __BPNode_H__
#define __BPNode_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <limits>
#include <ctime>
#include <cmath>
#include <map>
#include <boost/algorithm/string.hpp>
#include "BPEdge.h"

using namespace std;

class Edge_t;

class Node_t{
public:
    int Chr;
    int Position;
    int Length;
    int Support;
    double AvgDepth;
    vector<Edge_t*> HeadEdges;
    vector<Edge_t*> TailEdges;
public:
    Node_t(){};
    Node_t(int Chr, int Position, int Length): Chr(Chr), Position(Position), Length(Length){Support=0; AvgDepth=0;};
    Node_t(int Chr, int Position, int Length, int Support): Chr(Chr), Position(Position), Length(Length), Support(Support){AvgDepth=0;};
    Node_t(int Chr, int Position, int Length, int Support, double AvgDepth): Chr(Chr), Position(Position), Length(Length), Support(Support), AvgDepth(AvgDepth){};
    bool operator < (const Node_t& rhs) const{
        if(Chr!=rhs.Chr)
            return Chr<rhs.Chr;
        else if(Position!=rhs.Position)
            return Position<rhs.Position;
        else
            return Length<rhs.Length;
    };
    bool operator == (const Node_t& rhs) const{
        return ((Chr==rhs.Chr) && (Position==rhs.Position) && (Length==rhs.Length));
    };
    string Print(){
        return "{"+to_string(Chr)+","+to_string(Position)+","+to_string(Length)+"}";
    };
    string Print(const vector<string>& RefName){
        return "{"+RefName[Chr]+","+to_string(Position)+","+to_string(Length)+"}";
    };
};

#endif
