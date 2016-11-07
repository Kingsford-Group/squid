#ifndef __WRITEIO_H__
#define __WRITEIO_H__

#include "SegmentGraph.h"

using namespace std;

extern int concorddis, concorddisl, concordidx, concordidxl;

void WriteComponents(string outputfile, vector< vector<int> > Components);
void WriteBEDPE(string outputfile, SegmentGraph_t& SegmentGraph, vector< vector<int> >& Components, vector< pair<int, int> >& Node_NewChr, vector<string>& RefName);

#endif