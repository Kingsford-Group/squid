#ifndef __CONFIG_H__
#define __CONFIG_H__

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

using namespace std;

// parameter from BAM
extern uint16_t ReadLen;

// parameters from user input
	// read BAM file
extern bool UsingSTAR;
extern bool Phred_Type;
extern uint16_t Max_LowPhred_Len;
extern uint8_t Min_Phred;
extern uint16_t Min_MapQual;
	// building genome segment graph
extern int Concord_Dist_Pos;
extern int Concord_Dist_Idx;
extern int Min_Edge_Weight;
extern double DiscordantRatio; /*multiplier of discordant edge weight, related to normal-tumor cell ratio*/
	// input and output
extern string Input_BAM;
extern string Input_Chim_BAM;
extern string Input_FASTA;
extern string Output_Prefix;
extern bool Print_Graph;
extern bool Print_Components_Ordering;
extern bool Print_Total_Ordering;
extern bool Print_Rearranged_Genome;

// parsing
int print_help();
bool parse_arguments(int argc, char* argv[]);

#endif