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
#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include "api/BamWriter.h"
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace BamTools;

string inputbam="";
string inputchimbam="";
string prefix="Merged";

class NameStrand{
public:
	static bool Smaller(pair<string,bool> a, pair<string,bool> b){
		if(a.first!=b.first)
			return a.first<b.first;
		else
			return a.second<b.second;
	};
	static bool Equal(pair<string,bool> a, pair<string,bool> b){
		return a.first==b.first && a.second==b.second;
	};
};

void MergeBAM(string inputbam, string inputchimbam, string outputfile, vector< pair<string,bool> >& SplitName){
	vector<string> ChimName;
	SplitName.clear();

	BamReader bamreader;
	bamreader.Open(inputchimbam);
	SamHeader header=bamreader.GetHeader();
	RefVector refvector=bamreader.GetReferenceData();

	BamWriter bamwriter;
	bamwriter.Open(outputfile, header, refvector);
	
	if(bamreader.IsOpen()){
		BamAlignment record;
        while(bamreader.GetNextAlignment(record)){
        	bamwriter.SaveAlignment(record);
        	ChimName.push_back(record.Name);
        	if(!record.IsPrimaryAlignment())
        		SplitName.push_back(make_pair(record.Name, record.IsFirstMate()));
        }
        bamreader.Close();
	}
	sort(ChimName.begin(), ChimName.end());
	vector<string>::iterator endit=unique(ChimName.begin(), ChimName.end());
	ChimName.resize(distance(ChimName.begin(), endit));

	sort(SplitName.begin(), SplitName.end(), NameStrand::Smaller);
	vector< pair<string, bool> >::iterator endit2=unique(SplitName.begin(), SplitName.end(), NameStrand::Equal);
	SplitName.resize(distance(SplitName.begin(), endit2));

	bamreader.Open(inputbam);
	if(bamreader.IsOpen()){
		BamAlignment record;
        while(bamreader.GetNextAlignment(record)){
        	if(record.MapQuality==255){
        		if(!binary_search(ChimName.begin(), ChimName.end(), record.Name))
        			bamwriter.SaveAlignment(record);
        	}
        }
        bamreader.Close();
    }
    bamwriter.Close();
};

void FindSplitBAM(string inputbam, string inputchimbam, string outputfile, vector< pair<string,bool> >& SplitName){
	vector< pair<string, bool> > ChimName;

	BamReader bamreader;
	bamreader.Open(inputchimbam);
	SamHeader header=bamreader.GetHeader();
	RefVector refvector=bamreader.GetReferenceData();

	BamWriter bamwriter;
	bamwriter.Open(outputfile, header, refvector);

	if(bamreader.IsOpen()){
		BamAlignment record;
        while(bamreader.GetNextAlignment(record)){
        	if(binary_search(SplitName.begin(), SplitName.end(), make_pair(record.Name, record.IsFirstMate()), NameStrand::Smaller)){
        		bamwriter.SaveAlignment(record);
        		ChimName.push_back(make_pair(record.Name, record.IsFirstMate()));
        	}
        }
        bamreader.Close();
	}
	sort(ChimName.begin(), ChimName.end(), NameStrand::Smaller);
	vector< pair<string, bool> >::iterator endit=unique(ChimName.begin(), ChimName.end(), NameStrand::Equal);
	ChimName.resize(distance(ChimName.begin(), endit));

	bamreader.Open(inputbam);
	if(bamreader.IsOpen()){
		BamAlignment record;
        while(bamreader.GetNextAlignment(record)){
        	bool flag1=binary_search(SplitName.begin(), SplitName.end(), make_pair(record.Name, record.IsFirstMate()), NameStrand::Smaller);
        	bool flag2=binary_search(ChimName.begin(), ChimName.end(), make_pair(record.Name, record.IsFirstMate()), NameStrand::Smaller);
        	if(flag1 && !flag2)
        		bamwriter.SaveAlignment(record);
        }
        bamreader.Close();
    }
    bamwriter.Close();
};

int main(int argc, char* argv[]){
	vector< pair<string, bool> > SplitName;

	int opt;
    while((opt=getopt(argc,argv,"b:c:p:"))!=EOF){
        switch(opt){
            case 'b': inputbam=(string)optarg; break;
            case 'c': inputchimbam=(string)optarg; break;
            case 'p': prefix=(string)optarg; break;
        }
    }
    size_t stop=inputbam.find_last_of("/");
    prefix=inputbam.substr(0,stop+1)+prefix;

	MergeBAM(inputbam, inputchimbam, prefix+".bam", SplitName);
	//FindSplitBAM(inputbam, inputchimbam, prefix+".splitters.bam", SplitName);
}