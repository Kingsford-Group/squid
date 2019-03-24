/*
Part of SQUID transcriptomic structural variation detector
(c) 2017 by  Cong Ma, Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "SingleBamRec.h"
#include "ReadRec.h"
#include "BPNode.h"
#include "BPEdge.h"
#include "SegmentGraph.h"
#include "WriteIO.h"
#include "Config.h"

using namespace std;

int main(int argc, char* argv[]){
	time_t CurrentTime;
	string CurrentTimeStr;
	
	bool success=parse_arguments(argc, argv);

	if(success){
		map<string, int> RefTable;
		vector<string> RefName;
		vector<int> RefLength;

		BuildRefName(Input_BAM, RefName, RefTable, RefLength);
		for(map<string,int>::iterator it=RefTable.begin(); it!=RefTable.end(); it++)
			cout<<"Reference name "<<it->first<<"\t-->\t"<<it->second<<endl;

		SBamrecord_t Chimrecord;
		if(Input_Chim_BAM.size()!=0){
			BuildChimericSBamRecord(Chimrecord, RefName, Input_Chim_BAM);
		}

		SegmentGraph_t SegmentGraph(RefLength, Chimrecord, Input_BAM);
		SegmentGraph.ProcessFilter_step_edge();
		SegmentGraph.ProcessFilter_readnames(Chimrecord);

		SegmentGraph.OutputGraph(Output_Prefix + "_origraph.txt");
		SegmentGraph.OutputFilteredEdgeSupport(Output_Prefix + "_filtered_edges.txt");
		
		time(&CurrentTime);
		CurrentTimeStr=ctime(&CurrentTime);
		cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] Done."<<endl;

	}
}
