// module 1 (initial): generate fused genome junction for non-fusion-gene TSV
// module 2 (initial): generate temporary gtf file for STAR alignment
// module 3 (final): generate sequence for junction transcript

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
#include "GtfTrans.h"
#include "SV.h"

using namespace std;

void ReadBedpeSV(string filename, vector<SV_t>& SVs){
	SVs.clear();
	ifstream input(filename);
	string line;
	while(getline(input, line)){
		if(line[0]=='#')
			continue;
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t"));
		BP_t bp1(strs[0], stoi(strs[1]), stoi(strs[2]), (strs[8]=="-"));
		BP_t bp2(strs[3], stoi(strs[4]), stoi(strs[5]), (strs[9]=="-"));
		SV_t sv(bp1, bp2);
		SVs.push_back(sv);
	}
	input.close();
};

void ReadGenome(string filename, map<string, string>& Genome){
	ifstream input(filename);
	string line;
	string tmpseq="";
	string tmpname="";
	while(getline(input, line)){
		if(line[0]=='>'){
			if(tmpseq.size()!=0){
				Genome[tmpname]=tmpseq;
			}
			vector<string> strs;
			boost::split(strs, line, boost::is_any_of(" \t"));
			tmpseq="";
			tmpname=strs[0].substr(1);
		}
		else{
			tmpseq+=line;
		}
	}
	if(tmpseq.size()!=0){
		Genome[tmpname]=tmpseq;
	}
	input.close();
};

vector<Transcript_t> GeneBinarySearch(vector<Transcript_t>::iterator firstinterator, vector<Transcript_t>::iterator lastiterator, BP_t bp){
	vector<Transcript_t> genes;
	vector<Transcript_t>::iterator itbegin=firstinterator;
	vector<Transcript_t>::iterator itend=lastiterator;
	int bppos=(bp.IsLeft)?bp.StartPos:bp.EndPos;
	vector<Transcript_t>::iterator itmid;
	while(distance(itbegin, itend)>0){
		itmid=itbegin+distance(itbegin, itend)/2;
		if(itmid->Chr==bp.Chr && itmid->TxStart<=bppos && itmid->TxEnd>=bppos){
			itbegin=itmid;
			itend=itmid;
			break;
		}
		else if(itmid->Chr<bp.Chr || (itmid->Chr==bp.Chr && itmid->TxEnd<bppos))
			itbegin=itmid;
		else if(itmid->Chr>bp.Chr || (itmid->Chr==bp.Chr && itmid->TxStart>bppos))
			itend=itmid;
	}
	for(; distance(itbegin, itmid)<20 && distance(firstinterator, itbegin)>=0; itbegin--){
		if(itbegin->Chr==bp.Chr && itbegin->TxStart<=bppos && itbegin->TxEnd>=bppos)
			genes.push_back(*itbegin);
	}
	itend++;
	for(; distance(itend, itmid)<20 && distance(itend, lastiterator)>=0; itend++){
		if(itend->Chr==bp.Chr && itend->TxStart<=bppos && itend->TxEnd>=bppos)
			genes.push_back(*itend);
	}
	return genes;
};

void InitialJunction(vector<SV_t>& SVs, vector<Transcript_t>& vTrans, map<string, string>& Genome, string OutPrefix){
	int thresh=5;
	int countseq=0;
	ofstream fpfusion(OutPrefix+"_fusiongene.gtf", ios::out);
	ofstream fpseq(OutPrefix+"_juncseq.fa", ios::out);
	ofstream fptmpannot(OutPrefix+"_juncannot.gtf", ios::out);
	for(vector<SV_t>::iterator it=SVs.begin(); it!=SVs.end(); it++){
		vector<Transcript_t> genes1=GeneBinarySearch(vTrans.begin(), vTrans.end(), it->BP1);
		vector<Transcript_t> genes2=GeneBinarySearch(vTrans.begin(), vTrans.end(), it->BP2);
		bool is_fusion=false;
		if(genes1.size()!=0 && genes2.size()!=0){
			for(int i=0; i<genes1.size(); i++)
				for(int j=0; j<genes2.size(); j++){
					if((it->BP1.IsLeft==it->BP2.IsLeft)!=(genes1[i].Strand==genes2[j].Strand))
						is_fusion=true;
				}
		}
		// if fusion gene, output new transcript directly
		if(is_fusion){
			for(int i=0; i<genes1.size(); i++)
				for(int j=0; j<genes2.size(); j++){
					if((it->BP1.IsLeft==it->BP2.IsLeft)==(genes1[i].Strand==genes2[i].Strand))
						continue;
					bool isbp1first=false;
					bool isbp2first=false;
					if(it->BP1.IsLeft && genes1[i].Strand=='-')
						isbp1first=true;
					else if(!it->BP2.IsLeft && genes1[i].Strand=='+')
						isbp1first=true;
					if(it->BP2.IsLeft && genes2[j].Strand=='-')
						isbp2first=true;
					else if(!it->BP2.IsLeft && genes2[j].Strand=='+')
						isbp2first=true;
					assert(isbp1first^isbp2first);
					int bp1=(it->BP1.IsLeft)?it->BP1.StartPos:it->BP1.EndPos;
					int bp2=(it->BP2.IsLeft)?it->BP2.StartPos:it->BP2.EndPos;
					vector<Interval_t> tmpExons1, tmpExons2;
					for(vector<Interval_t>::iterator itexon=genes1[i].vExon.begin(); itexon!=genes1[i].vExon.end(); itexon++){
						if((it->BP1.IsLeft && itexon->EndPos>bp1) || (!it->BP1.IsLeft && itexon->StartPos<bp1))
							tmpExons1.push_back(*itexon);
						if(it->BP1.IsLeft && itexon->StartPos<bp1-thresh)
							tmpExons1.back().StartPos=bp1;
						if(!it->BP1.IsLeft && itexon->EndPos>bp1+thresh)
							tmpExons1.back().EndPos=bp1;
					}
					for(vector<Interval_t>::iterator itexon=genes2[j].vExon.begin(); itexon!=genes2[j].vExon.end(); itexon++){
						if((it->BP2.IsLeft && itexon->EndPos>bp2) || (!it->BP2.IsLeft && itexon->StartPos<bp2))
							tmpExons2.push_back(*itexon);
						if(it->BP2.IsLeft && itexon->StartPos<bp2-thresh)
							tmpExons2.back().StartPos=bp2;
						if(!it->BP2.IsLeft && itexon->EndPos>bp2+thresh)
							tmpExons2.back().EndPos=bp2;
					}
					if(it->BP1.IsLeft)
						reverse(tmpExons1.begin(), tmpExons1.end());
					if(it->BP2.IsLeft)
						reverse(tmpExons2.begin(), tmpExons2.end());
					vector<Interval_t> tmpExonAll; // fused transcript
					if(isbp1first){
						tmpExonAll.insert(tmpExonAll.end(), tmpExons1.begin(), tmpExons1.end());
						tmpExonAll.insert(tmpExonAll.end(), tmpExons2.begin(), tmpExons2.end());
					}
					else{
						tmpExonAll.insert(tmpExonAll.end(), tmpExons2.begin(), tmpExons2.end());
						tmpExonAll.insert(tmpExonAll.end(), tmpExons1.begin(), tmpExons1.end());
					}
					// write fused transcript to file
					fpfusion<<(tmpExons1[0].Chr)<<"\tfusedtrans\ttranscript\t"<<(tmpExons1[0].StartPos)<<"\t"<<(tmpExons1[0].EndPos)<<"\t.\t"<<(tmpExons1[0].Strand)<<"\t.\t";
					fpfusion<<"gene_id \"FG"<<"\"; transcript_id \"FGtrans\";\n";
					for(int i=0; i<tmpExonAll.size(); i++){
						const Interval_t exon=tmpExonAll[i];
						fpfusion<<(exon.Chr)<<"\tfusedtrans\texon\t"<<(exon.StartPos)<<"\t"<<(exon.EndPos)<<"\t.\t"<<(exon.Strand)<<"\t.\t";
						fpfusion<<"gene_id \"FG"<<"\"; transcript_id \"FGtrans\"; exon_number \""<<i<<"\"\n";
					}
				}
		}
		// if not fusion gene
		countseq++;
		int seg1=(it->BP1.IsLeft)?it->BP1.EndPos:it->BP1.StartPos; // other end of breakpoint segment
		int seg2=(it->BP2.IsLeft)?it->BP2.EndPos:it->BP2.StartPos;
		int bp1=(it->BP1.IsLeft)?it->BP1.StartPos:it->BP1.EndPos;
		int bp2=(it->BP2.IsLeft)?it->BP2.StartPos:it->BP2.EndPos;
		bool potential1=false, potential2=false;
		if(genes1.size()!=0)
			for(int i=0; i<genes1.size(); i++){
				if(it->BP1.IsLeft && genes1[i].TxEnd>seg1)
					seg1=genes1[i].TxEnd+50000;
				else if(!it->BP1.IsLeft && genes1[i].TxStart<seg1)
					seg1=genes1[i].TxStart-50000;
				for(vector<Interval_t>::iterator itexon=genes1[i].vExon.begin(); itexon!=genes1[i].vExon.end(); itexon++){
					if(it->BP1.IsLeft && genes1[i].Strand=='-' && abs(bp1-itexon->StartPos)<thresh){
						bp1=itexon->StartPos;
						potential1=true;
					}
					else if(!it->BP2.IsLeft && genes1[i].Strand=='+' && abs(bp1-itexon->EndPos)<thresh){
						bp1=itexon->EndPos;
						potential1=true;
					}
				}
			}
		if(genes2.size()!=0)
			for(int j=0; j<genes2.size(); j++){
				if(it->BP2.IsLeft && genes2[j].TxEnd>seg2)
					seg2=genes2[j].TxEnd+50000;
				else if(!it->BP2.IsLeft && genes2[j].TxStart<seg2)
					seg2=genes2[j].TxStart-50000;
				for(vector<Interval_t>::iterator itexon=genes2[j].vExon.begin(); itexon!=genes2[j].vExon.end(); itexon++){
					if(it->BP2.IsLeft && genes2[j].Strand=='-' && abs(bp2-itexon->StartPos)<thresh){
						bp2=itexon->StartPos;
						potential2=true;
					}
					else if(!it->BP2.IsLeft && genes2[j].Strand=='+' && abs(bp2-itexon->EndPos)<thresh){
						bp2=itexon->EndPos;
						potential2=true;
					}
				}
			}
		int start1=min(bp1, seg1), end1=max(bp1, seg1);
		int start2=min(bp2, seg2), end2=max(bp2, seg2);
		bool reversecomp1=(it->BP1.IsLeft), reversecomp2=(!it->BP2.IsLeft);
		string potential="none";
		if(potential1 && potential2)
			potential="either";
		else if(potential1)
			potential="bp1";
		else if(potential2)
			potential="bp2";
		// write sequence
		string tmp1=Genome[it->BP1.Chr].substr(start1, end1-start1);
		string tmp2=Genome[it->BP2.Chr].substr(start2, end2-start2);
		if(reversecomp1)
			reverse(tmp1.begin(), tmp1.end());
		if(reversecomp2)
			reverse(tmp2.begin(), tmp2.end());
		tmp1+=tmp2;
		fpseq<<">"<<countseq<<"\t"<<it->BP1.Chr<<" "<<start1<<" "<<end1<<"\t"<<it->BP2.Chr<<" "<<start2<<" "<<end2<<"\tpotential="<<potential<<endl;
		int nt=0;
		while(nt<tmp1.size()){
			int length=min(80, (int)tmp1.size()-nt);
			fpseq<<tmp1.substr(nt, nt+80)<<endl;
			nt+=length;
		}
		// write annotation
		for(int i=0; i<genes1.size(); i++){
			char Strand=genes1[i].Strand;
			if(reversecomp1)
				if(Strand=='-')
					Strand='+';
				else
					Strand='-';
			fptmpannot<<countseq<<"\tconverted\ttranscript\t"<<((start1<genes1[i].TxStart)?(genes1[i].TxStart-start1):0)<<"\t"<<((end1>genes1[i].TxEnd)?(genes1[i].TxEnd-start1):(end1-start1))<<"\t.\t";
			fptmpannot<<Strand<<"\t.\t"<<"gene_name \""<<genes1[i].GeneName<<"\"; transcript_id \""<<genes1[i].TransID<<"\";\n";
			int exonnumber=0;
			for(vector<Interval_t>::iterator itexon=genes1[i].vExon.begin(); itexon!=genes1[i].vExon.end(); it++){
				if(itexon->EndPos<start1 || itexon->StartPos>end1)
					continue;
				exonnumber++;
				fptmpannot<<countseq<<"\tconverted\texon\t"<<((start1<itexon->StartPos)?(itexon->StartPos-start1):0)<<"\t"<<((end1>itexon->EndPos)?(itexon->EndPos-start1):(end1-start1))<<"\t.\t";
				fptmpannot<<Strand<<"\t.\t"<<"gene_name \""<<genes1[i].GeneName<<"\"; transcript_id \""<<genes1[i].TransID<<"\"; exon_number \""<<exonnumber<<"\";\n";
			}
		}
		for(int j=0; j<genes2.size(); j++){
			int offset=end1-start1;
			char Strand=genes2[j].Strand;
			if(reversecomp2)
				if(Strand=='-')
					Strand='+';
				else
					Strand='-';
			fptmpannot<<countseq<<"\tconverted\ttranscript\t"<<((start2<genes2[j].TxStart)?(genes2[j].TxStart-start2+offset):(offset))<<"\t"<<((end2>genes2[j].TxEnd)?(genes2[j].TxEnd-start2+offset):(end2-start2+offset))<<"\t.\t";
			fptmpannot<<Strand<<"\t.\t"<<"gene_name \""<<genes2[j].GeneName<<"\"; transcript_id \""<<genes2[j].TransID<<"\";\n";
			int exonnumber=0;
			for(vector<Interval_t>::iterator itexon=genes2[j].vExon.begin(); itexon!=genes2[j].vExon.end(); itexon++){
				if(itexon->EndPos<start2 || itexon->StartPos>end2)
					continue;
				exonnumber++;
				fptmpannot<<countseq<<"\tconverted\texon\t"<<((start2<itexon->StartPos)?(itexon->StartPos-start2+offset):offset)<<"\t"<<((end2>itexon->EndPos)?(itexon->EndPos-start2+offset):(end2-start2+offset))<<"\t.\t";
				fptmpannot<<Strand<<"\t.\t"<<"gene_name \""<<genes2[j].GeneName<<"\"; transcript_id \""<<genes2[j].TransID<<"\"; exon_number \""<<exonnumber<<"\";\n";
			}
		}
	}
	fpfusion.close();
	fpseq.close();
	fptmpannot.close();
};


int main(int argc, char* argv[]){
	string genomefile;
	string SVpredfile;
	string annotationfile;
	string OutPrefix;

	if(argc==1)
		printf("SVTranscript <genomefile> <SVpredfile> <annotationfile> <OutPrefix>\n");
	else{
		genomefile=(string)argv[1];
		SVpredfile=(string)argv[2];
		annotationfile=(string)argv[3];
		OutPrefix=(string)argv[4];

		map<string,string> Genome;
		ReadGenome(genomefile, Genome);

		vector<SV_t> SVs;
		ReadBedpeSV(SVpredfile, SVs);

		vector<Transcript_t> vTrans;
		ReadGTF(annotationfile, vTrans);

		InitialJunction(SVs, vTrans, Genome, OutPrefix);
	}
}