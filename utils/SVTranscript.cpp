// module 1 (initial): generate fused genome junction for non-fusion-gene TSV
// module 2 (initial): generate temporary gtf file for STAR alignment
// module 3 (final): generate sequence for junction transcript

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
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
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "noshell.hpp"

using namespace std;

static std::map<char,char> Nucleotide={{'A','T'},{'C','G'},{'G','C'},{'T','A'},{'R','Y'},{'Y','R'},{'S','W'},{'W','S'},{'K','M'},{'M','K'},{'B','V'},{'V','B'},{'D','H'},{'H','D'}, {'N','N'}, {'.','.'},{'-','-'}};
void ReverseComplement(string::iterator itbegin, string::iterator itend){
	for(string::iterator it=itbegin; it!=itend; it++)
		*it=Nucleotide[toupper(*it)];
	std::reverse(itbegin, itend);
};

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

vector<Transcript_t> GeneBinarySearch(vector<Transcript_t>::iterator firstiterator, vector<Transcript_t>::iterator lastiterator, BP_t bp){
	vector<Transcript_t> genes;
	vector<Transcript_t>::iterator itbegin=firstiterator;
	vector<Transcript_t>::iterator itend=lastiterator;
	int bppos=(bp.IsLeft)?bp.StartPos:bp.EndPos;
	vector<Transcript_t>::iterator itmid;
	while(distance(itbegin, itend)>1){
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
	for(; distance(itbegin, itmid)<20 && distance(firstiterator, itbegin)>=0; itbegin--){
		if(itbegin->Chr==bp.Chr && itbegin->TxStart<=bppos && itbegin->TxEnd>=bppos)
			genes.push_back(*itbegin);
	}
	itend++;
	for(; distance(itmid, itend)<20 && itend!=lastiterator; itend++){
		if(itend->Chr==bp.Chr && itend->TxStart<=bppos && itend->TxEnd>=bppos)
			genes.push_back(*itend);
	}
	return genes;
};

bool IsIsoform(SV_t sv1, SV_t sv2, vector<Transcript_t>& genes1, vector<Transcript_t>& genes2, vector<Transcript_t>& tmpgenes1, vector<Transcript_t>& tmpgenes2){
	int thresh=5000;
	if(sv1.BP1.Chr==sv2.BP1.Chr && sv1.BP2.Chr==sv2.BP2.Chr && sv1.BP1.IsLeft==sv2.BP1.IsLeft && sv1.BP2.IsLeft==sv2.BP2.IsLeft){
		bool bp1near=false;
		bool bp2near=false;
		vector<string> bp1overlapgene;
		vector<string> bp2overlapgene;
		// check whether bp1 of svs has overlapping gene
		if(genes1.size()!=0 && tmpgenes1.size()!=0){
			vector<string> genename;
			vector<string> tmpgenename;
			for(vector<Transcript_t>::iterator ittrans=genes1.begin(); ittrans!=genes1.end(); ittrans++)
				genename.push_back(ittrans->TransID);
			for(vector<Transcript_t>::iterator ittrans=tmpgenes1.begin(); ittrans!=tmpgenes1.end(); ittrans++)
				tmpgenename.push_back(ittrans->TransID);
			sort(genename.begin(), genename.end());
			sort(tmpgenename.begin(), tmpgenename.end());
			set_intersection(genename.begin(), genename.end(), tmpgenename.begin(), tmpgenename.end(), back_inserter(bp1overlapgene));
		}
		// check whether bp2 of svs has overlapping gene
		if(genes2.size()!=0 && tmpgenes2.size()!=0){
			vector<string> genename;
			vector<string> tmpgenename;
			for(vector<Transcript_t>::iterator ittrans=genes2.begin(); ittrans!=genes2.end(); ittrans++)
				genename.push_back(ittrans->TransID);
			for(vector<Transcript_t>::iterator ittrans=tmpgenes2.begin(); ittrans!=tmpgenes2.end(); ittrans++)
				tmpgenename.push_back(ittrans->TransID);
			sort(genename.begin(), genename.end());
			sort(tmpgenename.begin(), tmpgenename.end());
			set_intersection(genename.begin(), genename.end(), tmpgenename.begin(), tmpgenename.end(), back_inserter(bp2overlapgene));
		}
		// whether TSV isoform
		if(!((sv1.BP1.StartPos<sv2.BP1.StartPos)==(sv1.BP1.EndPos<sv2.BP1.StartPos) && (sv1.BP1.StartPos<sv2.BP1.StartPos)==(sv1.BP1.StartPos<sv2.BP1.EndPos)))
			bp1near=true;
		else if(bp1overlapgene.size()!=0)
			bp1near=true;
		else if(genes1.size()==0 && tmpgenes1.size()==0 && (abs(sv1.BP1.StartPos-sv2.BP1.StartPos)<thresh || abs(sv1.BP1.EndPos-sv2.BP1.EndPos)<thresh))
			bp1near=true;
		if(!((sv1.BP2.StartPos<sv2.BP2.StartPos)==(sv1.BP2.EndPos<sv2.BP2.StartPos) && (sv1.BP2.StartPos<sv2.BP2.StartPos)==(sv1.BP2.StartPos<sv2.BP2.EndPos)))
			bp2near=true;
		else if(bp2overlapgene.size()!=0)
			bp2near=true;
		else if(genes2.size()==0 && tmpgenes2.size()==0 && (abs(sv1.BP2.StartPos-sv2.BP2.StartPos)<thresh || abs(sv1.BP2.EndPos-sv2.BP2.EndPos)<thresh))
			bp2near=true;
		// if isoform, only keep shared genes, discard bp specific genes
		if(bp1near && bp2near){
			vector<Transcript_t> tmpgenes;
			for(vector<Transcript_t>::iterator ittrans=genes1.begin(); ittrans!=genes1.end(); ittrans++)
				if(binary_search(bp1overlapgene.begin(), bp1overlapgene.end(), ittrans->TransID))
					tmpgenes.push_back(*ittrans);
			genes1=tmpgenes;
			tmpgenes.clear();
			for(vector<Transcript_t>::iterator ittrans=genes2.begin(); ittrans!=genes2.end(); ittrans++)
				if(binary_search(bp2overlapgene.begin(), bp2overlapgene.end(), ittrans->TransID))
					tmpgenes.push_back(*ittrans);
			genes2=tmpgenes;
		}
		return (bp1near && bp2near);
	}
	return false;
};

void WriteFusionGene(vector<SV_t>::iterator firstiterator, vector<SV_t>::iterator lastiterator, vector<Transcript_t>& genes1, vector<Transcript_t>& genes2, ofstream& fpfusion, string fuseID){
	int thresh=5;
	for(vector<SV_t>::iterator ittmp=firstiterator; ittmp!=lastiterator; ittmp++){
		for(int i=0; i<genes1.size(); i++)
			for(int j=0; j<genes2.size(); j++){
				if((ittmp->BP1.IsLeft==ittmp->BP2.IsLeft)==(genes1[i].Strand==genes2[j].Strand))
					continue;
				bool isbp1first=false;
				bool isbp2first=false;
				if(ittmp->BP1.IsLeft && genes1[i].Strand=='-')
					isbp1first=true;
				else if(!ittmp->BP1.IsLeft && genes1[i].Strand=='+')
					isbp1first=true;
				if(ittmp->BP2.IsLeft && genes2[j].Strand=='-')
					isbp2first=true;
				else if(!ittmp->BP2.IsLeft && genes2[j].Strand=='+')
					isbp2first=true;
				assert(isbp1first^isbp2first);
				int bp1=(ittmp->BP1.IsLeft)?ittmp->BP1.StartPos:ittmp->BP1.EndPos;
				int bp2=(ittmp->BP2.IsLeft)?ittmp->BP2.StartPos:ittmp->BP2.EndPos;
				vector<Interval_t> tmpExons1, tmpExons2;
				for(vector<Interval_t>::iterator itexon=genes1[i].vExon.begin(); itexon!=genes1[i].vExon.end(); itexon++){
					if((ittmp->BP1.IsLeft && itexon->EndPos>bp1) || (!ittmp->BP1.IsLeft && itexon->StartPos<bp1)){
						tmpExons1.push_back(*itexon);
						if(ittmp->BP1.IsLeft && itexon->StartPos<bp1-thresh)
							tmpExons1.back().StartPos=bp1;
						if(!ittmp->BP1.IsLeft && itexon->EndPos>bp1+thresh)
							tmpExons1.back().EndPos=bp1;
					}
				}
				for(vector<Interval_t>::iterator itexon=genes2[j].vExon.begin(); itexon!=genes2[j].vExon.end(); itexon++){
					if((ittmp->BP2.IsLeft && itexon->EndPos>bp2) || (!ittmp->BP2.IsLeft && itexon->StartPos<bp2)){
						tmpExons2.push_back(*itexon);
						if(ittmp->BP2.IsLeft && itexon->StartPos<bp2-thresh)
							tmpExons2.back().StartPos=bp2;
						if(!ittmp->BP2.IsLeft && itexon->EndPos>bp2+thresh)
							tmpExons2.back().EndPos=bp2;
					}
				}
				if(ittmp->BP1.IsLeft)
					reverse(tmpExons1.begin(), tmpExons1.end());
				if(ittmp->BP2.IsLeft)
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
				fpfusion<<"gene_id \"FG"<<fuseID<<"\"; transcript_id \"FGtrans"<<fuseID<<"\";\n";
				for(int i=0; i<tmpExonAll.size(); i++){
					const Interval_t exon=tmpExonAll[i];
					fpfusion<<(exon.Chr)<<"\tfusedtrans\texon\t"<<(exon.StartPos)<<"\t"<<(exon.EndPos)<<"\t.\t"<<(exon.Strand)<<"\t.\t";
					fpfusion<<"gene_id \"FG"<<fuseID<<"\"; transcript_id \"FGtrans"<<fuseID<<"\"; exon_number \""<<i<<"\"\n";
				}
			}
	}
};

SV_t NonFusionGeneJunction(vector<SV_t>::iterator firstiterator, vector<SV_t>::iterator lastiterator, vector<Transcript_t>& genes1, vector<Transcript_t>& genes2, bool& potential1, bool& potential2, map<string,string>& Genome){
	int thresh=5;
	int seg1=(firstiterator->BP1.IsLeft)?(firstiterator->BP1.EndPos):(firstiterator->BP1.StartPos); // other end of breakpoint segment
	int seg2=(firstiterator->BP2.IsLeft)?(firstiterator->BP2.EndPos):(firstiterator->BP2.StartPos); 
	int bp1=(firstiterator->BP1.IsLeft)?(firstiterator->BP1.StartPos):(firstiterator->BP1.EndPos); // breakpoint end
	int bp2=(firstiterator->BP2.IsLeft)?(firstiterator->BP2.StartPos):(firstiterator->BP2.EndPos);
	for(vector<SV_t>::iterator ittmp=firstiterator; ittmp!=lastiterator; ittmp++){
		if(ittmp->BP1.IsLeft){
			bp1=(bp1<ittmp->BP1.StartPos)?bp1:ittmp->BP1.StartPos;
			seg1=(seg1>ittmp->BP1.EndPos)?seg1:ittmp->BP1.EndPos;
		}
		else{
			bp1=(bp1>ittmp->BP1.EndPos)?bp1:ittmp->BP1.EndPos;
			seg1=(seg1<ittmp->BP1.StartPos)?seg1:ittmp->BP1.StartPos;
		}
		if(ittmp->BP2.IsLeft){
			bp2=(bp2<ittmp->BP2.StartPos)?bp2:ittmp->BP2.StartPos;
			seg2=(seg2>ittmp->BP2.EndPos)?seg2:ittmp->BP2.EndPos;
		}
		else{
			bp2=(bp2>ittmp->BP2.EndPos)?bp2:ittmp->BP2.EndPos;
			seg2=(seg2<ittmp->BP2.StartPos)?seg2:ittmp->BP2.StartPos;
		}
	}
	int rec1=seg1, rec2=seg2; // keep record of initial right or left most segment boundary.
	potential1=false;
	potential2=false;
	// move seg1 and seg2 to the end of overlapping gene
	if(genes1.size()!=0)
		for(int i=0; i<genes1.size(); i++){
			if(firstiterator->BP1.IsLeft && genes1[i].TxEnd>seg1)
				seg1=genes1[i].TxEnd;
			else if(!firstiterator->BP1.IsLeft && genes1[i].TxStart<seg1)
				seg1=genes1[i].TxStart;
			for(vector<Interval_t>::iterator itexon=genes1[i].vExon.begin(); itexon!=genes1[i].vExon.end(); itexon++){
				for(vector<SV_t>::iterator ittmp=firstiterator; ittmp!=lastiterator; ittmp++){
					int tmpbp=(ittmp->BP1.IsLeft)?(ittmp->BP1.StartPos):(ittmp->BP1.EndPos);
					if(ittmp->BP1.IsLeft && genes1[i].Strand=='-' && abs(tmpbp-itexon->StartPos)<thresh){
						ittmp->BP1.StartPos=itexon->StartPos;
						potential1=true;
						if(abs(bp1-itexon->StartPos)<thresh)
							bp1=itexon->StartPos;
					}
					else if(!ittmp->BP1.IsLeft && genes1[i].Strand=='+' && abs(tmpbp-itexon->EndPos)<thresh){
						ittmp->BP1.EndPos=itexon->EndPos;
						potential1=true;
						if(abs(bp1-itexon->EndPos)<thresh)
							bp1=itexon->EndPos;
					}
				}
			}
		}
	if(genes2.size()!=0)
		for(int j=0; j<genes2.size(); j++){
			if(firstiterator->BP2.IsLeft && genes2[j].TxEnd>seg2)
				seg2=genes2[j].TxEnd;
			else if(!firstiterator->BP2.IsLeft && genes2[j].TxStart<seg2)
				seg2=genes2[j].TxStart;
			for(vector<Interval_t>::iterator itexon=genes2[j].vExon.begin(); itexon!=genes2[j].vExon.end(); itexon++){
				for(vector<SV_t>::iterator ittmp=firstiterator; ittmp!=lastiterator; ittmp++){
					int tmpbp=(ittmp->BP2.IsLeft)?(ittmp->BP2.StartPos):(ittmp->BP2.EndPos);
					if(ittmp->BP2.IsLeft && genes2[j].Strand=='-' && abs(tmpbp-itexon->StartPos)<thresh){
						ittmp->BP2.StartPos=itexon->StartPos;
						potential2=true;
						if(abs(bp2-itexon->StartPos)<thresh)
							bp2=itexon->StartPos;
					}
					else if(!ittmp->BP2.IsLeft && genes2[j].Strand=='+' && abs(tmpbp-itexon->EndPos)<thresh){
						ittmp->BP2.EndPos=itexon->EndPos;
						potential2=true;
						if(abs(bp2-itexon->EndPos)<thresh)
							bp2=itexon->EndPos;
					}
				}
			}
		}
	// add more free space, in case of antisense fusion where codon is not preserved
	if(firstiterator->BP1.IsLeft && rec1+50000>seg1 && rec1+50000<Genome[firstiterator->BP1.Chr].size())
		seg1=rec1+50000;
	if(!firstiterator->BP1.IsLeft && rec1-50000<seg1 && rec1-50000>0)
		seg1=rec1-50000;
	if(firstiterator->BP2.IsLeft && rec2+50000>seg2 && rec2+50000<Genome[firstiterator->BP2.Chr].size())
		seg2=rec2+50000;
	if(!firstiterator->BP2.IsLeft && rec2-50000<seg2 && rec2-50000>0)
		seg2=rec2-50000;
	int start1=min(bp1, seg1), end1=max(bp1, seg1);
	int start2=min(bp2, seg2), end2=max(bp2, seg2);
	if(firstiterator->BP1.Chr==firstiterator->BP2.Chr && start2<end1){
		if(firstiterator->BP1.IsLeft && firstiterator->BP2.IsLeft)
			end1=(rec1+bp2)/2;
		else if(firstiterator->BP1.IsLeft && !firstiterator->BP2.IsLeft)
			end1=(rec1+rec2)/2;
		else if(!firstiterator->BP1.IsLeft && !firstiterator->BP2.IsLeft)
			start2=(bp1+rec2)/2;
	}
	bool reversecomp1=(firstiterator->BP1.IsLeft), reversecomp2=(!firstiterator->BP2.IsLeft);
	BP_t newbp1(firstiterator->BP1.Chr, start1, end1, reversecomp1); // IsLeft is used as IsReverseComplement
	BP_t newbp2(firstiterator->BP2.Chr, start2, end2, reversecomp2); // IsLeft is used as IsReverseComplement
	SV_t sv(newbp1, newbp2);
	return sv;
};

void WriteNonFusionGene(vector<SV_t>::iterator firstiterator, vector<SV_t>::iterator lastiterator, map<string, string>& Genome, vector<Transcript_t>& genes1, vector<Transcript_t>& genes2, SV_t& sv, bool potential1, bool potential2, ofstream& fpseq, ofstream& fptmpannot, string fuseID){
	int thresh=5;
	string potential="none";
	if(potential1 && potential2)
		potential="either";
	else if(potential1)
		potential="bp1";
	else if(potential2)
		potential="bp2";
	int start1=sv.BP1.StartPos, end1=sv.BP1.EndPos;
	int start2=sv.BP2.StartPos, end2=sv.BP2.EndPos;
	bool reversecomp1=sv.BP1.IsLeft, reversecomp2=sv.BP2.IsLeft;
	// write sequence
	string tmp1=Genome[firstiterator->BP1.Chr].substr(start1, end1-start1);
	string tmp2=Genome[firstiterator->BP2.Chr].substr(start2, end2-start2);
	if(reversecomp1)
		ReverseComplement(tmp1.begin(), tmp1.end());
	if(reversecomp2)
		ReverseComplement(tmp2.begin(), tmp2.end());
	tmp1+=tmp2;
	fpseq<<">"<<fuseID<<"\t"<<firstiterator->BP1.Chr<<" "<<start1<<" "<<end1<<"\t"<<firstiterator->BP2.Chr<<" "<<start2<<" "<<end2<<"\tpotential="<<potential<<endl;
	int nt=0;
	while(nt<tmp1.size()){
		int length=min(80, (int)tmp1.size()-nt);
		fpseq<<tmp1.substr(nt, 80)<<endl;
		nt+=length;
	}
	// write annotation for known genes within junction region
	for(int i=0; i<genes1.size(); i++){
		char Strand=genes1[i].Strand;
		if(reversecomp1)
			if(Strand=='-')
				Strand='+';
			else
				Strand='-';
		fptmpannot<<fuseID<<"\tconverted\ttranscript\t"<<((start1<genes1[i].TxStart)?(genes1[i].TxStart-start1):0)<<"\t"<<((end1>genes1[i].TxEnd)?(genes1[i].TxEnd-start1):(end1-start1))<<"\t.\t";
		fptmpannot<<Strand<<"\t.\t"<<"gene_name \""<<genes1[i].GeneName<<"\"; transcript_id \""<<genes1[i].TransID<<"\";\n";
		int exonnumber=0;
		for(vector<Interval_t>::iterator itexon=genes1[i].vExon.begin(); itexon!=genes1[i].vExon.end(); itexon++){
			if(itexon->EndPos<start1 || itexon->StartPos>end1)
				continue;
			exonnumber++;
			fptmpannot<<fuseID<<"\tconverted\texon\t"<<((start1<itexon->StartPos)?(itexon->StartPos-start1):0)<<"\t"<<((end1>itexon->EndPos)?(itexon->EndPos-start1):(end1-start1))<<"\t.\t";
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
		fptmpannot<<fuseID<<"\tconverted\ttranscript\t"<<((start2<genes2[j].TxStart)?(genes2[j].TxStart-start2+offset):(offset))<<"\t"<<((end2>genes2[j].TxEnd)?(genes2[j].TxEnd-start2+offset):(end2-start2+offset))<<"\t.\t";
		fptmpannot<<Strand<<"\t.\t"<<"gene_name \""<<genes2[j].GeneName<<"\"; transcript_id \""<<genes2[j].TransID<<"\";\n";
		int exonnumber=0;
		for(vector<Interval_t>::iterator itexon=genes2[j].vExon.begin(); itexon!=genes2[j].vExon.end(); itexon++){
			if(itexon->EndPos<start2 || itexon->StartPos>end2)
				continue;
			exonnumber++;
			fptmpannot<<fuseID<<"\tconverted\texon\t"<<((start2<itexon->StartPos)?(itexon->StartPos-start2+offset):offset)<<"\t"<<((end2>itexon->EndPos)?(itexon->EndPos-start2+offset):(end2-start2+offset))<<"\t.\t";
			fptmpannot<<Strand<<"\t.\t"<<"gene_name \""<<genes2[j].GeneName<<"\"; transcript_id \""<<genes2[j].TransID<<"\"; exon_number \""<<exonnumber<<"\";\n";
		}
	}
	// write isoform TSV junction that are not directly adjacent in sequence
	for(vector<SV_t>::iterator ittmp=firstiterator; ittmp!=lastiterator; ittmp++){
		int bp1=(ittmp->BP1.IsLeft)?ittmp->BP1.StartPos:ittmp->BP1.EndPos;
		int bp2=(ittmp->BP2.IsLeft)?ittmp->BP2.StartPos:ittmp->BP2.EndPos;
		int junction1=(sv.BP1.IsLeft)?sv.BP1.StartPos:sv.BP1.EndPos;
		int junction2=(!sv.BP2.IsLeft)?sv.BP2.StartPos:sv.BP2.EndPos;
		if(abs(bp1-junction1)<=thresh && abs(bp2-junction2)<=thresh)
			continue;
		if(sv.BP1.IsLeft)
			bp1=end1-bp1;
		else
			bp1-=start1;
		if(sv.BP2.IsLeft)
			bp2=end1-start1+end2-bp2;
		else
			bp2=end1-start1+bp2-start2;
		string newtransID=fuseID+"00"+to_string(distance(firstiterator, ittmp));
		fptmpannot<<fuseID<<"\tfusion\ttranscript\t"<<max(0, bp1-100)<<"\t"<<min(bp2+100, end1-start1+end2-start2)<<"\t";
		fptmpannot<<"+"<<"\t.\t"<<"gene_name \"neofusion"<<newtransID<<"\"; transcript_id \"neotrans"<<newtransID<<"\"; exon_number \"1\";\n";
		// exon1
		fptmpannot<<fuseID<<"\tfusion\texon\t"<<max(0, bp1-100)<<"\t"<<bp1<<"\t";
		fptmpannot<<"+"<<"\t.\t"<<"gene_name \"neofusion"<<newtransID<<"\"; transcript_id \"neotrans"<<newtransID<<"\" ;\n";
		// exon2
		fptmpannot<<fuseID<<"\tfusion\texon\t"<<bp2<<"\t"<<min(bp2+100, end1-start1+end2-start2)<<"\t";
		fptmpannot<<"+"<<"\t.\t"<<"gene_name \"neofusion"<<newtransID<<"\"; transcript_id \"neotrans"<<newtransID<<"\"; exon_number \"2\";\n";
	}
};

void InitialJunction(vector<SV_t>& SVs, vector<Transcript_t>& vTrans, map<string, string>& Genome, string OutPrefix, vector<SV_t>& JunctionRegion, vector<bool>& Potent1, vector<bool>& Potent2){
	int countseq=0;
	int countfg=0;
	JunctionRegion.clear();
	Potent1.clear();
	Potent2.clear();
	ofstream fpfusion(OutPrefix+"/fusiongene.gtf", ios::out);
	ofstream fpseq(OutPrefix+"/juncseq.fa", ios::out);
	ofstream fptmpannot(OutPrefix+"/juncannot.gtf", ios::out);
	vector<SV_t>::iterator itbegin=SVs.begin();
	vector<Transcript_t> genes1=GeneBinarySearch(vTrans.begin(), vTrans.end(), itbegin->BP1);
	vector<Transcript_t> genes2=GeneBinarySearch(vTrans.begin(), vTrans.end(), itbegin->BP2);
	for(vector<SV_t>::iterator it=SVs.begin()+1; ; it++){
		vector<Transcript_t> tmpgenes1, tmpgenes2;
		if(it!=SVs.end()){
			tmpgenes1=GeneBinarySearch(vTrans.begin(), vTrans.end(), it->BP1);
			tmpgenes2=GeneBinarySearch(vTrans.begin(), vTrans.end(), it->BP2);
		}
		if(it!=SVs.end() && IsIsoform(*itbegin, *it, genes1, genes2, tmpgenes1, tmpgenes2))
			continue;
		bool is_fusion=false;
		if(genes1.size()!=0 && genes2.size()!=0){
			for(int i=0; i<genes1.size(); i++)
				for(int j=0; j<genes2.size(); j++){
					if((itbegin->BP1.IsLeft==itbegin->BP2.IsLeft)!=(genes1[i].Strand==genes2[j].Strand))
						is_fusion=true;
				}
		}
		// if fusion gene, output new transcript directly
		if(is_fusion){
			countfg++;
			ostringstream fuseID;
			fuseID<<setfill('0')<<setw(4)<<countfg;
			WriteFusionGene(itbegin, it, genes1, genes2, fpfusion, fuseID.str());
		}
		// if not fusion gene
		else{
			countseq++;
			bool potential1, potential2;
			SV_t sv=NonFusionGeneJunction(itbegin, it, genes1, genes2, potential1, potential2, Genome);
			JunctionRegion.push_back(sv);
			Potent1.push_back(potential1);
			Potent2.push_back(potential2);
			WriteNonFusionGene(itbegin, it, Genome, genes1, genes2, sv, potential1, potential2, fpseq, fptmpannot, to_string(countseq));
		}
		genes1=tmpgenes1;
		genes2=tmpgenes2;
		itbegin=it;
		if(it==SVs.end())
			break;
	}
	fpfusion.close();
	fpseq.close();
	fptmpannot.close();
};

void WriteFastq(string FqPrefix, string OutPrefix, vector<string>& Rnames){
	sort(Rnames.begin(), Rnames.end());
	vector<string>::iterator stringit=unique(Rnames.begin(), Rnames.end());
	Rnames.resize(distance(Rnames.begin(), stringit));

	ifstream input(FqPrefix+"_1.fastq");
	ofstream output(OutPrefix+"_1.fastq", ios::out);
	string line;
	int count=0;
	bool flag=false;
	while(getline(input, line)){
		if(count%4==0){
			vector<string> strs;
			boost::split(strs, line, boost::is_any_of(" "));
			string qname=strs[0].substr(1);
			if(qname.substr(qname.size()-2, 2)=="/1" || qname.substr(qname.size()-2, 2)=="/2")
				qname=qname.substr(0, qname.size()-2);
			if(binary_search(Rnames.begin(), Rnames.end(), qname))
				flag=true;
			else
				flag=false;
		}
		if(flag)
			output<<line<<endl;
		count++;
	}
	input.close();
	output.close();

	input.open(FqPrefix+"_2.fastq");
	output.open(OutPrefix+"_2.fastq", ios::out);
	count=0;
	flag=false;
	while(getline(input, line)){
		if(count%4==0){
			vector<string> strs;
			boost::split(strs, line, boost::is_any_of(" "));
			string qname=strs[0].substr(1);
			if(qname.substr(qname.size()-2, 2)=="/1" || qname.substr(qname.size()-2, 2)=="/2")
				qname=qname.substr(0, qname.size()-2);
			if(binary_search(Rnames.begin(), Rnames.end(), qname))
				flag=true;
			else
				flag=false;
		}
		if(flag)
			output<<line<<endl;
		count++;
	}
	input.close();
	output.close();
};

int RegionalReadWritter(string BamPrefix, string ConcordBam, string ChimericBam, string FqPrefix, string OutPrefix, vector<SV_t>& JunctionRegion){
	vector<string> Rnames;
	samFile * bamfile=sam_open((BamPrefix+"/"+ConcordBam).c_str(), "r");
	hts_idx_t *idx=NULL;
	bam_hdr_t *header=NULL;
	bam1_t *b=NULL;
	hts_itr_t *iter=NULL;

	// Read concordant bam file
	if(bamfile==NULL){
		cout<<"cannot open bam file\n";
		return -1;
	}
	if((header=sam_hdr_read(bamfile))==0){
		cout<<"cannot open header\n";
		return -1;
	}
	idx=sam_index_load(bamfile, (BamPrefix+"/"+ConcordBam).c_str());
	if(idx==NULL){
		cout<<"cannot load index\n";
		return -1;
	}
	for(vector<SV_t>::iterator it=JunctionRegion.begin(); it!=JunctionRegion.end(); it++){
		// for bp1 region
		int tid=-1;
		for(int i=0; i<header->n_targets; i++)
			if((string)header->target_name[i]==it->BP1.Chr)
				tid=i;
		assert(tid>=0);
		iter=sam_itr_queryi(idx, tid, it->BP1.StartPos, it->BP1.EndPos);
		b=bam_init1();
		while(sam_itr_next(bamfile, iter, b)>=0){
			// check if read and mate within region, and add read name to list.
			if(b->core.pos>=it->BP1.StartPos && b->core.pos<=it->BP1.EndPos && b->core.mtid==tid && b->core.mpos>=it->BP1.StartPos && b->core.mpos<=it->BP1.EndPos)
				Rnames.push_back((string)bam_get_qname(b));
		}
		// for bp2 region
		tid=-1;
		for(int i=0; i<header->n_targets; i++)
			if((string)header->target_name[i]==it->BP2.Chr)
				tid=i;
		assert(tid>=0);
		iter=sam_itr_queryi(idx, tid, it->BP2.StartPos, it->BP2.EndPos);
		b=bam_init1();
		while(sam_itr_next(bamfile, iter, b)>=0){
			// check if read and mate within region, and add read name to list.
			if(b->core.pos>=it->BP2.StartPos && b->core.pos<=it->BP2.EndPos && b->core.mtid==tid && b->core.mpos>=it->BP2.StartPos && b->core.mpos<=it->BP2.EndPos)
				Rnames.push_back((string)bam_get_qname(b));
		}
	}
	hts_itr_destroy(iter);
	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(bamfile);

	// Read chimeric bam file
	bamfile=sam_open((BamPrefix+"/"+ChimericBam).c_str(), "r");
	if(bamfile==NULL){
		cout<<"cannot open bam file\n";
		return -1;
	}
	if((header=sam_hdr_read(bamfile))==0){
		cout<<"cannot open header\n";
		return -1;
	}
	b=bam_init1();
	while(sam_read1(bamfile, header, b)>=0){
		// check if read and mate within region, and add read name to list.
		bool flag1=false, flag2=false;
		for(vector<SV_t>::iterator it=JunctionRegion.begin(); it!=JunctionRegion.end(); it++){
			if((string)header->target_name[b->core.tid]==it->BP1.Chr && b->core.pos>=it->BP1.StartPos && b->core.pos<=it->BP1.EndPos)
				flag1=true;
			else if((string)header->target_name[b->core.tid]==it->BP2.Chr && b->core.pos>=it->BP2.StartPos && b->core.pos<=it->BP2.EndPos)
				flag1=true;
			if((string)header->target_name[b->core.mtid]==it->BP1.Chr && b->core.mpos>=it->BP1.StartPos && b->core.mpos<=it->BP1.EndPos)
				flag2=true;
			else if((string)header->target_name[b->core.mtid]==it->BP2.Chr && b->core.mpos>=it->BP2.StartPos && b->core.mpos<=it->BP2.EndPos)
				flag2=true;
		}
		if(flag1 && flag2)
			Rnames.push_back((string)bam_get_qname(b));
	}
	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(bamfile);

	// Write fastq file from query name vector
	WriteFastq(FqPrefix, OutPrefix, Rnames);
	return 0;
};

int Aligning(string STARpath, string STARindex, int starthreads, string ReadPrefix, string OutPrefix){
	noshell::Exit e1=noshell::C("mkdir", "-p", STARindex);
	noshell::Exit e2=noshell::C(STARpath, "--runThreadN", starthreads, "--runMode", "genomeGenerate", "--genomeDir", STARindex, 
		"--genomeFastaFiles", OutPrefix+"/juncseq.fa", "--sjdbGTFfile", OutPrefix+"/juncannot.gtf");
	if(!e2.success()){
		cout<<"STAR index cannot be generated successfully\n";
		return -1;
	}
	noshell::Exit e3=noshell::C("mv", "Log.out", STARindex);
	noshell::Exit e4=noshell::C("mkdir", "-p", OutPrefix+"/StarAlign/");
	noshell::Exit e5=noshell::C(STARpath, "--runThreadN", starthreads, "--genomeDir", STARindex, "--readFilesIn", ReadPrefix+"_1.fastq", 
		ReadPrefix+"_2.fastq", "--outFileNamePrefix", OutPrefix+"/StarAlign/", "--outSAMtype", "BAM", "SortedByCoordinate", "--outSAMstrandField", "intronMotif");
	if(!e5.success()){
		cout<<"STAR alignment fail\n";
		return -1;
	}
	return 0;
};

int Assembly(string execpath, string bamfile, string OutPrefix){
	if(execpath.find("stringtie")!=string::npos){
		noshell::Exit e=noshell::C(execpath, bamfile, "-o", OutPrefix+"/StarAlign/genes.gtf");
		if(!e.success()){
			cout<<"scallop failed\n";
			return -1;
		}
	}
	else if(execpath.find("scallop")!=string::npos){
		noshell::Exit e=noshell::C(execpath, "-i", bamfile, "-o", OutPrefix+"/StarAlign/genes.gtf", "--verbose", 0);
		if(!e.success()){
			cout<<"stringtie failed\n";
			return -1;
		}
	}
	return 0;
};

int ConvertCoord(Transcript_t& trans, Transcript_t& newtrans, SV_t sv){
	sort(trans.vExon.begin(), trans.vExon.end());
	if(trans.Strand=='-')
		reverse(trans.vExon.begin(), trans.vExon.end());
	int range1=sv.BP1.EndPos-sv.BP1.StartPos;
	newtrans.TransID=trans.TransID;
	newtrans.GeneName=trans.GeneName;
	newtrans.Chr=trans.Chr;
	newtrans.Strand=trans.Strand;
	newtrans.TxStart=trans.TxStart;
	newtrans.TxEnd=trans.TxEnd;
	int sep=0;
	for(vector<Interval_t>::iterator it=trans.vExon.begin(); it!=trans.vExon.end(); it++){
		if(it->EndPos<=range1){
			Interval_t tmpexon(sv.BP1.Chr, sv.BP1.StartPos+it->StartPos, sv.BP1.StartPos+it->EndPos, it->Strand, it->TransID);
			if(sv.BP1.IsLeft){
				tmpexon.EndPos=sv.BP1.StartPos+range1-it->StartPos;
				tmpexon.StartPos=sv.BP1.StartPos+range1-it->EndPos;
				if(tmpexon.Strand=='+')
					tmpexon.Strand='-';
				else if(tmpexon.Strand=='-')
					tmpexon.Strand='+';
			}
			newtrans.vExon.push_back(tmpexon);
			if(it!=trans.vExon.begin() && (it-1)->StartPos>=range1)
				sep=(int)newtrans.vExon.size();
		}
		else if(it->StartPos>=range1){
			Interval_t tmpexon(sv.BP2.Chr, sv.BP2.StartPos+it->StartPos-range1, sv.BP2.StartPos+it->EndPos-range1, it->Strand, it->TransID);
			if(sv.BP2.IsLeft){
				tmpexon.EndPos=sv.BP2.EndPos-(it->StartPos-range1);
				tmpexon.StartPos=sv.BP2.EndPos-(it->EndPos-range1);
				if(tmpexon.Strand=='+')
					tmpexon.Strand='-';
				else if(tmpexon.Strand=='-')
					tmpexon.Strand='+';
			}
			newtrans.vExon.push_back(tmpexon);
			if(it!=trans.vExon.begin() && (it-1)->EndPos<=range1)
				sep=(int)newtrans.vExon.size();
		}
		else{
			Interval_t tmpexon(sv.BP1.Chr, sv.BP1.StartPos+it->StartPos, sv.BP1.EndPos, it->Strand, it->TransID);
			if(sv.BP1.IsLeft){
				tmpexon.EndPos=sv.BP1.StartPos+range1-it->StartPos;
				tmpexon.StartPos=sv.BP1.StartPos;
				if(tmpexon.Strand=='+')
					tmpexon.Strand='-';
				else if(tmpexon.Strand=='-')
					tmpexon.Strand='+';
			}
			Interval_t tmpexon2(sv.BP2.Chr, sv.BP2.StartPos, sv.BP2.StartPos+it->EndPos-range1, it->Strand, it->TransID);
			if(sv.BP2.IsLeft){
				tmpexon2.EndPos=sv.BP2.EndPos;
				tmpexon2.StartPos=sv.BP2.EndPos-(it->EndPos-range1);
				if(tmpexon2.Strand=='+')
					tmpexon2.Strand='-';
				else if(tmpexon2.Strand=='-')
					tmpexon2.Strand='+';
			}
			if(it->StartPos<range1){
				newtrans.vExon.push_back(tmpexon);
				sep=(int)newtrans.vExon.size();
				newtrans.vExon.push_back(tmpexon2);
			}
			else{
				newtrans.vExon.push_back(tmpexon2);
				sep=(int)newtrans.vExon.size();
				newtrans.vExon.push_back(tmpexon);
			}
		}
	}
	return sep;
};

void Change2ClosestExon(vector<Interval_t>::iterator firstiterator, vector<Interval_t>::iterator lastiterator, int bppos, vector<Transcript_t>& genes){
	int thresh=50;
	for(vector<Interval_t>::iterator it=firstiterator; it!=lastiterator; it++){
		bool flag=false;
		Interval_t exon;
		double ratio=0;
		if(it->StartPos!=bppos && it->EndPos!=bppos){
			for(int i=0; i<genes.size(); i++)
				for(vector<Interval_t>::iterator itexon=genes[i].vExon.begin(); itexon!=genes[i].vExon.end(); itexon++){
					int overlaplength=0;
					if(it->Chr==itexon->Chr && it->StartPos<=itexon->StartPos)
						overlaplength=min(itexon->EndPos-itexon->StartPos, it->EndPos-itexon->StartPos);
					else if(it->Chr==itexon->Chr && it->StartPos>itexon->StartPos)
						overlaplength=min(itexon->EndPos-it->StartPos, it->EndPos-it->StartPos);
					if(overlaplength<0)
						overlaplength=0;
					int tmpratio=1.0*overlaplength/(itexon->EndPos-itexon->StartPos)+1.0*overlaplength/(it->EndPos-it->StartPos);
					if(tmpratio>ratio){
						exon=(*itexon);
						flag=true;
					}
				}
			// check before and after exons, if overlap with modified one, modify them instead.
			if(!flag)
				continue;
			it->StartPos=exon.StartPos;
			it->EndPos=exon.EndPos;
			if(it!=firstiterator && (it-1)->Chr==it->Chr && !((it->StartPos<(it-1)->StartPos)==(it->EndPos<(it-1)->StartPos) && (it->StartPos<(it-1)->StartPos)==(it->StartPos<(it-1)->EndPos))){
				if((it-1)->StartPos>=it->StartPos && (it-1)->EndPos<=it->EndPos)
					(it-1)->EndPos=it->StartPos;
				else if((it-1)->StartPos<it->StartPos)
					(it-1)->EndPos=it->StartPos;
				else
					(it-1)->StartPos=it->EndPos;
			}
			if((it+1)!=lastiterator && (it+1)->Chr==it->Chr && !((it->StartPos<(it+1)->StartPos)==(it->EndPos<(it+1)->StartPos) && (it->StartPos<(it+1)->StartPos)==(it->StartPos<(it+1)->EndPos))){
				if((it+1)->StartPos>=it->StartPos && (it+1)->EndPos<=it->EndPos)
					(it+1)->EndPos=it->StartPos;
				else if((it+1)->StartPos<it->StartPos)
					(it+1)->EndPos=it->StartPos;
				else
					(it+1)->StartPos=it->EndPos;
			}
		}
		else{
			for(int i=0; i<genes.size(); i++)
				for(vector<Interval_t>::iterator itexon=genes[i].vExon.begin(); itexon!=genes[i].vExon.end(); itexon++){
					if(itexon->Chr==it->Chr && abs(itexon->StartPos-it->StartPos)<thresh && it->StartPos!=bppos){
						it->StartPos=itexon->StartPos;
						// check whether changing conflict with other exons
						if(it!=firstiterator && (it-1)->Chr==it->Chr && !((it->StartPos<(it-1)->StartPos)==(it->EndPos<(it-1)->StartPos) && (it->StartPos<(it-1)->StartPos)==(it->StartPos<(it-1)->EndPos)))
							(it-1)->EndPos=it->StartPos;
						if(it!=firstiterator && (it+1)->Chr==it->Chr && !((it->StartPos<(it+1)->StartPos)==(it->EndPos<(it+1)->StartPos) && (it->StartPos<(it+1)->StartPos)==(it->StartPos<(it+1)->EndPos)))
							(it+1)->EndPos=it->StartPos;
					}
					if(itexon->Chr==it->Chr && abs(itexon->EndPos-it->EndPos)<thresh && it->EndPos!=bppos){
						it->EndPos=itexon->EndPos;
						// check whether changing conflict with other exons
						if(it!=firstiterator && (it-1)->Chr==it->Chr && !((it->StartPos<(it-1)->StartPos)==(it->EndPos<(it-1)->StartPos) && (it->StartPos<(it-1)->StartPos)==(it->StartPos<(it-1)->EndPos)))
							(it-1)->StartPos=it->EndPos;
						if(it!=firstiterator && (it+1)->Chr==it->Chr && !((it->StartPos<(it+1)->StartPos)==(it->EndPos<(it+1)->StartPos) && (it->StartPos<(it+1)->StartPos)==(it->StartPos<(it+1)->EndPos)))
							(it+1)->StartPos=it->EndPos;
					}
				}
		}
	}
};

void DetermineGTF(Transcript_t& trans, int sep, SV_t& sv, bool potential1, bool potential2, vector<Transcript_t>& vTrans, ofstream& output){
	if(potential1){
		vector<Transcript_t> genes1=GeneBinarySearch(vTrans.begin(), vTrans.end(), sv.BP1);
		int bppos=(sv.BP1.IsLeft)?sv.BP1.StartPos:sv.BP1.EndPos;
		Change2ClosestExon(trans.vExon.begin(), trans.vExon.begin()+sep, bppos, genes1);
	}
	if(potential2){
		vector<Transcript_t> genes2=GeneBinarySearch(vTrans.begin(), vTrans.end(), sv.BP2);
		int bppos=(!sv.BP2.IsLeft)?sv.BP2.StartPos:sv.BP2.EndPos;
		Change2ClosestExon(trans.vExon.begin(), trans.vExon.begin()+sep, bppos, genes2);
	}
	Transcript_t tmptrans(trans.TransID, trans.GeneName, trans.Chr, trans.Strand, trans.TxStart, trans.TxEnd);
	for(vector<Interval_t>::iterator it=trans.vExon.begin(); it!=trans.vExon.end(); it++)
		if(it->EndPos-it->StartPos>3)
			tmptrans.vExon.push_back(*it);
	output<<(tmptrans.vExon[0].Chr)<<"\tTSVtrans\ttranscript\t"<<(tmptrans.vExon[0].StartPos)<<"\t"<<(tmptrans.vExon[0].EndPos)<<"\t.\t"<<(tmptrans.vExon[0].Strand)<<"\t.\t";
	output<<"gene_id \"FG"<<"\"; transcript_id \"FGtrans\";\n";
	for(int i=0; i<tmptrans.vExon.size(); i++){
		const Interval_t exon=tmptrans.vExon[i];
		output<<(exon.Chr)<<"\tTSVtrans\texon\t"<<(exon.StartPos)<<"\t"<<(exon.EndPos)<<"\t.\t"<<(exon.Strand)<<"\t.\t";
		output<<"gene_id \"FG"<<"\"; transcript_id \"FGtrans\"; exon_number \""<<i<<"\"\n";
	}
};

int FinalJunction(vector<Transcript_t>& vTrans, vector<SV_t>& JunctionRegion, vector<bool>& Potent1, vector<bool>& Potent2, string OutPrefix){
	string newtransfile=OutPrefix+"/StarAlign/genes.gtf";
	ifstream input(newtransfile);
	ofstream output(OutPrefix+"/nonfg.gtf", ios::out);
	string line;
	Transcript_t tmptrans;
	bool keeptrans=false;
	while(getline(input, line)){
		if(line[0]=='>')
			continue;
		vector<string> strs;
		boost::split(strs, line, boost::is_any_of("\t"));
		if(strs[2]=="transcript"){
			if(tmptrans.Chr!=""){
				Transcript_t tmptransconvert;
				int sep=ConvertCoord(tmptrans, tmptransconvert, JunctionRegion[stoi(tmptrans.Chr)-1]);
				DetermineGTF(tmptransconvert, sep, JunctionRegion[stoi(tmptrans.Chr)-1], Potent1[stoi(tmptrans.Chr)-1], Potent2[stoi(tmptrans.Chr)-1], vTrans, output);
			}
			tmptrans.Chr="";
			int thischr=stoi(strs[0])-1;
			int bp1length=JunctionRegion[thischr].BP1.EndPos-JunctionRegion[thischr].BP1.StartPos;
			if(stoi(strs[3])<=bp1length && stoi(strs[4])>=bp1length){
				keeptrans=true;
				tmptrans.Chr=strs[0];
				tmptrans.TxStart=stoi(strs[3]);
				tmptrans.TxEnd=stoi(strs[4]);
				tmptrans.Strand=strs[6][0];
				tmptrans.vExon.clear();
			}
			else
				keeptrans=false;
		}
		else if(keeptrans){
			string TID=Transcript_t::GetFeature(line, "transcript_id");
			Interval_t tmpexon(strs[0], stoi(strs[3])-1, stoi(strs[4]), strs[6][0], TID);
			tmptrans.vExon.push_back(tmpexon);
		}
	}
	input.close();
	output.close();
	return 0;
};

string genomefile;
string SVpredfile;
string annotationfile;
string OutPrefix;
string BamPrefix;
string FqPrefix;
string ConcordBam="Aligned.sortedByCoord.out.bam";
string ChimericBam="Chimeric.out.bam";
string STARpath="STAR";
int starthreads=8;
string execpath="scallop";

int ParseArguments(int argc, char* argv[]){
	int i=1;
	while(i<argc){
		if(string(argv[i])=="-g" && i<argc-1){
			genomefile=(string)argv[i+1];
			i+=2;
		}
		else if(string(argv[i])=="-p" && i<argc-1){
			SVpredfile=(string)argv[i+1];
			i+=2;
		}
		else if(string(argv[i])=="-a" && i<argc-1){
			annotationfile=(string)argv[i+1];
			i+=2;
		}
		else if(string(argv[i])=="-o" && i<argc-1){
			OutPrefix=(string)argv[i+1];
			i+=2;
		}
		else if(string(argv[i])=="-b" && i<argc-1){
			BamPrefix=(string)argv[i+1];
			i+=2;
		}
		else if(string(argv[i])=="-q" && i<argc-1){
			FqPrefix=(string)argv[i+1];
			i+=2;
		}
		else
			return -1;
	}
	return 0;
};

int main(int argc, char* argv[]){

	if(argc==1)
		printf("SVTranscript -g <genomefile> -p <SVpredfile> -a <annotationfile> -o <OutPrefix> -b <BamPrefix> -q <FqPrefix>\n");
	else{
		if(ParseArguments(argc, argv)!=0)
			return -1;

		noshell::Exit e=noshell::C("mkdir", "-p", OutPrefix);
		if(!e.success()){
			cout<<"fail to create directory "<<OutPrefix<<endl;
			return -1;
		}

		map<string,string> Genome;
		ReadGenome(genomefile, Genome);

		vector<SV_t> SVs;
		ReadBedpeSV(SVpredfile, SVs);

		vector<Transcript_t> vTrans;
		ReadGTF(annotationfile, vTrans);

		vector<SV_t> JunctionRegion;
		vector<bool> Potent1;
		vector<bool> Potent2;
		InitialJunction(SVs, vTrans, Genome, OutPrefix, JunctionRegion, Potent1, Potent2);

		int status=RegionalReadWritter(BamPrefix, ConcordBam, ChimericBam, FqPrefix, OutPrefix+"/JuncReads", JunctionRegion);
		if(status!=0){
			cout<<"something wrong in writing reads within junction region\n";
			return -1;
		}

		string STARindex=OutPrefix+"/STARindex";
		status=Aligning(STARpath, STARindex, starthreads, OutPrefix+"/JuncReads", OutPrefix);
		if(status!=0)
			return -1;
		status=Assembly(execpath, OutPrefix+"/StarAlign/"+ConcordBam, OutPrefix);
		if(status!=0)
			return -1;

		status=FinalJunction(vTrans, JunctionRegion, Potent1, Potent2, OutPrefix);
	}
	return 0;
}