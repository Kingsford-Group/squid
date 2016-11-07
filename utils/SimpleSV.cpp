#include "SimpleSV.h"

using namespace std;

void SimpleSV_t::DecideSwap(){
    if(StartPos>EndPos){
        int tmp=StartPos;
        StartPos=EndPos; EndPos=tmp;
    }
};

pair<int,int> SimpleSV_t::UpdatePoint(pair<int,int> BP){
    if(Type==INS && BP.first==RefID && BP.second>=StartPos)
        return make_pair(BP.first, BP.second+EndPos);
    else if(Type==INV && BP.first==RefID && BP.second>=StartPos && BP.second<EndPos)
        return make_pair(BP.first, StartPos+EndPos-BP.second);
    else if(Type==DEL && RefID==BP.first && BP.second>=EndPos)
        return make_pair(BP.first, BP.second-(EndPos-StartPos));
    else
        return BP;
};

void SimpleSV_t::UpdateSimpleSV(vector<SimpleSV_t>::iterator itssv){
    if(!(itssv->StartPos==StartPos && itssv->EndPos==EndPos)){
        itssv->StartPos=(UpdatePoint(make_pair(itssv->RefID, itssv->StartPos))).second;
        if(itssv->Type!=INS)
            itssv->EndPos=(UpdatePoint(make_pair(itssv->RefID, itssv->EndPos))).second;
    }
};

void SimpleSV_t::EditnReverse(map<int, int>& RefLength){
    if(Type==INS){
        RefLength[RefID]+=EndPos;
        Type=DEL;
        EndPos+=StartPos;
    }
    else if(Type==DEL){
        RefLength[RefID]-=EndPos-StartPos;
        Type=INS;
        EndPos-=StartPos;
    }
};

void SimpleSV_t::UpdateTRA(vector<TRA_t>::iterator ittra){
    ittra->Pos1=(UpdatePoint(make_pair(ittra->Ref1,ittra->Pos1))).second;
    ittra->Pos2=(UpdatePoint(make_pair(ittra->Ref2,ittra->Pos2))).second;
};

void TRA_t::UpdateSimpleSV(map<int, int>& RefLength, vector<SimpleSV_t>::iterator itssv){
    itssv->StartPos=UpdatePoint(RefLength, make_pair(itssv->RefID,itssv->StartPos),(DirType)0).second;
    if(itssv->Type!=INS)
        itssv->EndPos=UpdatePoint(RefLength, make_pair(itssv->RefID,itssv->EndPos),(DirType)0).second;
};
