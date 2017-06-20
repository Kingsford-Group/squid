#include "TRA.h"

using namespace std;

void TRA_t::DecideSwap(){
    if(Ref1>Ref2){
        int tmp[2]={Ref1, Pos1};
        DirType tmpdt=DT1;
        Ref1=Ref2; Pos1=Pos2; DT1=DT2;
        Ref2=tmp[0]; Pos2=tmp[1]; DT2=tmpdt;
    }
};

pair<int,int> TRA_t::UpdatePoint(map<int, int>& RefLength, pair<int,int> BP, DirType DT){ //BP is the apri of reference and position, DT is the part supported by reads(DT1 or !DT2)
    if(DT1==1 && DT2==0){
        if(BP.first==Ref1 && (BP.second<Pos1 || (BP.second==Pos1 && DT==0)))
            return make_pair(Ref2, Pos2+Pos1-BP.second-1);
        else if(BP.first==Ref1 && (BP.second>Pos1 || (BP.second==Pos1 && DT==1)))
            return make_pair(Ref1, RefLength[Ref2]-Pos2+BP.second-Pos1);
        else if(BP.first==Ref2 && (BP.second>Pos2 || (BP.second==Pos2 && DT==0)))
            return make_pair(Ref1, RefLength[Ref2]-BP.second-1);
        else
            return BP;
    }
    else if(DT1==1 && DT2==1){
        if(BP.first==Ref1 && (BP.second<Pos1 || (BP.second==Pos1 && DT==0)))
            return make_pair(Ref2, BP.second);
        else if(BP.first==Ref1 && (BP.second>Pos1 || (BP.second==Pos1 && DT==1)))
            return make_pair(Ref1, Pos2+BP.second-Pos1);
        else if(BP.first==Ref2 && (BP.second<Pos2 || (BP.second==Pos2 && DT==0)))
            return make_pair(Ref1, BP.second);
        else if(BP.first==Ref2 && (BP.second>Pos2 || (BP.second==Pos2 && DT==1)))
            return make_pair(Ref2, Pos1+BP.second-Pos2);
        else
            return BP;
    }
    else if(DT1==0 && DT2==0){
        if(BP.first==Ref1 && (BP.second>Pos1 || (BP.second==Pos1 && DT==1)))
            return make_pair(Ref2, Pos2+BP.second-Pos1);
        else if(BP.first==Ref2 && (BP.second>Pos2 || (BP.second==Pos2 && DT==1)))
            return make_pair(Ref1, Pos1+BP.second-Pos2);
        else
            return BP;
    }
    else{
        if(BP.first==Ref1 && (BP.second>Pos1 || (BP.second==Pos1 && DT==1)))
            return make_pair(Ref2, RefLength[Ref1]-BP.second-1);
        else if(BP.first==Ref2 && (BP.second<Pos2 || (BP.second==Pos2 && DT==0)))
            return make_pair(Ref1, Pos1+Pos2-BP.second-1);
        else if(BP.first==Ref2 && (BP.second>Pos2 || (BP.second==Pos2 && DT==1)))
            return make_pair(Ref2, RefLength[Ref1]-Pos1+BP.second-Pos2);
        else return BP;
    }
};

void TRA_t::UpdateTRA(map<int, int>& RefLength, vector<TRA_t>::iterator ittra){
    if(!(ittra->Ref1==Ref1 && ittra->Pos1==Pos1 && ittra->Ref2==Ref2 && ittra->Pos2==Pos2)){
        pair<int,int> tmp1=UpdatePoint(RefLength, make_pair(ittra->Ref1, ittra->Pos1), ittra->DT1);
        pair<int,int> tmp2=UpdatePoint(RefLength, make_pair(ittra->Ref2, ittra->Pos2), (DirType)(((int)ittra->DT2+1)%2));
        ittra->Ref1=tmp1.first; ittra->Pos1=tmp1.second;
        ittra->Ref2=tmp2.first; ittra->Pos2=tmp2.second;
        ittra->DecideSwap();
    }
};

void TRA_t::EditnReverse(map<int, int>& RefLength){
    if(DT1==1 && DT2==0){
        int partlength[4]={Pos1, RefLength[Ref1]-Pos1, Pos2, RefLength[Ref2]-Pos2};
        Pos1=partlength[3]; Pos2=partlength[2];
        RefLength[Ref1]=partlength[1]+partlength[3];
        RefLength[Ref2]=partlength[0]+partlength[2];
    }
    else if(DT1==1 && DT2==1){
        int partlength[4]={Pos1, RefLength[Ref1]-Pos1, Pos2, RefLength[Ref2]-Pos2};
        Pos1=partlength[2]; Pos2=partlength[0];
        RefLength[Ref1]=partlength[1]+partlength[2];
        RefLength[Ref2]=partlength[0]+partlength[3];
    }
    else if(DT1==0 && DT2==0){
        int partlength[4]={Pos1, RefLength[Ref1]-Pos1, Pos2, RefLength[Ref2]-Pos2};
        RefLength[Ref1]=partlength[0]+partlength[3];
        RefLength[Ref2]=partlength[1]+partlength[2];
    }
    else{
        int partlength[4]={Pos1, RefLength[Ref1]-Pos1, Pos2, RefLength[Ref2]-Pos2};
        Pos1=partlength[0]; Pos2=partlength[1];
        RefLength[Ref1]=partlength[0]+partlength[2];
        RefLength[Ref2]=partlength[1]+partlength[3];
    }
};
