#include "McTruthInfo.hh"

#include <iostream>
//#include <cassert>
#include "classes/DelphesClasses.h"

using namespace std;

McTruthInfo::McTruthInfo() :
  genParticles_(0)
{
}

McTruthInfo::~McTruthInfo() {
}

bool McTruthInfo::checkMom(int index, int PidToLookFor) {

  //get mom index
  int momIndex1 = ((GenParticle*)genParticles_->At(index))->M1;
  if (momIndex1 >=0 ) {
    if ( std::abs(((GenParticle*)genParticles_->At(momIndex1))->PID)==PidToLookFor) return true;
    else  return checkMom( momIndex1,PidToLookFor);
  }

  //if momIndex is <0
  return false;

}

int McTruthInfo::GetTtbarDecayCode(TClonesArray* genParticles) {

  if (genParticles!=0)     genParticles_=genParticles;
  
  int ne=0, nmu=0, ntau=0, nq=0,nb=0;
  //  cout<<" ~~~~~~~~~~~"<<endl;
  for (int k = 0 ; k<genParticles_->GetEntries(); k++) {
    GenParticle * c =(GenParticle*) genParticles_->At(k);
    int pid = std::abs(c->PID);
    //    cout<<c->PID<<" "<<c->M1<<endl;
    if (pid == 11 || pid == 13 || pid==15 || pid==1 ||pid==2 ||pid==3||pid==4 ||pid==5) { //find ultimate top->W decay product
      bool      topIsUltimateMom = checkMom(k, 6);

      if (topIsUltimateMom) {
	if      (pid ==11) ne++;
	else if (pid ==13) nmu++;
	else if (pid ==15) ntau++;
	else if (pid==5) nb++;
	else  nq++;
      }
    }
  }

  nq/=2;

  //  cout<<ne<<" "<<nmu<<" "<<ntau<<" "<<nq<<endl;

  //here's my old coding scheme -- guess i'll stick with it!
  //  enum TopDecayCategory {kTTbarUnknown=0,kAllLeptons=1,kAllHadronic=2,kOneElectron=3,kOneMuon=4,kOneTauE=5,kOneTauMu=6,kOneTauHadronic=7,kAllTau=8,kTauPlusLepton=9, nTopCategories=10};

  if       ( ne+nmu==2 ) return 1;
  else if  ( nq == 2)    return 2;
  else if  ( ne==1 && (ne+nmu+ntau==1) ) return 3;
  else if  ( nmu==1 && (ne+nmu+ntau==1)) return 4;
  else if  ( ntau==1&& (ne+nmu+ntau==1)) return 5;
  //no codes 6,7 because i don't have tau decay info right now (code 5 and 6 and 7 are merged)
  else if  ( ntau==2) return 8;
  else if  ( ntau==1 && ne+nmu==1) return 9;

  return 0;

}
