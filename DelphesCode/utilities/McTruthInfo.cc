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

void McTruthInfo::Dump(TClonesArray* genParticles) {
  if (genParticles!=0)     genParticles_=genParticles;

  cout<<" ~~~~~~~~~~~"<<endl;
  for (int k = 0 ; k<genParticles_->GetEntries(); k++) {
    GenParticle * c =(GenParticle*) genParticles_->At(k);
    if (c==0) continue;
    cout<<c->PID<<" "<<c->M1<<endl;
  }

}

bool McTruthInfo::isSusy(int pid) {
  //is this particle a SUSY particle or not?
  pid = std::abs(pid); //just in case

  if (pid>=1000001 && pid<=1000039) return true;
  if (pid>=2000001 && pid<=2000039) return true;
  //indices only go up to 2000015 but going to 39 seems ok too

  return false;
}

int McTruthInfo::getSusyProductionProcess(TClonesArray* genParticles) {

  if (genParticles!=0)     genParticles_=genParticles;

  vector<int> susyMoms = findSusyMoms(genParticles);

  int n_slepton = 0;
  int n_ewkino = 0;
  int n_gluino=0;
  int n_stop = 0;
  int n_sbottom = 0;
  int n_squark = 0;
  int n_other=0;

  for (size_t i=0; i<susyMoms.size(); i++) {
    int pid = std::abs(susyMoms[i]);
    if ( pid >= 1000001 && pid <= 1000004) n_squark++;
    else if ( pid == 1000005 || pid==2000005) n_sbottom++;
    else if ( pid == 1000006 || pid==2000006) n_stop++;
    else if ( pid >= 1000011 && pid <= 1000016) n_slepton++;
    else if ( pid >= 2000011 && pid <= 2000016) n_slepton++;
    else if ( pid >= 2000001 && pid <= 2000004) n_squark++;
    else if ( pid==1000021) n_gluino++;
    else if ( pid>=1000012 && pid <= 1000037) n_ewkino++;
    else n_other++;
  }

  return   n_slepton*1000000
    +    n_ewkino*100000
    +    n_gluino*10000
    +    n_stop *1000
    +    n_sbottom*100
    +    n_squark *10
    +    n_other*1;
  
}

vector<int> McTruthInfo::findSusyMoms(TClonesArray* genParticles) {
  //look for SUSY particles whose moms are not SUSY particles
  //these are the produced particles in the hard-scatter
  vector<int> susyMoms;

  if (genParticles!=0)     genParticles_=genParticles;

  for (int k = 0 ; k<genParticles_->GetEntries(); k++) {
    GenParticle * c =(GenParticle*) genParticles_->At(k);
    if (c==0) continue; //try to prevent crashes....
    int pid = std::abs(c->PID);

    int momIndex1 = c->M1;
    if (momIndex1<0) continue;
    GenParticle * theMom = (GenParticle*)genParticles_->At(momIndex1);
    int momPid = std::abs(theMom->PID);
    if (isSusy(pid) && !isSusy(momPid)) {
      susyMoms.push_back(c->PID);
    }
    //could save time by cutting off the loop after vector size is 2
  }
  return susyMoms;

}

bool McTruthInfo::checkMom(int index, int PidToLookFor) {
  //recursive check of whether the mom or any grandmom of the particle at index
  //is of PID PidToLoopFor

  //get mom index
  GenParticle * theCand =(GenParticle*) genParticles_->At(index);
  if (theCand==0) return false; //TClonesArray::At returns 0 if index is not found
  int momIndex1 = theCand->M1;
  if (momIndex1 >=0 ) {
    GenParticle * theMom = (GenParticle*)genParticles_->At(momIndex1);
    if (theMom==0) return false; //TClonesArray::At returns 0 if index is not found
    if ( std::abs(theMom ->PID)==PidToLookFor) return true;
    else  return checkMom( momIndex1,PidToLookFor);
  }

  //if momIndex is <0
  return false;

}

int McTruthInfo::GetTtbarDecayCode(TClonesArray* genParticles) {
  //classify ttbar events using my old coding scheme
  if (genParticles!=0)     genParticles_=genParticles;
  
  int ne=0, nmu=0, ntau=0, nq=0,nb=0;
  //  cout<<" ~~~~~~~~~~~"<<endl;
  for (int k = 0 ; k<genParticles_->GetEntries(); k++) {
    GenParticle * c =(GenParticle*) genParticles_->At(k);
    if (c==0) continue; //try to prevent crashes....
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

  return 0; //uncategorized

}
