#include "McTruthInfo.hh"

#include <iostream>
#include <cassert>
#include "classes/DelphesClasses.h"

#include "external/mt2analysis/Utilities.hh"


using namespace std;

McTruthInfo::McTruthInfo() :
  genParticles_(0)
{
}

McTruthInfo::~McTruthInfo() {
}

void McTruthInfo::Dump() {

  cout<<" ~~~~~~~~~~~"<<endl;
  for (int k = 0 ; k<genParticles_->GetEntries(); k++) {
    GenParticle * c =(GenParticle*) genParticles_->At(k);
    if (c==0) continue;
    cout<<k<<"\t"<<c->PID<<" "<<c->M1<<" "<<" "<<c->D1<<"\t"<<c->Mass<<endl;
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

/* look for a given particle in SUSY cascades */
int McTruthInfo::findPinSusy(int pidToFind) {

  int np=0;

  for (int k = 0 ; k<genParticles_->GetEntries(); k++) {
    GenParticle * c =(GenParticle*) genParticles_->At(k);
    if (c==0) continue; //try to prevent crashes....
    int pid = std::abs(c->PID);

    if (pid==pidToFind) { 
      //does this particle come from SUSY?
      if (  checkMom(k,1000000) ) np++;
    }
  }
  return np;

}

int McTruthInfo::countTrueLeptons(leptonType lt) {

  int n=0;
  //easy as pie, except that depending on what is stored, we might actually
  //want to restrict ourselves to status 3 or the Pythia 8 equivalent
  for (int k = 0 ; k<genParticles_->GetEntries(); k++) {
    GenParticle * c =(GenParticle*) genParticles_->At(k);
    if (c==0) continue; //try to prevent crashes....
    int pid = std::abs(c->PID);

    if (pid==11 && lt==kElMu) n++;
    else if (pid==13 && lt==kElMu) n++;
    else if (pid==15 && lt==kTau) n++;
  }

  return n;
}

/*
look for chi_20 -> slepton lepton -> lepton chi_10 lepton
code is 10*nStaus + nSElectron+Smuon
should allow for a simple cut that removes events with taus or
multiple edge decays while still retaining that info in case it is of interest
 */
int McTruthInfo::findChi2ToChi1() {

  int code = 0;

  for (int k = 0 ; k<genParticles_->GetEntries(); k++) {
    GenParticle * c =(GenParticle*) genParticles_->At(k);
    if (c==0) continue; //try to prevent crashes....
    int pid = std::abs(c->PID);

    if (pid==1000022) { //chi10
      //if the mom a slepton?
      bool el = checkMom(k,1000011,0);
      bool mu=checkMom(k,1000013,0);
      bool tau = checkMom(k,1000015,0);
      if (el||mu||tau ) {
	//is the grandmom a chi_20?
	if (checkMom(k,1000023,1) ) {
	  if (el) ++code;
	  else if (mu) ++code;
	  else if (tau) code += 10;
	  else assert(0);
	}
      }
    }
  }
  return code;
}

int McTruthInfo::getSusyProductionProcess() {

  vector<int> susyMoms = findSusyMoms();

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

vector<int> McTruthInfo::findSusyMoms() {
  //look for SUSY particles whose moms are not SUSY particles
  //these are the produced particles in the hard-scatter
  vector<int> susyMoms;

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

bool McTruthInfo::checkMom(int index, int PidToLookFor, int recursionsLeft) {
  //recursive (or not) check of whether the mom or any grandmom (if recursive)
  // of the particle at index is of PID PidToLoopFor

  //special case: PidToLookFor of 1000000 means 'any SUSY particle'

  //get mom index
  GenParticle * theCand =(GenParticle*) genParticles_->At(index);
  if (theCand==0) return false; //TClonesArray::At returns 0 if index is not found
  int momIndex1 = theCand->M1;
  if (momIndex1 >=0 ) {
    GenParticle * theMom = (GenParticle*)genParticles_->At(momIndex1);
    if (theMom==0) return false; //TClonesArray::At returns 0 if index is not found
    if ( std::abs(theMom ->PID)==PidToLookFor) return true;
    else if (PidToLookFor==1000000 && isSusy(theMom->PID)) return true;
    else  {
      if (recursionsLeft>0) return checkMom( momIndex1,PidToLookFor,recursionsLeft-1);
      else return false;
    }
  }

  //if momIndex is <0
  return false;

}

int McTruthInfo::GetTtbarDecayCode() {
  //classify ttbar events using my old coding scheme
  
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


// *** DY truth info
//find Z->ll
//return: gen-level P4, flavor
vector< pair< TLorentzVector, int> > McTruthInfo::GetDYTruth() {

  //this method assumes Z comes first in list (which it should)
  //we ignore Z decays other than e,mu

  int Zindex=-99;
  vector<  pair< TLorentzVector, int> > Zdaughters;
  for (int k = 0 ; k<genParticles_->GetEntries(); k++) {
    GenParticle * c =(GenParticle*) genParticles_->At(k);
    if (c==0) continue;
    if ( c->PID ==23 ) {    //look for Z
      Zindex = k;
    }
    if ( c->M1 == Zindex) {    //look for daughters of Z
      int pid = std::abs(c->PID);
      if (pid == 11 || pid == 13) {
	pair<TLorentzVector,int> thisDaughter = make_pair( c->P4(), c->PID);
	Zdaughters.push_back(thisDaughter);
      }
    }

    if (Zdaughters.size()==2) break;
  }
  return Zdaughters;
}

//input: (1) list of reco e and mu
// (2) output of GetDYTruth()
//output: vector of is matched or not, given as indices on the e/mu list
vector<int> McTruthInfo::MatchDYRecoGen(const vector< pair< TLorentzVector, int> > & genDY,
					TClonesArray * els,TClonesArray * mus) {

  vector<int> recomatches;

  //loop over gen-level particles
  for (size_t igen = 0; igen<genDY.size(); ++igen) {

    TLorentzVector  genP4 = genDY[igen].first;
    //check if particle is e or mu
    int  flavor = std::abs(genDY[igen].second);

    //loop over e/mu list
    int nr=0;
    if (flavor==11) nr = els->GetEntries();
    else if (flavor==13) nr=mus->GetEntries();
    else assert(0);

    int matchindex=-1;
    for (int ireco=0;ireco<nr;ireco++) {
      //look for DR (and pT?) match
      TLorentzVector * recoP4=0;
      if (flavor==11) recoP4 = new TLorentzVector( ((Electron*)els->At(ireco))->P4());
      else if (flavor==13) recoP4 = new TLorentzVector( ((Muon*)mus->At(ireco))->P4());

      double dr = Util::GetDeltaR(genP4.Eta(), recoP4->Eta(), genP4.Phi(), recoP4->Phi());
      //it appears that there is ~no resolution smearing in Delphes for
      //phi and eta, so dr match is sufficient, and can be very tight
      //leave loose pt match, but not really necessary
      //(it appears that pT tends to be well within 5%, so this is very loose)
      if (dr<0.04 && std::abs(1- genP4.Pt()/recoP4->Pt() )<0.3) {
	matchindex=ireco;
	//debug
	//	cout<<"[DY match] "<<flavor<<"\t"<<dr<<" "<< 1- genP4.Pt()/recoP4->Pt() <<endl;
	delete recoP4;
	break;
      }

      delete recoP4;
    }
    //if found, push_back index on e/mu list to output vector
    //if not found, pushback dummy index
    recomatches.push_back(matchindex);

  }
  return recomatches;

}

