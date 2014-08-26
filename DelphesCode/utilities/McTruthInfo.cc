#include "McTruthInfo.hh"

#include <iostream>
#include <cassert>

//#include "external/mt2analysis/Utilities.hh"


using namespace std;

McTruthInfo::McTruthInfo() :
  genParticles_(0),
  nGenParticles_(0)
{
  //  setNloKFactor(20000,1.7);
  //can hard-code these here or let the user of the class set them
}

McTruthInfo::~McTruthInfo() {
}

void McTruthInfo::Dump() {

  cout<<" ~~~~~~~~~~~"<<endl;
  for (int k = 0 ; k<nGenParticles_; k++) {
    GenParticle * c =(GenParticle*) genParticles_->At(k);
    if (c==0) continue;
    TString out;
    out.Form("%d  \t %d \t %d \t %d \t %.1f \t %d",k,c->PID,c->M1,c->D1,c->Mass,c->Status);
    cout<<out<<endl;
  }

}

double McTruthInfo::getNloKFactor(int code) {

  //need to return kfactors_[code] but more efficiently
  map<int,double>::iterator it = kfactors_.find(code);
  if ( it != kfactors_.end()) return it->second;
  return 1;
}


void McTruthInfo::setNloKFactor_fromfileName(TString fname){

  if(fname.Contains("naturalModel1")) setNloKFactor_NM1();
  else if(fname.Contains("naturalModel2")) setNloKFactor_NM2();
  else if(fname.Contains("naturalModel3")) setNloKFactor_NM3();
  else if(fname.Contains("stc")) setNloKFactor_STC();
  else if(fname.Contains("stoc")) setNloKFactor_STOC();

}



void McTruthInfo::setNloKFactor_NM1(){

  cout<<" Setting k-factors for NM1"<<endl;

  // the actual vaules should match with https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSYSignalSamplesForFutureStudies
  setNloKFactor(2000000,1.17);//  2000000 == slepton pair production
  setNloKFactor(2000001,1.19);//  2000001 == slepton sneutrino production
  setNloKFactor(2000002,1.17);//  2000002 == sneutrino pair production
    // split the ewkino process
  setNloKFactor(200002,1.20); //  200002  == N1N1 production
  setNloKFactor(200004,1.20); //  200004  == C1N1 production
  setNloKFactor(200006,1.19); //  200004  == C1C1 production
  setNloKFactor(200010,1.1); //  200010  == N2N1 production
  setNloKFactor(200012,1.20); //  200012  == N2C1 production
  //from here only default values will be updated
  setNloKFactor(200018,1.20); //  200018  == N2N2 production
  setNloKFactor(200028,1.20); //  200028  == N3N1 production
  setNloKFactor(200030,1.20); //  200030  == N3C1 production
  setNloKFactor(200036,1.20); //  200036  == N3N2 production
  setNloKFactor(200054,1.20); //  200054  == N3N3 production
  setNloKFactor(200082,1.20); //  200082  == N4N1 production
  setNloKFactor(200084,1.20); //  200084  == N4C1 production
  setNloKFactor(200090,1.20); //  200090  == N4N2 production
  setNloKFactor(200108,1.20); //  200108  == N4N3 production
  setNloKFactor(200162,1.20); //  200162  == N4N4 production
  setNloKFactor(200244,1.20); //  200244  == C2N1 production
  setNloKFactor(200246,1.20); //  200246  == C2C1 production
  setNloKFactor(200252,1.20); //  200252  == C2N2 production
  setNloKFactor(200270,1.20); //  200270  == C2N3 production
  setNloKFactor(200324,1.20); //  200324  == C2N4 production
  setNloKFactor(200486,1.20); //  200486  == C2C2 production
  /// end of ewkino block
  setNloKFactor(20000,1.6);  //  20000   == gluino pair production
  setNloKFactor(2000,1.77);   //  2000    == stop pair production
  setNloKFactor(200,1.84);    //  200     == sbottom pair production
  setNloKFactor(10010,2.0);  //  10010   == gluino-squark production
  //setNloKFactor(20,);     //  20      == squark-squark production 
}

void McTruthInfo::setNloKFactor_NM2(){
  cout<<" Setting k-factors for NM2"<<endl;

  // setNloKFactor(2000000,);//  2000000 == slepton pair production
  //setNloKFactor(2000001,);//  2000001 == slepton sneutrino production
  //setNloKFactor(2000002,);//  2000002 == sneutrino pair production
    // split the ewkino process
  setNloKFactor(200002,1.2); //  200002  == N1N1 production
  setNloKFactor(200004,1.23); //  200004  == C1N1 production
  setNloKFactor(200006,1.21); //  200004  == C1C1 production
  setNloKFactor(200010,1.15); //  200010  == N2N1 production
  setNloKFactor(200012,1.21); //  200012  == N2C1 production
  //from here only default values will be updated
  setNloKFactor(200018,1.2); //  200018  == N2N2 production
  setNloKFactor(200028,1.2); //  200028  == N3N1 production
  setNloKFactor(200030,1.2); //  200030  == N3C1 production
  setNloKFactor(200036,1.2); //  200036  == N3N2 production
  setNloKFactor(200054,1.2); //  200054  == N3N3 production
  setNloKFactor(200082,1.2); //  200082  == N4N1 production
  setNloKFactor(200084,1.2); //  200084  == N4C1 production
  setNloKFactor(200090,1.2); //  200090  == N4N2 production
  setNloKFactor(200108,1.2); //  200108  == N4N3 production
  setNloKFactor(200162,1.2); //  200162  == N4N4 production
  setNloKFactor(200244,1.2); //  200244  == C2N1 production
  setNloKFactor(200246,1.2); //  200246  == C2C1 production
  setNloKFactor(200252,1.2); //  200252  == C2N2 production
  setNloKFactor(200270,1.2); //  200270  == C2N3 production
  setNloKFactor(200324,1.2); //  200324  == C2N4 production
  setNloKFactor(200486,1.2); //  200486  == C2C2 production
  /// end of ewkino block
  setNloKFactor(20000,1.6);  //  20000   == gluino pair production
  setNloKFactor(2000,1.77);   //  2000    == stop pair production
  setNloKFactor(200,1.84);    //  200     == sbottom pair production
  setNloKFactor(10010,2.0);  //  10010   == gluino-squark production
  //setNloKFactor(20,);     //  20      == squark-squark production 
}


void McTruthInfo::setNloKFactor_NM3(){
  cout<<" Setting k-factors for NM3"<<endl;

  // setNloKFactor(2000000,);//  2000000 == slepton pair production
  //setNloKFactor(2000001,);//  2000001 == slepton sneutrino production
  //setNloKFactor(2000002,);//  2000002 == sneutrino pair production
    // split the ewkino process
  setNloKFactor(200002,1.2); //  200002  == N1N1 production
  setNloKFactor(200004,1.28); //  200004  == C1N1 production
  setNloKFactor(200006,1.27); //  200004  == C1C1 production
  setNloKFactor(200010,1.28); //  200010  == N2N1 production
  setNloKFactor(200012,1.28); //  200012  == N2C1 production
  //from here only default values will be updated
  setNloKFactor(200018,1.2); //  200018  == N2N2 production
  setNloKFactor(200028,1.2); //  200028  == N3N1 production
  setNloKFactor(200030,1.2); //  200030  == N3C1 production
  setNloKFactor(200036,1.2); //  200036  == N3N2 production
  setNloKFactor(200054,1.2); //  200054  == N3N3 production
  setNloKFactor(200082,1.2); //  200082  == N4N1 production
  setNloKFactor(200084,1.2); //  200084  == N4C1 production
  setNloKFactor(200090,1.2); //  200090  == N4N2 production
  setNloKFactor(200108,1.2); //  200108  == N4N3 production
  setNloKFactor(200162,1.2); //  200162  == N4N4 production
  setNloKFactor(200244,1.2); //  200244  == C2N1 production
  setNloKFactor(200246,1.2); //  200246  == C2C1 production
  setNloKFactor(200252,1.2); //  200252  == C2N2 production
  setNloKFactor(200270,1.2); //  200270  == C2N3 production
  setNloKFactor(200324,1.2); //  200324  == C2N4 production
  setNloKFactor(200486,1.2); //  200486  == C2C2 production
  /// end of ewkino block
  setNloKFactor(20000,1.6);  //  20000   == gluino pair production
  setNloKFactor(2000,1.8);   //  2000    == stop pair production
  setNloKFactor(200,1.84);    //  200     == sbottom pair production
  setNloKFactor(10010,2.0);  //  10010   == gluino-squark production
  //setNloKFactor(20,);     //  20      == squark-squark production 
}


void McTruthInfo::setNloKFactor_STC(){
  cout<<" Setting k-factors for STC"<<endl;

  setNloKFactor(2000000,1.2);//  2000000 == slepton pair production
  setNloKFactor(2000001,1.22);//  2000001 == slepton sneutrino production
  setNloKFactor(2000002,1.21);//  2000002 == sneutrino pair production
    // split the ewkino process
  setNloKFactor(200002,1.2); //  200002  == N1N1 production
  setNloKFactor(200004,1.27); //  200004  == C1N1 production
  setNloKFactor(200006,1.27); //  200004  == C1C1 production
  setNloKFactor(200010,1.21); //  200010  == N2N1 production
  setNloKFactor(200012,1.27); //  200012  == N2C1 production
  //from here only default values will be updated
  setNloKFactor(200018,1.2); //  200018  == N2N2 production
  setNloKFactor(200028,1.2); //  200028  == N3N1 production
  setNloKFactor(200030,1.2); //  200030  == N3C1 production
  setNloKFactor(200036,1.2); //  200036  == N3N2 production
  setNloKFactor(200054,1.2); //  200054  == N3N3 production
  setNloKFactor(200082,1.2); //  200082  == N4N1 production
  setNloKFactor(200084,1.2); //  200084  == N4C1 production
  setNloKFactor(200090,1.2); //  200090  == N4N2 production
  setNloKFactor(200108,1.2); //  200108  == N4N3 production
  setNloKFactor(200162,1.2); //  200162  == N4N4 production
  setNloKFactor(200244,1.2); //  200244  == C2N1 production
  setNloKFactor(200246,1.2); //  200246  == C2C1 production
  setNloKFactor(200252,1.2); //  200252  == C2N2 production
  setNloKFactor(200270,1.2); //  200270  == C2N3 production
  setNloKFactor(200324,1.2); //  200324  == C2N4 production
  setNloKFactor(200486,1.2); //  200486  == C2C2 production
  /// end of ewkino block
  setNloKFactor(20000,1.3);  //  20000   == gluino pair production
  setNloKFactor(2000,1.7);   //  2000    == stop pair production
  setNloKFactor(200,1.75);    //  200     == sbottom pair production
  setNloKFactor(10010,2.24);  //  10010   == gluino-squark production
  //setNloKFactor(20,);     //  20      == squark-squark production 
}


void McTruthInfo::setNloKFactor_STOC(){
  cout<<" Setting k-factors for STOC"<<endl;

  //setNloKFactor(2000000,);//  2000000 == slepton pair production
  //setNloKFactor(2000001,);//  2000001 == slepton sneutrino production
  //setNloKFactor(2000002,);//  2000002 == sneutrino pair production
    // split the ewkino process
  setNloKFactor(200002,1.2); //  200002  == N1N1 production
  setNloKFactor(200004,1.22); //  200004  == C1N1 production
  setNloKFactor(200006,1.17); //  200004  == C1C1 production
  setNloKFactor(200010,1.16); //  200010  == N2N1 production
  setNloKFactor(200012,1.15); //  200012  == N2C1 production
  //from here only default values will be updated
  setNloKFactor(200018,1.2); //  200018  == N2N2 production
  setNloKFactor(200028,1.2); //  200028  == N3N1 production
  setNloKFactor(200030,1.2); //  200030  == N3C1 production
  setNloKFactor(200036,1.2); //  200036  == N3N2 production
  setNloKFactor(200054,1.2); //  200054  == N3N3 production
  setNloKFactor(200082,1.2); //  200082  == N4N1 production
  setNloKFactor(200084,1.2); //  200084  == N4C1 production
  setNloKFactor(200090,1.2); //  200090  == N4N2 production
  setNloKFactor(200108,1.2); //  200108  == N4N3 production
  setNloKFactor(200162,1.2); //  200162  == N4N4 production
  setNloKFactor(200244,1.2); //  200244  == C2N1 production
  setNloKFactor(200246,1.2); //  200246  == C2C1 production
  setNloKFactor(200252,1.2); //  200252  == C2N2 production
  setNloKFactor(200270,1.2); //  200270  == C2N3 production
  setNloKFactor(200324,1.2); //  200324  == C2N4 production
  setNloKFactor(200486,1.2); //  200486  == C2C2 production
  /// end of ewkino block
  setNloKFactor(20000,1.5);  //  20000   == gluino pair production
  setNloKFactor(2000,1.55);   //  2000    == stop pair production
  setNloKFactor(200,2.43);    //  200     == sbottom pair production
  setNloKFactor(10010,2.64);  //  10010   == gluino-squark production
  //setNloKFactor(20,);     //  20      == squark-squark production 
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

  for (int k = 0 ; k< nGenParticles_; k++) {
    GenParticle * c =(GenParticle*) genParticles_->At(k);
    if (c==0) continue; //try to prevent crashes....
    if (c->Status!=3) continue; //status 3 [pythia 6] particles only!
    int pid = std::abs(c->PID);

    if (pid==pidToFind) { 
      //does this particle come from SUSY?
      if (  checkMom(k,1000000) ) np++;
    }
  }
  return np;

}

int McTruthInfo::countTrueLeptons(leptonType lt) {

  //filling this lepton list was added after the fact
  //would actually be more efficient to fill the lepton list once per event
  //and use that list to do things like counting each lepton type
  //but not important enough to worry about
  leptons_.clear();
  neutrinos_.clear();

  int n=0;
  //easy as pie, except that depending on what is stored, we might actually
  //want to restrict ourselves to status 3 or the Pythia 8 equivalent
  for (int k = 0 ; k< nGenParticles_; k++) {
    GenParticle * c =(GenParticle*) genParticles_->At(k);
    if (c==0) continue; //try to prevent crashes....
    if (c->Status!=3) continue; //status 3 [pythia 6] particles only!
    int pid = std::abs(c->PID);

    if (pid==11 && lt==kElMu) n++;
    else if (pid==13 && lt==kElMu) n++;
    else if (pid==15 && lt==kTau) n++;

    if (pid==11||pid==13||pid==15) leptons_.push_back(c);
    else if (pid==12||pid==14||pid==16) neutrinos_.push_back(c);
  }

  return n;
}

std::pair<float,float> McTruthInfo::getGenMet() { //return the |vector sum pt| of the status==3 neutrinos in the event (and phi)
  if (neutrinos_.size()==0) return make_pair(0,-99);
  
  TLorentzVector metvec; //init to 0 by default

  for (size_t k=0;k<neutrinos_.size();++k) {
    metvec += neutrinos_.at(k)->P4();
  }

  pair<float,float>  mymet = make_pair( metvec.Vect().Perp() , -(metvec.Vect().Phi() ) );
  return mymet;
}
/*
look for chi_20 -> slepton lepton -> lepton chi_10 lepton
code is 10*nStaus + nSElectron+Smuon
should allow for a simple cut that removes events with taus or
multiple edge decays while still retaining that info in case it is of interest

name of the function is now a misnomer as it works for chi2 or chi4
 */
int McTruthInfo::findChi2ToChi1(int chiToFind) {

  vector<GenParticle*> * my_edge_lep1;
  vector<GenParticle*> * my_edge_lep2;

  if (chiToFind==2) {
    my_edge_lep1 = &edge_l1_;
    my_edge_lep2 = &edge_l2_;
  }
  else if (chiToFind==4) {
    my_edge_lep1 = &edge4_l1_;
    my_edge_lep2 = &edge4_l2_;
  }
  else assert(0);

  my_edge_lep1->clear();
  my_edge_lep2->clear();

  int code = 0;

  int grandParentPid = 1000023;
  if (chiToFind==4) grandParentPid=1000035;

  GenParticle * edge_l1 = 0;

  for (int k = 0 ; k< nGenParticles_; k++) {
    GenParticle * c =(GenParticle*) genParticles_->At(k);
    if (c==0) continue; //try to prevent crashes....
    int pid = std::abs(c->PID);

    if (pid==1000022) { //chi10
      //is the mom a slepton?
      bool el = checkMom(k,1000011,0);
      bool mu=checkMom(k,1000013,0);
      bool tau = checkMom(k,1000015,0);
      if (el||mu||tau ) {
	//is the grandmom a chi_20?
	if (checkMom(k,grandParentPid,1) ) {
	  if (el) ++code;
	  else if (mu) ++code;
	  else if (tau) code += 10;
	  else assert(0);
	}
      }
    }
    //not needed for primary logic of finding decay chain, but i also want to find leptons themselves
    else if (pid==11 || pid==13) {
      //look for lepton from slepton, with gmom of chi20
      if ( checkMom(k,pid+1000000,0) && checkMom(k,grandParentPid,1) ) {
	my_edge_lep2->push_back( c);
	assert(edge_l1);
	my_edge_lep1->push_back( edge_l1);
      }
      //look for lepton from chi20
      else if (checkMom(k,grandParentPid,0)) edge_l1 = c;
	// edge_l1_.push_back( c);
      //vector structure only needed for the case where there are 2 N2->N1 decays
    }

  }

  //  if (chiToFind==2 && code!=0) cout<<"[code] "<<chiToFind<<" "<<code<<"\t"<<edge_l1_.size()<<" "<<edge_l2_.size()<<endl;
  //if (chiToFind==4 && code!=0) cout<<"[code] "<<chiToFind<<" "<<code<<"\t"<<edge4_l1_.size()<<" "<<edge4_l2_.size()<<endl;

  return code;
}

float McTruthInfo::getGenMll(int index) {

  if (index==1 && edge_l1_.size()>=1) {
    TLorentzVector mll = edge_l1_[0]->P4() + edge_l2_[0]->P4();
    return mll.M();
  }
  else if (index==2 && edge_l1_.size()>=2) {
    TLorentzVector mll = edge_l1_[1]->P4() + edge_l2_[1]->P4();
    return mll.M();
  }
  return -1;
}

float McTruthInfo::getGenEdgeLeptonPt(int index) { //forget about the possibility of a second edge decay

  if ( edge_l1_.size()>=1) {
    if (index==1) return edge_l1_[0]->P4().Pt();
    if (index==2) return edge_l2_[0]->P4().Pt();
  }

  return -1;
}

float McTruthInfo::getIsolationOfMatch(const unsigned int ilep,const TClonesArray* els,const TClonesArray* mus) {
  const float dummy=-1e9;
  if (ilep >= leptons_.size()) return dummy;

  GenParticle * g = leptons_.at(ilep);
  TLorentzVector genP4=g->P4();

  int  flavor = std::abs(g->PID);
  int nr=0;
  if ( flavor==11) nr=els->GetEntries();
  else if (flavor==13) nr=mus->GetEntries();
  else     return dummy; //tau
  
  for (int ireco=0; ireco<nr;ireco++) {
    TLorentzVector * recoP4=0;
    if (flavor==11) recoP4 = new TLorentzVector( ((Electron*)els->At(ireco))->P4());
    else if (flavor==13) recoP4 = new TLorentzVector( ((Muon*)mus->At(ireco))->P4());
    //    double dr = Util::GetDeltaR(genP4.Eta(), recoP4->Eta(), genP4.Phi(), recoP4->Phi());
    double dr = recoP4->DeltaR(genP4);



    if (dr<0.03 && std::abs(1- genP4.Pt()/recoP4->Pt() )<0.3) {
      delete recoP4;
      float iso = (flavor==11) ?  ((Electron*)els->At(ireco))->IsolationVar :  ((Muon*)mus->At(ireco))->IsolationVar;
      return iso;
    }
    delete recoP4; 
  }
  return dummy; //did not find a match
}

bool McTruthInfo::matchesChi2ToChi1Gen(const TLorentzVector & l1, const TLorentzVector & l2,int l1_flavor,int l2_flavor,int chiToMatch) {
  //check if l1 and l2 are DR matches to the gen-level edge leptons
  vector<GenParticle*> * edge_l1;
  vector<GenParticle*> * edge_l2;
  if (chiToMatch==2) {
    edge_l1 = &edge_l1_;
    edge_l2 = &edge_l2_;
  }
  else if (chiToMatch==4) {
    edge_l1 = &edge4_l1_;
    edge_l2 = &edge4_l2_;
  }
  else assert(0);

  const double max = 0.05;

  bool match=false;

  if (!(edge_l2->size() == edge_l1->size())) {
    cout<<"bad! "<<edge_l2->size()<<" "<<edge_l1->size()<<endl;
    assert(0);
  }
  //cout<<"~~ "<<edge_l2->size()<<endl;
  for (size_t j = 0; j<edge_l2->size() ; j++) {

    //check flavor
    //    cout<<"flavors "<<l1_flavor<<" "<<edge_l1_[j]->PID<<endl;
    if (std::abs(l1_flavor)!= std::abs( (*edge_l1)[j]->PID)) continue;

    TLorentzVector g_l1 = (*edge_l1)[j]->P4();
    TLorentzVector g_l2 = (*edge_l2)[j]->P4();
    //associate via charge
    if (l1_flavor == - (*edge_l2)[j]->PID) { //tricky! in l1_flavor I used e- is negative. PID uses e- is positive
      g_l1 = (*edge_l2)[j]->P4();
      g_l2 = (*edge_l1)[j]->P4();
    }

    //    l1.Print();
    //    g_l1.Print();
    //    l2.Print();
    //    g_l2.Print();

    //    double dr1 =  Util::GetDeltaR(l1.Eta(),g_l1.Eta(),l1.Phi(), g_l1.Phi());
    double dr1 = l1.DeltaR(g_l1);  

    //    double dr2 =  Util::GetDeltaR(l2.Eta(),g_l2.Eta(),l2.Phi(), g_l2.Phi());
    double dr2 =  l2.DeltaR(g_l2);




    if (dr1 < max && dr2<max) {
      //      cout<<"Is match"<<endl;
      match=true;
      break;
    }
    //    else cout<<"No match"<<endl;
  }
  return match;
}

int McTruthInfo::getSusyProductionProcess() {

  vector<int> susyMoms = findSusyMoms();

  int n_slepton = 0;
  // distinguish snu and slep
  int n_sneutrino = 0;

  int n_ewkino = 0;
  // in order to distinguish ewkino
  int n_ewkino_c1 = 0;
  int n_ewkino_c2 = 0;
  int n_ewkino_n1 = 0;
  int n_ewkino_n2 = 0;
  int n_ewkino_n3 = 0;
  int n_ewkino_n4 = 0;

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
    else if ( pid >= 1000011 && pid <= 1000016){
      n_slepton++;   
      if(pid%2==0) n_sneutrino++;
	}
    else if ( pid >= 2000011 && pid <= 2000016){
      n_slepton++;
      if(pid%2==0) n_sneutrino++;
	}
    else if ( pid >= 2000001 && pid <= 2000004) n_squark++;
    else if ( pid==1000021) n_gluino++;
    else if ( pid>=1000012 && pid <= 1000037) {
      // distinguish individual ewkinos 
      n_ewkino++;
      if (pid==1000022) n_ewkino_n1++;
      else if (pid==1000023) n_ewkino_n2++;
      else if (pid==1000024) n_ewkino_c1++;
      else if (pid==1000025) n_ewkino_n3++;
      else if (pid==1000035) n_ewkino_n4++;
      else if (pid==1000037) n_ewkino_c2++;
    }
    else n_other++;
  }
  return   n_slepton*1000000 + n_sneutrino
    +    n_ewkino*100000 + (n_ewkino_n1 + 3*(n_ewkino_n2 + 3*(n_ewkino_c1 + 3*(n_ewkino_n3 + 3*(n_ewkino_n4 + 3*(n_ewkino_c2))))))
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


  for (int k = 0 ; k< nGenParticles_; k++) {
    GenParticle * c =(GenParticle*) genParticles_->At(k);
    if (c==0) continue; //try to prevent crashes....
    int pid = std::abs(c->PID);
    int momIndex1 = c->M1;
    if (momIndex1<0 || momIndex1>=nGenParticles_) continue;
    GenParticle * theMom = (GenParticle*)genParticles_->At(momIndex1);
    if (theMom==0) continue; //saw momIndex1 be invalid once
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
  int momIndex1 = theCand->M1;// -3;// test OFFSET to compensate for a bug in sample. doesn't work anyway
  if (momIndex1 >=0 && momIndex1< nGenParticles_ ) { //sanity check on momIndex
    GenParticle * theMom = (GenParticle*)genParticles_->At(momIndex1);
    if (theMom==0)       return false; //TClonesArray::At returns 0 if index is not found
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

int McTruthInfo::GetTtbarDecayCode(float & genmll) {
  //classify ttbar events using my old coding scheme (now modified a bit - see below)

  std::vector<GenParticle*> leptons; //was tempted to use the already defined global vector edge_l1_ for this
  //but that seems sloppy

  int ne=0, nmu=0, ntau=0, nq=0,nb=0;
  //  cout<<" ~~~~~~~~~~~"<<endl;
  for (int k = 0 ; k< nGenParticles_; k++) {
    GenParticle * c =(GenParticle*) genParticles_->At(k);
    if (c==0) continue; //try to prevent crashes....
    if (c->Status!=3) continue; //status 3 check
    int pid = std::abs(c->PID);
    //    cout<<c->PID<<" "<<c->M1<<endl;
    if (pid == 11 || pid == 13 || pid==15 || pid==1 ||pid==2 ||pid==3||pid==4 ||pid==5) { //find ultimate top->W decay product
      bool      topIsUltimateMom = checkMom(k, 6);

      if (topIsUltimateMom) {
	//nb - using the edge_l1_
	if      (pid ==11) {ne++;  leptons.push_back(c);}
	else if (pid ==13) {nmu++; leptons.push_back(c);}
	else if (pid ==15) ntau++;
	else if (pid==5) nb++;
	else  nq++;
      }
    }
  }

  nq/=2;

  genmll=-1;
  if (leptons.size()==2) {
    TLorentzVector a=leptons[0]->P4();
    TLorentzVector b=leptons[1]->P4();
    genmll= (a+b).M();
  }

  //updated coding scheme (a small perturbation on the old one)
  /*
0  = uncategorized
11 = dilepton (ee)
12 = dilepton (emu)
13 = dilepton (mumu)
2  = all-hadronic
3  = semi-lep e
4  = semi-lep mu
5  = semi-lep tau
8  = dilepton (tau tau)
9  = dilepton (tau+e/mu)
  */

  if (ne==2) return 11;
  else if (nmu==2) return 13;
  else if (ne==1 && nmu==1) return 12;
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
  for (int k = 0 ; k< nGenParticles_; k++) {
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

      //      double dr = Util::GetDeltaR(genP4.Eta(), recoP4->Eta(), genP4.Phi(), recoP4->Phi());
      double dr = recoP4->DeltaR(genP4);

      //it appears that there is ~no resolution smearing in Delphes for
      //phi and eta, so dr match is sufficient, and can be very tight
      //leave loose pt match, but not really necessary
      //(it appears that pT tends to be well within 5%, so this is very loose)
      if (dr<0.03 && std::abs(1- genP4.Pt()/recoP4->Pt() )<0.3) {
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

//  LocalWords:  setNloKFactor
