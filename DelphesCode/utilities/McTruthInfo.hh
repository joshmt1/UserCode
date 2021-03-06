// -*- C++ -*-
//container for flat root tree

#ifndef MCTRUTHINFO_H
#define MCTRUTHINFO_H

#include "TClonesArray.h"
#include "TLorentzVector.h"
#include <vector>
#include <utility>
#include "classes/DelphesClasses.h"
#include <map>

class McTruthInfo {

public:
  enum leptonType {kElMu,kTau};

  McTruthInfo();
  ~McTruthInfo();
  int GetTtbarDecayCode(float & genmll); //for ttbar

  //for SUSY
  std::vector<int> findSusyMoms();
  int getSusyProductionProcess();
  void setNloKFactor(int code,double k) {kfactors_[code]=k;}
  void setNloKFactor_fromfileName(TString fname);
  void setNloKFactor_NM1();
  void setNloKFactor_NM2();
  void setNloKFactor_NM3();
  void setNloKFactor_STC();
  void setNloKFactor_STOC();
  double getNloKFactor(int code);
  int findChi2ToChi1(int chiToFind=2); //return 10*nStaus + nSElectron+Smuon
  float getGenMll(int index); //index is which edge decay (in rare case of 2)
  float getGenEdgeLeptonPt(int index); //index is which lepton; only supports first edge decay in the event 
  bool matchesChi2ToChi1Gen(const TLorentzVector & l1, const TLorentzVector & l2,int l1_flavor,int l2_flavor,int chiToMatch=2);
  int findPinSusy(int pidToFind);

  //for anything
  void Set(TClonesArray* g) {genParticles_=g; nGenParticles_=g->GetEntries();}
  int countTrueLeptons(leptonType lt);
  std::pair<float,float> getGenMet(); //return the |vector sum pt| of the status==3 neutrinos in the event (and phi)
  int countNeutrinos() {return (int)neutrinos_.size();}
  std::vector<GenParticle*> getGenLeptons() {return leptons_;}
  float getIsolationOfMatch(const unsigned int ilep,const TClonesArray* els,const TClonesArray* mus);
  void Dump();

  //DY MC truth
  std::vector< std::pair< TLorentzVector, int> > GetDYTruth() ;
  std::vector<int> MatchDYRecoGen(const std::vector< std::pair< TLorentzVector, int> > & genDY,
				  TClonesArray * els,TClonesArray * mus);


private:
  TClonesArray * genParticles_;
  int nGenParticles_;
  std::vector<GenParticle*> edge_l1_; //leptons from chi2->chi1
  std::vector<GenParticle*> edge_l2_;
  std::vector<GenParticle*> edge4_l1_;//leptons from chi4->chi1
  std::vector<GenParticle*> edge4_l2_;
  std::vector<GenParticle*> leptons_;
  std::vector<GenParticle*> neutrinos_;
  std::map<int,double> kfactors_;  
  bool checkMom(int index, int PidToLookFor,int recursionsLeft=999);

  bool isSusy(int pid) ;


};

#endif
