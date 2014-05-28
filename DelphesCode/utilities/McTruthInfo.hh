// -*- C++ -*-
//container for flat root tree

#ifndef MCTRUTHINFO_H
#define MCTRUTHINFO_H

#include "TClonesArray.h"
#include "TLorentzVector.h"
#include <vector>
#include <utility>
#include "classes/DelphesClasses.h"

class McTruthInfo {

public:
  enum leptonType {kElMu,kTau};

  McTruthInfo();
  ~McTruthInfo();
  int GetTtbarDecayCode(); //for ttbar

  //for SUSY
  std::vector<int> findSusyMoms();
  int getSusyProductionProcess();
  int findChi2ToChi1(); //return 10*nStaus + nSElectron+Smuon
  bool matchesChi2ToChi1Gen(const TLorentzVector & l1, const TLorentzVector & l2,int l1_flavor,int l2_flavor);
  int findPinSusy(int pidToFind);

  //for anything
  void Set(TClonesArray* g) {genParticles_=g;}
  int countTrueLeptons(leptonType lt);
  void Dump();

  //DY MC truth
  std::vector< std::pair< TLorentzVector, int> > GetDYTruth() ;
  std::vector<int> MatchDYRecoGen(const std::vector< std::pair< TLorentzVector, int> > & genDY,
				  TClonesArray * els,TClonesArray * mus);


private:
  TClonesArray * genParticles_;
  std::vector<GenParticle*> edge_l1_;
  std::vector<GenParticle*> edge_l2_;
  bool checkMom(int index, int PidToLookFor,int recursionsLeft=999);

  bool isSusy(int pid) ;


};

#endif
