// -*- C++ -*-
//compute m(ll) for Edge search

#ifndef MLLCOMPUTER_H
#define MLLCOMPUTER_H

#include "TClonesArray.h"
#include "TLorentzVector.h"

#include <set>
#include <utility>

class MllComputer {

  enum LeptonFlavor { kElectron=0,kMuon=1};

public:
  MllComputer(TClonesArray* el, TClonesArray* mu);
  ~MllComputer();
  float GetMll();
  float GetMaxEta() {return maxEta_;}
  bool isSF() {return isSF_;}
  TLorentzVector* GetLeptonP4(int index);

  float GetMee_Test() ;//special test

private:
  TClonesArray* el_;
  TClonesArray* mu_;

  bool foundGood_;
  //  std::vector<int> good_el_;
  //  std::vector<int> good_mu_;
  void findGoodLeptons() ;
  float maxEta_;
  bool isSF_;
  TLorentzVector* p4_l1_;
  TLorentzVector* p4_l2_;

  // pair(pT, pair(flavor,index))
  std::set< std::pair<float, std::pair<LeptonFlavor,int> > > lep_pt_;
  int GetCharge( std::set< std::pair<float, std::pair<LeptonFlavor,int> > >::reverse_iterator & rit);
  TLorentzVector Get4Vector( std::set< std::pair<float, std::pair<LeptonFlavor,int> > >::reverse_iterator & rit) ;

};

#endif
