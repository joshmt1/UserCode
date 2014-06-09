// -*- C++ -*-
//compute m(ll) for Edge search

#ifndef MLLCOMPUTER_H
#define MLLCOMPUTER_H

#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

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
  int GetLeptonFlavor(int index);
  float GetLeptonIsolation(int index);//absolute isolation!
  int GetNExtraLeptons() {return extra_leptons_.size();}
  TLorentzVector GetExtraLeptonP4(int index);

  float GetMee_Test() ;//special test
  float GetMll_random() ;

  //cuts -- making these public out of laziness
  float minpt_;
  float maxetacut_;//cut value
  bool removegap_;
  //algorithm customizations
  bool randomizeLeptons_;
  int seed_;

private:
  TClonesArray* el_;
  TClonesArray* mu_;

  TRandom3 * rand_;

  bool foundGood_;
  //  std::vector<int> good_el_;
  //  std::vector<int> good_mu_;
  void findGoodLeptons() ;
  float maxEta_;//largest value found in event
  bool isSF_;
  TLorentzVector* p4_l1_;
  TLorentzVector* p4_l2_;
  float l1_iso_;
  float l2_iso_;
  int l1_flavor_;
  int l2_flavor_;

  // pair(pT, pair(flavor,index))
  std::set< std::pair<float, std::pair<LeptonFlavor,int> > > lep_pt_;
  std::set< std::pair<float, std::pair<LeptonFlavor,int> > > extra_leptons_;
  int GetCharge( std::set< std::pair<float, std::pair<LeptonFlavor,int> > >::reverse_iterator & rit);
  float GetIsolation( std::set< std::pair<float, std::pair<LeptonFlavor,int> > >::reverse_iterator & rit);
  TLorentzVector Get4Vector( std::set< std::pair<float, std::pair<LeptonFlavor,int> > >::reverse_iterator & rit) ;

};

#endif
