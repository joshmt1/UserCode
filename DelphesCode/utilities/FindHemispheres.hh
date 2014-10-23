// -*- C++ -*-

/*
implement the MT2 hemisphere finding algorithm

This implementation by Josh Thompson (Cornell) based on the description in CMS AN-2013/215 (ETH Zurich)
Gives mostly identical results to ETH code, but is slower
*/

#ifndef FINDHEMISPHERES_H
#define FINDHEMISPHERES_H

#include "TLorentzVector.h"
#include "TClonesArray.h"

#include <set>
#include <vector>

class FindHemispheres {

public:
  FindHemispheres( TClonesArray * jets);
  ~FindHemispheres();

  double Px(int hemi);
  double Py(int hemi);
  std::vector<int> getGrouping();

  int maxIter_;

private:
  void find() ;
  TLorentzVector getAxis( const std::set<int> & jets );
  double getLundDistance( const TLorentzVector & ax, int ijet) ;
  double getDijetInvMass( int i, int j);
  bool isGoodJet(int jj);

  TClonesArray * jets_;
  bool done_;

  //indices of jets in hemispheres 1 and 2
  std::set<int> hemi1_;
  std::set<int> hemi2_;
  //pseudojets corresponding to the hemispheres
  TLorentzVector pjet1_;
  TLorentzVector pjet2_;

};


#endif

