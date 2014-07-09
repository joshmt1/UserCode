#ifndef JETLEPTONCLEANER_H
#define JETLEPTONCLEANER_H

/*
This class contains 'database's of sample info.
Mainly this is the cross-section for each sample.
Also it classifies each sample into a process enum

 */

#include "TClonesArray.h"

class JetLeptonCleaner {

public:

  JetLeptonCleaner(TClonesArray* jets,TClonesArray* els,TClonesArray* mus,TClonesArray* photons);
  ~JetLeptonCleaner();

private:
  template <class T> void applyIsolation(TClonesArray* collection) ;
  template <class T> void cleanJets(TClonesArray* jets,TClonesArray* leps);
  float isolationCut_;
  float drmin_;

};

#endif
