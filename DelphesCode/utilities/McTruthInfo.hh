// -*- C++ -*-
//container for flat root tree

#ifndef MCTRUTHINFO_H
#define MCTRUTHINFO_H

#include "TClonesArray.h"

class McTruthInfo {

public:
  McTruthInfo();
  ~McTruthInfo();
  int GetTtbarDecayCode(TClonesArray* genParticles=0) ;

  void Dump(TClonesArray* genParticles=0);

private:
  TClonesArray * genParticles_;
  bool checkMom(int index, int PidToLookFor);

};

#endif
