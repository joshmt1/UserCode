// -*- C++ -*-
//container for flat root tree

#ifndef MCTRUTHINFO_H
#define MCTRUTHINFO_H

#include "TClonesArray.h"

class McTruthInfo {

public:
  McTruthInfo();
  ~McTruthInfo();
  int GetTtbarDecayCode();
  std::vector<int> findSusyMoms();
  int getSusyProductionProcess();
  void Dump();
  int findChi2ToChi1(); //return 10*nStaus + nSElectron+Smuon
  int findZinSusy();



  void Set(TClonesArray* g) {genParticles_=g;}

private:
  TClonesArray * genParticles_;
  bool checkMom(int index, int PidToLookFor,int recursionsLeft=999);

  bool isSusy(int pid) ;


};

#endif
