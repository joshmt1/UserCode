#ifndef CROSSSECTIONS_H
#define CROSSSECTIONS_H

/*
This class contains 'database's of sample info.
Mainly this is the cross-section for each sample.
Also it classifies each sample into a process enum

 */

#include "TString.h"

class CrossSections {

public:
  enum ProcessType { kTop, kBoson, kRare, kSignal, kNone };

  CrossSections(const TString samplename);
  ~CrossSections();

  double Get() {return xs_;}
  ProcessType GetProcess() {return proc_;}

private:
  //  TString samplename_;
  void SetProc(TString name);

  double xs_;
  ProcessType proc_;

};

#endif
