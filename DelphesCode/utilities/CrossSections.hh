#ifndef CROSSSECTIONS_H
#define CROSSSECTIONS_H

#include "TString.h"

class CrossSections {

public:
  CrossSections(const TString samplename);
  ~CrossSections();

  double Get() {return xs_;}

private:
  //  TString samplename_;
  double xs_;

};

#endif
