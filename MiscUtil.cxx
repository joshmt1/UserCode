/*
a place for very simple functions that don't depend on anything else

can be easily #included into other code
*/

#include "TString.h"

namespace jmt {

  //======misc utilities======
  //gets rid of = > < from cuts in order to be better included in file names
  TString fortranize(TString cut) {

    cut.ReplaceAll("==","eq");
    cut.ReplaceAll(">=","gte");
    cut.ReplaceAll("<=","lte");
    
    cut.ReplaceAll(">","gt");
    cut.ReplaceAll("<","lt");
    
    return cut;
  }

}
