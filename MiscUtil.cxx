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

  //======== container for run, lumisection, event number =========
  class eventID {
  public:
    eventID();    
    eventID(ULong64_t RunNumber, ULong64_t LumiSection, ULong64_t EventNumber);
    ~eventID();
    
    bool operator< (const eventID & id) const;
    bool operator== (const eventID & id) const;
    bool operator!= (const eventID & id) const;
    
    ULong64_t run;
    ULong64_t ls;
    ULong64_t ev;
    
  };
  
  eventID::eventID() : run(0),ls(0),ev(0) {}
  eventID::eventID(ULong64_t RunNumber, ULong64_t LumiSection, ULong64_t EventNumber) : 
    run(RunNumber),ls(LumiSection),ev(EventNumber) {}
  eventID::~eventID() {}
  
  bool eventID::operator== (const eventID & id) const {

    if (ev == id.ev &&
	ls == id.ls &&
	run == id.run) return true;
    
    return false;
  }

  bool eventID::operator!= (const eventID & id) const {
    return !(*this == id);
  }

  bool eventID::operator< (const eventID & id) const {
    
    if (run < id.run) return true;
    else if (run > id.run) return false;    
    else { //if run is equal

      //now compare lumi section
      if ( ls < id.ls ) return true;
      else if (ls > id.ls) return false;    
      else { //if ls is equal
	
	if ( ev < id.ev ) return true;
	else if (ev > id.ev) return false;    
	
	else { //these are the same event!
	  return false;
	}
      }
    }

    //this shouldn't happen
    assert(0);
    return false;
  }
  
} //end of namespace


