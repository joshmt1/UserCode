/*
a place for very simple functions that don't depend on anything else

can be easily #included into other code
*/

#include "TString.h"
#include <iostream>

namespace jmt {

  //======misc utilities======
  //gets rid of = > < from cuts in order to be better included in file names
  TString fortranize(TString cut) {

    cut.ReplaceAll("==","eq");
    cut.ReplaceAll(">=","gte");
    cut.ReplaceAll("<=","lte");
    
    cut.ReplaceAll(">","gt");
    cut.ReplaceAll("<","lt");

    cut.ReplaceAll("/","Over");
    cut.ReplaceAll("+","Plus");
    cut.ReplaceAll("*","Times");

    cut.ReplaceAll("*","Times");
    
    //this is pretty ugly
    cut.ReplaceAll("(","L");
    cut.ReplaceAll(")","R");

    return cut;
  }

  void weightedMean(Int_t n, Double_t *a, Double_t *e) {

    Double_t numerator=0;
    Double_t denominator=0;

    for (Int_t i=0; i<n; ++i) {
      Double_t esq = e[i]*e[i];
      if (esq == 0) {cout<<"Error is zero!"<<endl; return;}
      numerator   += a[i] / esq;
      denominator += 1.0  / esq;
    }

    cout<<"mean = "<<numerator / denominator<<" +/- "<<sqrt(1/denominator)<<endl;

  }

  //implemented in TMath in later versions of root, so i'm using the TMath name
  bool AreEqualAbs( const double & a, const double & b, const double epsilon=0.001) {

    return ( fabs(a-b) < epsilon ) ;

  }

  // == error propagation ==

  //because not all versions of ROOT have it built in
  double errOnIntegral(const TH1D* h, int binlow=1, int binhigh=-1) {
    
    if (binhigh == -1) binhigh = h->GetNbinsX();
    
    double err=0;
    for (int i = binlow; i<= binhigh; i++) {
      err += h->GetBinError(i) * h->GetBinError(i);
    }
    
    if (err<0) return err;
    return sqrt(err);
  }

  //return error on a/b
  double errAoverB(double a, double aerr, double b, double berr) {
    if (b==0 || a==0) {
      cout<<"Warning in errAoverB -- a or b is zero!"<<endl;
      return -1;
    }
    return (a/b)*sqrt( pow( aerr/a, 2) + pow( berr/b, 2) );
  }

  //return error on a*b
  double errAtimesB(double a, double aerr, double b, double berr) {
    if (b==0 || a==0) {
      cout<<"Warning in errAtimesB -- a or b is zero!"<<endl;
      return -1;
    }
    
    return a*b*sqrt( pow( aerr/a, 2) + pow( berr/b, 2) );
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


