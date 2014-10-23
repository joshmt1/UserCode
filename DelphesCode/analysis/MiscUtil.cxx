/*
a place for very simple functions that don't depend on anything else

can be easily #included into other code
*/

#include "TH1D.h"
#include "TString.h"
#include <iostream>
#include "TMath.h"

#include <cassert>

//for test
#include "RooStats/NumberCountingUtils.h"

namespace jmt {

  //======misc utilities======
  //gets rid of = > < from cuts in order to be better included in file names
  TString fortranize(TString cut,const TString options="") {

    cut.ReplaceAll(" ",""); //remove all spaces

    cut.ReplaceAll("==","eq");
    cut.ReplaceAll(">=","gte");
    cut.ReplaceAll("<=","lte");
    
    cut.ReplaceAll(">","gt");
    cut.ReplaceAll("<","lt");

    cut.ReplaceAll("&&","and");
    cut.ReplaceAll("||","or");

    cut.ReplaceAll("/","Over");
    cut.ReplaceAll("+","Plus");
    cut.ReplaceAll("*","Times");

    //this is pretty ugly
    cut.ReplaceAll("(","L");
    cut.ReplaceAll(")","R");

    //can't decide if this one should be a default or only enabled by an option
    cut.ReplaceAll(".","p");

    if (options.Contains("dash")) cut.ReplaceAll("-","_");

    return cut;
  }

  //utility function for making output more readable
  TString format_nevents(double n,double e, const bool moreDigits=false,const bool latexMode_=true) {
    const TString pm = latexMode_ ?  "\\pm ": " +/- ";

    const int eCutoff = moreDigits ? 10 : 1;
    const int extraDigits = moreDigits ? 1:0;

    TString mathmode = latexMode_ ? "$" : "";
  
    char out[100];
    if (e >= eCutoff || e < 0.00001) { //show whole numbers only
      sprintf(out,"%s%.0f%s%.0f%s",mathmode.Data(),n,pm.Data(),e,mathmode.Data());
    }
    else {
      int nfig = ceil(fabs(log10(e))) + extraDigits;
      TString form="%s%.";
      form+=nfig; form+="f%s%.";
      form+=nfig; form+="f%s";
      sprintf(out,form.Data(),mathmode.Data(),n,pm.Data(),e,mathmode.Data());
    }
    return TString(out);
  }


  Double_t zbi(double n_on, double mu_b_hat, double sigma_b) {
    //implementation of the Zbi significance, courtesy of Dominick Olivito
    
    //double n_on     = 140.;                         // total events in signal region (S+B)
    //double mu_b_hat = 83.33;                        // mean num of BG events expected in sig. region
    //double sigma_b  = 8.333;                        // uncertainty of mu_b_hat  

    //turns out this is already implemented in RooStats
    //works identically except doesn't suffer from precision problem for high significances
    return RooStats::NumberCountingUtils::BinomialObsZ(n_on,mu_b_hat,sigma_b/mu_b_hat);

    /* comment out the implementation from Dominick
    double tau      = mu_b_hat / (sigma_b*sigma_b); // scale factor to corresp. Noff/Non              
    double n_off    = tau*mu_b_hat;
    double P_Bi     = TMath::BetaIncomplete(1./(1.+tau), n_on, n_off+1);
    double Z_Bi=0;
    if (P_Bi >0 && log10(P_Bi)<-15) Z_Bi = 8;
    else                            Z_Bi     = sqrt(2.)*TMath::ErfInverse(1 - 2.*P_Bi);           
   
    if (verbose) {
      cout  <<"  total events in signal region (S+B)               - n_on     " <<n_on      <<endl
	    <<"  mean num of BG events expected in sig. region     - mu_b_hat " <<mu_b_hat  <<endl
	    <<"  uncertainty of mu_b_hat                           - sigma_b  " <<sigma_b   <<endl
	    <<"  scale factor to corresp. Noff/Non                 - tau      " <<tau       <<endl
	    <<"  tau*mu_b_hat                                      - n_off    " <<n_off     <<endl
	    <<"  TMath::BetaIncomplete(1./(1.+tau), n_on, n_off+1) - P_Bi     " <<P_Bi      <<endl
	    <<"  sqrt(2.)*TMath::ErfInverse(1 - 2.*P_Bi)           - Z_Bi     " <<Z_Bi      <<endl;
      
    }

    return Z_Bi;
    */
  }

  Double_t weightedMean(Int_t n, Double_t *a, Double_t *e, const bool returnError=false) {
    using namespace std;

    Double_t numerator=0;
    Double_t denominator=0;

    for (Int_t i=0; i<n; ++i) {
      Double_t esq = e[i]*e[i];
      if (esq == 0) {cout<<"Error is zero!"<<endl; return -1e9;}
      numerator   += a[i] / esq;
      denominator += 1.0  / esq;
    }

    cout<<"mean = "<<numerator / denominator<<" +/- "<<sqrt(1/denominator)<<endl;
    return returnError ? sqrt(1/denominator) : numerator / denominator;
  }

  //implemented in TMath in later versions of root, so i'm using the TMath name
  bool AreEqualAbs( const double & a, const double & b, const double epsilon=0.001) {

    return ( fabs(a-b) < epsilon ) ;

  }

  //deal with the annoyance of logical booleans that are stored in an ntuple as double
  bool doubleToBool( double a) {
    using namespace std;

    int i = TMath::Nint(a);
    if (i==0) return false;
    else if (i==1) return true;
    
    cout<<"[doubleToBool] Something weird going on! "<<a<<"\t"<<i<<endl;
    if (i>1) return true;
    return false;
  }

  //delta phi calculations are pretty commonplace, so put it here
  double deltaPhi(double phi1, double phi2) { 
    double result = phi1 - phi2;
    while (result > TMath::Pi()) result -= 2*TMath::Pi();
    while (result <= -TMath::Pi()) result += 2*TMath::Pi();
    return fabs(result);
  }

  double deltaR2(double eta1, double phi1, double eta2, double phi2) {
    double deta = eta1 - eta2;
    double dphi = deltaPhi(phi1, phi2);
    return deta*deta + dphi*dphi;
  }

  double deltaR(double eta1, double phi1, double eta2, double phi2) {
    return TMath::Sqrt(deltaR2 (eta1, phi1, eta2, phi2));
  }

  int signOf(double a) {
    return (a>= 0.0) ? 1 : -1;
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
    using namespace std;
    if (b==0) {
      cout<<"Warning in errAoverB -- b is zero!"<<endl;
      return -1;
    }
    return (1/b)*sqrt( aerr*aerr + a*a*berr*berr/(b*b) );
  }

  //return error on a*b
  double errAtimesB(double a, double aerr, double b, double berr) {
    return sqrt( b*b*aerr*aerr + a*a*berr*berr);
  }

  //simple routines for addition in quadrature
  double addInQuad(double a, double b) {return sqrt(a*a + b*b);  }
  double addInQuad(double a, double b, double c) {return sqrt(a*a + b*b + c*c);  }
  double addInQuad(double a, double b, double c, double d) {return sqrt(a*a + b*b + c*c +d*d);  }

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
  
  
  void printHist(const TH1D* h, int binlow=1, int binhigh=-1, bool withErrors=true) {
    if (binhigh == -1) binhigh = h->GetNbinsX();

    std::cout << h->GetName() << ": " ;
    for (int i = binlow; i<= binhigh; i++) {
      if(withErrors) std::cout << format_nevents(h->GetBinContent(i), h->GetBinError(i)) << ", ";
      else std::cout << " " << h->GetBinContent(i);
    }
    std::cout << std::endl;
  }

  
} //end of namespace


