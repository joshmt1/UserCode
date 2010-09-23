#ifndef _PlotUtil_h_
#define _PlotUtil_h_

#include "TLegend.h"

#include "HistHolder.h"

class PlotUtil {
public :

  PlotUtil(HistHolder* hh);
  virtual ~PlotUtil();

  //methods for doing the same thing to a bunch of samples
  void addSample(TString sampleid,TTree* sampletree,TString sampleName="",UInt_t kcolor=1);
  void createHistos(TString varname,int nbins, float hmin, float hmax) ;
  //if drawcommand is empty, then varname is used
  void fillHistos(TString varname, TString cut, TString drawcommand="");
  //    This function takes the histogram called _qcd and adds the _qcdnnn to it
  void addQCD(TString varname);
  //if no varname is specified, the last used varname is used
  void fillLegend(TLegend * leg, TString varname="");
  //will use whatever the current canvas is
  void drawPlots(TString varname="");

  //======methods that calculate things======
  //this is built in to later versions of ROOT, but not 5.22
  Double_t ErrorOnIntegral(const TH1D* h, const Int_t lowbin=1, Int_t highbin=0) ;

  //=====misc settings=====
  void setDebug(bool debug) {debug_=debug;}

 private:
  void setLineColors(TString varname); //called by fillLegend()
  TString getFullName(TString varname,TString sampleid);

  //the TString is the sampleid
  std::map< TString, TTree*> samples_;
  std::map< TString, TString> sampleNames_;
  std::map< TString, UInt_t> sampleColors_;

  HistHolder * hh_;
  bool debug_;
  TString lastVarname_;
  
};
#endif
