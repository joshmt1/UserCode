/*
====== this is the nominal code for drawing RA2b data/MC comparison plots =======
-- drawPlots() -- main plotting routine to draw pretty stacks of MC with data on top
-- drawSimple() -- very simple function to draw exactly one sample and put it in a file
-- drawR() -- draws r(MET). a bit more kludgely, but works.

Various utility functions are at the top (should probably be moved elsewhere for better readability)
Samples to plot are currently hard-coded in loadSamples().

Various routines at the bottom call the above functions. Some of them (e.g. drawSomething() ) I usually
use in a quasi-interactive mode.

-- drawOwen() -- uses drawPlots() and drawSimple() to make a file with the histograms needed for
toys and fits to data.

Details:
This code replaces drawCutflowPlots.C (which had replaced drawBasicPlots.C). 
The input files are the reducedTrees.

Potential improvements:
 -- drawTotalSM_,drawTotalSMSusy_,drawSusyOnly_ don't have a setter function at the moment
 -- renormalizeBins_ doesn't have a setter function at the moment
 -- same with drawMarkers_
 -- the calls to renormBins() currently have a kludge that hard-codes the reference bin to 2
 -- The samples to be plotted are currently defined at compile time.
Could make it so that all samples are always loaded, but there is an
independent list that controls which samples are plotted. This would
allow the user to switch the plotted samples on the fly.
[any mechanism like this would also have to allow the user to control
the order in which the samples are stack. This is currently defined
by the order of the push_backs in loadSamples()]
 -- The cuts are defined by the user settings the selection_ string
directly. This can be error prone. A better interface would allow the user
to set cuts more intutitively and with less change of e.g. accidentally 
forgotten cuts (while still preserving the current flexibility).

 -- someday I should test whether I can get rid of duplicate
functionality for TH1F and TH1D e.g. the case of addOverflowBin()
*/

#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TString.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include "TCut.h"

//not clear yet whether I want to use HistHolder.h or not
#include <fstream>

#include <iostream>
#include <map>

TH1D* hinteractive=0;

//holds a list of the sample names (using same code as file name)
std::vector<TString> samples_;
std::map<TString, TFile*> files_;
std::map<TString, TH1D*> histos_;
std::map<TString, UInt_t> sampleColor_;
std::map<TString, TString> sampleOwenName_;
std::map<TString, TString> sampleLabel_;
std::map<TString, UInt_t> sampleMarkerStyle_;

TString inputPath = "/cu2/joshmt/V00-03-01_2/";
//TString cutdesc = "Baseline0_PF_JERbias_pfMEThigh_PFLep0e0mu_minDP_MuonEcalCleaning";
TString cutdesc = "Baseline0_PF_JERbias_pfMEThigh_PFLepRA20e0mu_minDP_MuonCleaning";
//TString cutdesc = "Baseline0_PF_pfMEThigh_PFLepRA20e0mu_minDP_MuonEcalCleaning";
//TString cutdesc = "Baseline0_PF_pfMEThigh_PFLep0e0mu_minDP_MuonEcalCleaning";
TString selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && cutCleaning==1";

float leg_x1 = 0.696, leg_x2=0.94, leg_y1=0.5, leg_y2=0.92;

bool quiet_=false;
bool doRatio_=false;
bool logy_=false;
bool dostack_=true;
bool doleg_=true;
bool dodata_=true;
bool addOverflow_=true;
//bool doSubtraction_=false;
bool drawQCDErrors_=false;
bool renormalizeBins_=false;

bool normalized_=false;

bool drawTotalSM_=false;
bool drawTotalSMSusy_=false;
bool drawSusyOnly_=false;
bool drawMarkers_=true;

bool doVerticalLine_=false;
double verticalLinePosition_=0;

bool doCustomPlotMax_=false;
double customPlotMax_=0;

bool doCustomPlotMin_=false;
double customPlotMin_=0;

float maxScaleFactor_ = 1.05;

TCanvas* thecanvas=0;
//TCanvas* cratio=0;
TLegend* leg=0;
THStack* thestack=0;
TH1D* totalsm=0;
TH1D* totalsmsusy=0;
TH1D* totalewk=0;
TH1D* totalqcdttbar=0;
TH1D* totalnonttbar=0;
TH1D* totalnonqcd=0;
TH1D* ratio=0; float ratioMin=0; float ratioMax=2;
TGraphErrors* qcderrors=0;
bool loaded_=false; //bookkeeping

void setQuiet(bool q) {
  quiet_ = q;
}

void setQCDErrorMode(bool drawErrors) {
  drawQCDErrors_=drawErrors;
}

void doRatioPlot(bool doIt) {
  doRatio_=doIt;
}

void setPlotMaximum(double max) {
  customPlotMax_=max;
  doCustomPlotMax_=true;
}

void resetPlotMaximum() {
  doCustomPlotMax_=false;
}

void setPlotMinimum(double min) {
  customPlotMin_=min;
  doCustomPlotMin_=true;
}

void resetPlotMinimum() {
  doCustomPlotMin_=false;
}

void enableVerticalLine(double position) {
  doVerticalLine_=true;
  verticalLinePosition_ =position;
}

// void showDataMinusMC(bool dosub) {
//   doSubtraction_=dosub;
// }

void resetVerticalLine() {
  doVerticalLine_=false;
}

void doOverflowAddition(bool doOv) {
  addOverflow_ = doOv;
}

void setLogY(bool dolog) {
  logy_=dolog;
  if (logy_) maxScaleFactor_=3;
  else maxScaleFactor_=1.05;
}

void setStackMode(bool dostack, bool normalized=false) {
  dostack_=dostack;
  normalized_=normalized;
}

void doData(bool dodata) {
  dodata_=dodata;
}

void drawLegend(bool doleg) {
  doleg_=doleg;
}

int mainpadWidth; int mainpadHeight;
int ratiopadHeight = 250;
// TPad* mainPad=0;
// TPad* ratioPad=0;
void renewCanvas(const TString opt="") {
  if (thecanvas!=0) delete thecanvas;

  int canvasWidth = mainpadWidth;
  int canvasHeight = opt.Contains("ratio") ? mainpadHeight+ratiopadHeight : mainpadHeight;

  thecanvas= new TCanvas("thecanvas","the canvas",canvasWidth,canvasHeight);
  thecanvas->cd()->SetRightMargin(0.04);

  if (opt.Contains("ratio")) {
    thecanvas->Divide(1,2);
    const float padding=0.01; const float ydivide=0.2;
    thecanvas->GetPad(1)->SetPad( padding, ydivide + padding, 1-padding, 1-padding);
    thecanvas->GetPad(2)->SetPad( padding, padding, 1-padding, ydivide-padding);
    if (!quiet_)  cout<< thecanvas->GetPad(1)->GetXlowNDC() <<"\t"
		      << thecanvas->GetPad(1)->GetWNDC() <<"\t"
		      << thecanvas->GetPad(1)->GetYlowNDC() <<"\t"
		      << thecanvas->GetPad(1)->GetHNDC() <<endl;
    if (logy_) thecanvas->GetPad(1)->SetLogy();
  }
  else { if (logy_) thecanvas->SetLogy(); }


  int cdarg = opt.Contains("ratio") ? 1 : 0;
  thecanvas->cd(cdarg);

}

void resetPadDimensions() {
  mainpadWidth = 600; 
  mainpadHeight=550;
}

void setPadDimensions(int x, int y) {

  mainpadWidth = x; 
  mainpadHeight= y;
}

void resetHistos() {
  for ( std::map<TString, TH1D*>::iterator i = histos_.begin(); i!=histos_.end(); ++i) {
    if (i->second != 0) {
      delete  i->second;
      i->second= 0;
    }
  }
}

double findOverallMax(const TH1D* hh) {

  double max=-1e9;

  for (int i=1; i<= hh->GetNbinsX(); i++) {
    double val = hh->GetBinContent(i) + hh->GetBinError(i);
    if (val>max) max=val;
  }
  return max;
}

//code largely lifted from Owen
//returned string is the y title of the renormalized histo
TString renormBins( TH1D* hp, int refbin ) {

  if ( hp==0 ) return "PROBLEM";

  double refbinwid = hp->GetBinLowEdge( refbin+1 ) - hp->GetBinLowEdge( refbin ) ;
  if (!quiet_)  printf(" reference bin: [%6.1f,%6.1f], width = %6.3f\n",  hp->GetBinLowEdge( refbin ), hp->GetBinLowEdge( refbin+1 ), refbinwid ) ;
  
  for ( int bi=1; bi<= hp->GetNbinsX(); bi++ ) {
    double binwid = hp->GetBinLowEdge( bi+1 ) - hp->GetBinLowEdge( bi ) ;
    double sf = refbinwid / binwid ;
    if (!quiet_)    printf("  bin %d : width= %6.2f, sf=%7.3f\n", bi, binwid, sf ) ;
    hp->SetBinContent( bi, sf*(hp->GetBinContent( bi )) ) ;
    hp->SetBinError( bi, sf*(hp->GetBinError( bi )) ) ;
  } // bi.

  TString ytitle;
  ytitle.Form("(Events / bin) * (%5.1f / bin width)", refbinwid );

  return ytitle;
}

TString getCutString(TString extraSelection="") {
  TString weightedcut="weight"; 
  if (selection_!="") {
    weightedcut += "*(";
    weightedcut+=selection_;
    if (extraSelection != "") {
      weightedcut += " && ";
      weightedcut +=extraSelection;
    }
    weightedcut+=")";
  }
  else if (extraSelection !="") {
    weightedcut += "*(";
    weightedcut +=extraSelection;
    weightedcut+=")";
  }
  if (!quiet_)  cout<<weightedcut<<endl;
  return weightedcut;
}

void addOverflowBin(TH1D* theHist) {
  //this code was written for when there was a customizable plot range (post-histo creation)
  //it could be made a lot simpler now

  int lastVisibleBin = theHist->GetNbinsX();
  //  cout<<theHist<<"  "<<lastVisibleBin<<"\t";

  //in case there is no custom range, the code should just add the overflow bin to the last bin of the histo
  double lastBinContent = theHist->GetBinContent(lastVisibleBin);
  double lastBinError = pow(theHist->GetBinError(lastVisibleBin),2); //square in prep for addition

  //  cout<<"Overflow addition: "<<lastBinContent<<" +/- "<<sqrt(lastBinError)<<" --> ";

  //now loop over the bins that aren't being shown at the moment (including the overflow bin)
  for (int ibin = lastVisibleBin+1; ibin <= 1 + theHist->GetNbinsX() ; ++ibin) {
    lastBinContent += theHist->GetBinContent( ibin);
    lastBinError += pow(theHist->GetBinError( ibin),2);
  }
  lastBinError = sqrt(lastBinError);

  theHist->SetBinContent(lastVisibleBin,lastBinContent);
  theHist->SetBinError(lastVisibleBin,lastBinError);
  if (!quiet_)  cout<<lastBinContent<<" +/- "<<lastBinError<<endl;
}

void addOverflowBin(TH1F* theHist) {
  //this code was written for when there was a customizable plot range (post-histo creation)
  //it could be made a lot simpler now
  //this one is copied from the function of the same name that takes a TH1D

  int lastVisibleBin = theHist->GetNbinsX();
  //  cout<<theHist<<"  "<<lastVisibleBin<<"\t";

  //in case there is no custom range, the code should just add the overflow bin to the last bin of the histo
  double lastBinContent = theHist->GetBinContent(lastVisibleBin);
  double lastBinError = pow(theHist->GetBinError(lastVisibleBin),2); //square in prep for addition

  //  cout<<"Overflow addition: "<<lastBinContent<<" +/- "<<sqrt(lastBinError)<<" --> ";

  //now loop over the bins that aren't being shown at the moment (including the overflow bin)
  for (int ibin = lastVisibleBin+1; ibin <= 1 + theHist->GetNbinsX() ; ++ibin) {
    lastBinContent += theHist->GetBinContent( ibin);
    lastBinError += pow(theHist->GetBinError( ibin),2);
  }
  lastBinError = sqrt(lastBinError);

  theHist->SetBinContent(lastVisibleBin,lastBinContent);
  theHist->SetBinError(lastVisibleBin,lastBinError);
  if (!quiet_)  cout<<lastBinContent<<" +/- "<<lastBinError<<endl;
}

void drawVerticalLine() {
  if (thecanvas==0) return;

  //this is a fine example of ROOT idiocy
  TVirtualPad* thePad = thecanvas->GetPad(0); //needs fixing for ratio plots
  double xmin,ymin,xmax,ymax;
  thePad->GetRangeAxis(xmin,ymin,xmax,ymax);
  //for academic interest, can get the same numbers using e.g. thePad->GetUymax()
  if (logy_) {
    ymax = pow(10, ymax);
    ymin = pow(10, ymin);
  }
  TLine theLine(verticalLinePosition_,ymin,verticalLinePosition_,ymax);
  theLine.SetLineColor(kBlue);
  theLine.SetLineWidth(3);

  theLine.DrawClone();

}

TFile* fdata=0;
TH1D* hdata=0;
//TH1D* hdataSubtracted=0;
void loadSamples() {
  if (loaded_) return;
  loaded_=true;

  resetPadDimensions();

  //this block controls what samples will enter your plot
  //order of this vector controls order of samples in stack

  //careful -- QCD must have 'QCD' in its name somewhere.
  //samples_.push_back("QCD"); //madgraph
  //samples_.push_back("PythiaQCD");
  samples_.push_back("PythiaPUQCD");
  samples_.push_back("TTbarJets");

  //flip this bool to control whether SingleTop is loaded as one piece or 3
  if (true) samples_.push_back("SingleTop");
  else {
    samples_.push_back("SingleTop-sChannel");
    samples_.push_back("SingleTop-tChannel");
    samples_.push_back("SingleTop-tWChannel");
  }
  samples_.push_back("WJets");
  samples_.push_back("ZJets");
  samples_.push_back("Zinvisible");
  //  samples_.push_back("LM13");

  //these blocks are just a "dictionary"
  //no need to ever comment these out
  if (false) {
    sampleColor_["LM13"] = kRed-9;//kGray; //borrowed from a different sample
    sampleColor_["QCD"] = kYellow;
    sampleColor_["PythiaQCD"] = kYellow;
    sampleColor_["PythiaPUQCD"] = kYellow;
    sampleColor_["PythiaPUQCDFlat"] = kYellow;
    sampleColor_["TTbarJets"]=kRed+1;
    sampleColor_["SingleTop"] = kMagenta;
    sampleColor_["WJets"] = kGreen-3;
    sampleColor_["ZJets"] = kAzure-2;
    sampleColor_["Zinvisible"] = kOrange-3;
    sampleColor_["SingleTop-sChannel"] = kMagenta+1; //for special cases
    sampleColor_["SingleTop-tChannel"] = kMagenta+2; //for special cases
    sampleColor_["SingleTop-tWChannel"] = kMagenta+3; //for special cases
    sampleColor_["TotalSM"] = kBlue+2;
    sampleColor_["Total"] = kGreen+3;
  }
  else { //alternate color scheme requested by Owen
    sampleColor_["LM13"] = kBlue+2;
    sampleColor_["QCD"] = 2;
    sampleColor_["PythiaQCD"] = 2;
    sampleColor_["PythiaPUQCD"] =2;
    sampleColor_["PythiaPUQCDFlat"] =2;
    sampleColor_["TTbarJets"]=4;
    sampleColor_["SingleTop"] = kMagenta;
    sampleColor_["WJets"] = kOrange;
    sampleColor_["ZJets"] = 7;
    sampleColor_["Zinvisible"] = kOrange+7;
    sampleColor_["SingleTop-sChannel"] = kMagenta+1; //for special cases
    sampleColor_["SingleTop-tChannel"] = kMagenta+2; //for special cases
    sampleColor_["SingleTop-tWChannel"] = kMagenta+3; //for special cases
    sampleColor_["TotalSM"] =kGreen+2; //owen requested 3
    sampleColor_["Total"] = 6;
  }

  sampleLabel_["LM13"] = "LM13";
  sampleLabel_["QCD"] = "QCD";
  sampleLabel_["PythiaQCD"] = "QCD (Pythia Z2)";
  sampleLabel_["PythiaPUQCDFlat"] = "QCD (Pileup)"; 
  sampleLabel_["PythiaPUQCD"] = "QCD (Pileup)";
  sampleLabel_["TTbarJets"]="t#bar{t}";
  sampleLabel_["SingleTop"] = "Single-Top";
  sampleLabel_["WJets"] = "W#rightarrowl#nu";
  sampleLabel_["ZJets"] = "Z/#gamma*#rightarrowl^{+}l^{-}";
  sampleLabel_["Zinvisible"] = "Z#rightarrow#nu#nu";
  sampleLabel_["SingleTop-sChannel"] = "Single-Top (s)";
  sampleLabel_["SingleTop-tChannel"] = "Single-Top (t)";
  sampleLabel_["SingleTop-tWChannel"] = "Single-Top (tW)";
  sampleLabel_["TotalSM"] = "SM";
  sampleLabel_["Total"] = "SM + LM13"; //again, this is a hack

  sampleMarkerStyle_["LM13"] = kFullStar;
  sampleMarkerStyle_["QCD"] = kFullCircle;
  sampleMarkerStyle_["PythiaQCD"] = kOpenCircle;
  sampleMarkerStyle_["PythiaPUQCDFlat"] = kOpenCircle;  
  sampleMarkerStyle_["PythiaPUQCD"] = kOpenCircle;
  sampleMarkerStyle_["TTbarJets"]= kFullSquare;
  sampleMarkerStyle_["SingleTop"] = kOpenSquare;
  sampleMarkerStyle_["WJets"] = kMultiply;
  sampleMarkerStyle_["ZJets"] = kFullTriangleUp;
  sampleMarkerStyle_["Zinvisible"] = kFullTriangleDown;
  sampleMarkerStyle_["SingleTop-sChannel"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-tChannel"] = kOpenSquare;
  sampleMarkerStyle_["SingleTop-tWChannel"] = kOpenSquare;
  sampleMarkerStyle_["TotalSM"] = kOpenCross; //FIXME?
  sampleMarkerStyle_["Total"] = kDot; //FIXME?

  sampleOwenName_["LM13"] = "lm13";
  sampleOwenName_["QCD"] = "qcd";
  sampleOwenName_["PythiaQCD"] = "qcd";
  sampleOwenName_["PythiaPUQCDFlat"] = "qcd"; 
  sampleOwenName_["PythiaPUQCD"] = "qcd";
  sampleOwenName_["TTbarJets"]="ttbar";
  sampleOwenName_["SingleTop"] = "singletop";
  sampleOwenName_["WJets"] = "wjets";
  sampleOwenName_["ZJets"] = "zjets";
  sampleOwenName_["Zinvisible"] = "zinvis";
  sampleOwenName_["SingleTop-sChannel"] = "singletops";
  sampleOwenName_["SingleTop-tChannel"] = "singletopt";
  sampleOwenName_["SingleTop-tWChannel"] = "singletoptw";
  sampleOwenName_["TotalSM"] = "totalsm";
  sampleOwenName_["Total"] = "total";  

  for (unsigned int isample=0; isample<samples_.size(); isample++) {
    TString fname="reducedTree.";
    fname+=cutdesc;
    fname+=".";
    fname+=samples_[isample];
    fname+=".root";
    fname.Prepend(inputPath);
    files_[samples_[isample]] = new TFile(fname);
    if (files_[samples_[isample]]->IsZombie() ) cout<<"file error with "<<samples_[isample]<<endl;
    else { if (!quiet_)    cout<<"Added sample: "<<samples_[isample]<<endl;}
  }

  //load data file too
  TString dname="reducedTree.";
  dname+=cutdesc;
  dname+=".data.root";
  dname.Prepend(inputPath);
  if (dname.Contains("JERbias")) dname.ReplaceAll("JERbias_",""); //JERbias not relevant for data
  if ( dodata_) {
    fdata = new TFile(dname);
    if (fdata->IsZombie()) cout<<"Problem with data file! "<<dname<<endl;
  }

}

void renewLegend() {

  if (leg!=0) delete leg;
  leg = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2);
  leg->SetBorderSize(0);
  leg->SetLineStyle(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);

}

//if something is passed to varbins, then low and high will be ignored
float drawSimple(const TString var, const int nbins, const float low, const float high, const TString filename, 
		 const TString histname , const TString samplename, const float* varbins=0) {

  loadSamples();

//I would rather implement this functionality via drawPlots(), but I think it will be simpler
//to just write something simple

//no presentation, just fill the histogram and save
  TTree* tree=0;
  if (samplename=="data") {
    tree = (TTree*) fdata->Get("reducedTree");
  }
  else {
    for (unsigned int isample=0; isample<samples_.size(); isample++) {
      if ( samples_[isample] == samplename) {
	if (!quiet_) cout <<samples_[isample]<<endl;
	tree = (TTree*) files_[samples_[isample]]->Get("reducedTree");
      }
    }
  }
  if (tree==0) {cout<<"Something went wrong finding your sample!"<<endl; return 0;}
  gROOT->cd();
  
  //when owen explicitly uses a histo type, it is a TH1F
  TH1F* hh=0;
  if (varbins==0) {
    hh = new TH1F(histname,histname,nbins,low,high);
  }
  else {
    hh = new TH1F(histname,histname,nbins,varbins);
  }
  hh->Sumw2();
   
  tree->Project(histname,var,getCutString().Data());
  float theIntegral = hh->Integral(0,nbins+1);

  if (addOverflow_)  addOverflowBin( hh ); //manipulates the TH1F

  //at this point i've got a histogram. what more could i want?
  TFile fout(filename,"UPDATE");
  hh->Write();
  fout.Close();
  delete hh; //deleting ROOT objects can be dangerous...but i've tried carefully to avoid a double deletion here by doing gROOT->cd() before the creation of hh
  return theIntegral;
}

float drawSimple(const TString var, const int nbins, const float* varbins, const TString filename, 
		 const TString histname , const TString samplename) {
  return drawSimple(var, nbins, 0, 1, filename, histname, samplename, varbins);
}


void drawPlots(const TString var, const int nbins, const float low, const float high, const TString xtitle, TString ytitle, TString filename="", const float* varbins=0) {
  //  cout<<"[drawPlots] var = "<<var<<endl;

  loadSamples();

  if (filename=="") filename=var;

  //  TH1D* thestackH=0;

  gROOT->SetStyle("CMS");
  //gStyle->SetHatchesLineWidth(1);

  TString canvasOpt = doRatio_ ? "ratio" : "";
  const int mainPadIndex = doRatio_ ? 1 : 0;
  renewCanvas(canvasOpt);

  thecanvas->cd(mainPadIndex);
  renewLegend();

  if (dostack_) {
    if (thestack!= 0 ) delete thestack;
    thestack = new THStack("thestack","--");
    if (doRatio_) {
      if (ratio!=0) delete ratio;
      ratio = (varbins==0) ? new TH1D("ratio","data/(SM MC)",nbins,low,high) : new TH1D("ratio","",nbins,varbins);
      ratio->Sumw2();
    }
  }
  if (totalsm!=0) delete totalsm;
  totalsm = (varbins==0) ? new TH1D("totalsm","",nbins,low,high) : new TH1D("totalsm","",nbins,varbins);
  totalsm->Sumw2();
  if (totalsmsusy!=0) delete totalsmsusy;
  totalsmsusy = (varbins==0) ? new TH1D("totalsmsusy","",nbins,low,high) : new TH1D("totalsmsusy","",nbins,varbins);
  totalsmsusy->Sumw2();
  if (totalewk!=0) delete totalewk;
  totalewk = (varbins==0) ? new TH1D("totalewk","",nbins,low,high) : new TH1D("totalewk","",nbins,varbins);
  totalewk->Sumw2();
  if (totalqcdttbar!=0) delete totalqcdttbar;
  totalqcdttbar = (varbins==0) ? new TH1D("totalqcdttbar","",nbins,low,high) : new TH1D("totalqcdttbar","",nbins,varbins);
  totalqcdttbar->Sumw2();
  if (totalnonttbar!=0) delete totalnonttbar;
  totalnonttbar = (varbins==0) ? new TH1D("totalnonttbar","",nbins,low,high) : new TH1D("totalnonttbar","",nbins,varbins);
  totalnonttbar->Sumw2();
  if (totalnonqcd!=0) delete totalnonqcd;
  totalnonqcd = (varbins==0) ? new TH1D("totalnonqcd","",nbins,low,high) : new TH1D("totalnonqcd","",nbins,varbins);
  totalnonqcd->Sumw2();

  totalsm->SetMarkerColor(sampleColor_["TotalSM"]);
  totalsm->SetLineColor(sampleColor_["TotalSM"]);
  totalsm->SetLineWidth(2);
  totalsm->SetMarkerStyle(sampleMarkerStyle_["TotalSM"]);
  if (!drawMarkers_)  totalsm->SetMarkerSize(0); //no marker for this one

  totalsmsusy->SetMarkerColor(sampleColor_["Total"]);
  totalsmsusy->SetLineColor(sampleColor_["Total"]);
  totalsmsusy->SetLineWidth(2);
  totalsmsusy->SetMarkerStyle(sampleMarkerStyle_["Total"]);
  if (!drawMarkers_)  totalsmsusy->SetMarkerSize(0); //no marker for this one

  //here is the part that is really different from the previous implementation
  //need to make new histograms
  resetHistos(); //delete existing histograms
  TString opt="hist e";
  double histMax=-1e9;
  for (unsigned int isample=0; isample<samples_.size(); isample++) {
    if (!quiet_)   cout <<samples_[isample]<<endl;

    gROOT->cd();
    //should each histo have a different name? maybe
    TString hname = var; hname += "_"; hname += samples_[isample];
    histos_[samples_[isample]] = (varbins==0) ? new TH1D(hname,"",nbins,low,high) : new TH1D(hname,"",nbins,varbins);
    histos_[samples_[isample]]->Sumw2();

    //qcd reweighting not implemented yet

    TTree* tree = (TTree*) files_[samples_[isample]]->Get("reducedTree");
    gROOT->cd();
    tree->Project(hname,var,getCutString().Data());
    //now the histo is filled
    
    if (renormalizeBins_) ytitle=renormBins(histos_[samples_[isample]],2 ); //manipulates the TH1D //FIXME hard-coded "2"
    if (addOverflow_)  addOverflowBin( histos_[samples_[isample]] ); //manipulates the TH1D
    histos_[samples_[isample]]->SetXTitle(xtitle);
    histos_[samples_[isample]]->SetYTitle(ytitle);

    //if we're going to draw QCD errors, create a TGraphErrors from the QCD histogram
    if (drawQCDErrors_ && samples_[isample].Contains("QCD")) {
      if (qcderrors!=0) delete qcderrors;
      qcderrors = new TGraphErrors(histos_[samples_[isample]]);
      qcderrors->SetFillStyle(3353);
      qcderrors->SetFillColor(1);
    }
    if (!samples_[isample].Contains("LM")) {
      totalsm->Add(histos_[samples_[isample]]);
      if (!quiet_)    cout << "totalsm: " << samples_[isample] << endl;
    }
    if (!samples_[isample].Contains("LM") && !samples_[isample].Contains("QCD") && !samples_[isample].Contains("TTbar")) {
      totalewk->Add(histos_[samples_[isample]]);
      if (!quiet_) cout << "totalewk: " << samples_[isample] << endl;
    }
    if (samples_[isample].Contains("QCD") || samples_[isample].Contains("TTbar")){
      totalqcdttbar->Add(histos_[samples_[isample]]);
      if (!quiet_) cout << "totalqcdttbar: " << samples_[isample] << endl;
    }
    if (!samples_[isample].Contains("TTbar") && !samples_[isample].Contains("LM")){
      totalnonttbar->Add(histos_[samples_[isample]]);
      if (!quiet_) cout << "totalnonttbar: " << samples_[isample] << endl;
    }
    if (!samples_[isample].Contains("QCD") && !samples_[isample].Contains("LM")){
       totalnonqcd->Add(histos_[samples_[isample]]);
      if (!quiet_) cout << "totalnonqcd: " << samples_[isample] << endl;
    }
    totalsmsusy->Add(histos_[samples_[isample]]); //add everything!

    //now just do a bunch of histogram formatting
    if (!dostack_) {
      //set line color instead of fill color for this type of plot
      histos_[samples_[isample]]->SetLineColor(sampleColor_[samples_[isample]]);
      histos_[samples_[isample]]->SetMarkerStyle(sampleMarkerStyle_[samples_[isample]]);
      histos_[samples_[isample]]->SetMarkerColor(sampleColor_[samples_[isample]]);
      if (!drawMarkers_) histos_[samples_[isample]]->SetMarkerSize(0);

      //ad hoc additions
      histos_[samples_[isample]]->SetLineWidth(2);
    }
    else {
      histos_[samples_[isample]]->SetFillColor(sampleColor_[samples_[isample]]);
      histos_[samples_[isample]]->SetMarkerSize(0);
    }

    if (dostack_) { //add histo to stack
      leg->AddEntry(histos_[samples_[isample]], sampleLabel_[samples_[isample]]);
      thestack->Add(histos_[samples_[isample]] );
    }
    else { //draw non-stacked histo
      //normalize
      if ( normalized_ && histos_[samples_[isample]]->Integral() >0)  histos_[samples_[isample]]->Scale( 1.0 / histos_[samples_[isample]]->Integral() );
      if (!drawSusyOnly_ || samples_[isample].Contains("LM")) { //drawSusyOnly_ means don't draw SM
	//set max
	if ( findOverallMax( histos_[samples_[isample]]) > histMax) histMax = findOverallMax(histos_[samples_[isample]]);
	
	leg->AddEntry(histos_[samples_[isample]], sampleLabel_[samples_[isample]]);

	histos_[samples_[isample]]->Draw(opt);
	if (!opt.Contains("same")) opt+=" same";
      }
    }
  } //loop over samples and fill histograms

  if (drawTotalSM_) leg->AddEntry(totalsm, sampleLabel_["TotalSM"]);
  if (drawTotalSMSusy_) leg->AddEntry(totalsmsusy, sampleLabel_["Total"]);

  if (!dostack_) {
    //this is all a re-implemenataion of stuff done is HistHolder. Oh well.

    if (drawTotalSM_) histMax = totalsm->GetMaximum();
    for (unsigned int isample=0; isample<samples_.size(); isample++) {
      double pmx= doCustomPlotMax_ ? customPlotMax_ : histMax*maxScaleFactor_;
      histos_[samples_[isample]]->SetMaximum(pmx);
      if (doCustomPlotMin_) histos_[samples_[isample]]->SetMinimum(customPlotMin_);
    }

    if (drawTotalSM_) { 
      totalsm->Draw(opt); 
      if (doCustomPlotMax_) totalsm->SetMaximum(customPlotMax_);
    }
    if (drawTotalSMSusy_) {
      totalsmsusy->Draw(opt); 
      if (doCustomPlotMax_) totalsmsusy->SetMaximum(customPlotMax_);
    }
  }
  else {
    thestack->Draw("hist");
    thestack->GetHistogram()->GetXaxis()->SetTitle(xtitle);
    thestack->GetHistogram()->GetYaxis()->SetTitle(ytitle);

    if (doVerticalLine_) drawVerticalLine(); //i want to draw the data last

    if (drawQCDErrors_) qcderrors->Draw("2 same");

    if (doCustomPlotMax_) thestack->SetMaximum(customPlotMax_);
    if (doCustomPlotMin_) thestack->SetMinimum(customPlotMin_);
  } //if doStack_

  if (dodata_) {
    gROOT->cd();
    if (!quiet_)     cout<<"Drawing data!"<<endl;
    if (hdata != 0) delete hdata;
    TString hname = var; hname += "_"; hname += "data";
    hdata = (varbins==0) ? new TH1D(hname,"",nbins,low,high) : new TH1D(hname,"",nbins,varbins);
    hdata->Sumw2();
    TTree* dtree = (TTree*) fdata->Get("reducedTree");
    gROOT->cd();
    dtree->Project(hname,var,selection_.Data());
    //now the histo is filled
    
    hdata->UseCurrentStyle(); //maybe not needed anymore
    hdata->SetMarkerColor(kBlack);
    hdata->SetLineWidth(2);
    hdata->SetMarkerStyle(kFullCircle);
    hdata->SetMarkerSize(1);
    if (renormalizeBins_) renormBins(hdata,2 ); //manipulates the histogram //FIXME hard-coded "2"
    if (addOverflow_)     addOverflowBin(hdata); // manipulates the histogram!
    leg->AddEntry(hdata,"Data");

    hdata->Draw("SAME");
    if (!doCustomPlotMax_) {
      double mymax = dostack_ ? thestack->GetMaximum() : findOverallMax(totalsm); //these probably return near-identical values, in fact
      if (findOverallMax(hdata) > mymax) {
	if (dostack_) thestack->SetMaximum( maxScaleFactor_*findOverallMax(hdata));
	else { //i don't like repeating this loop constantly; at a minimum it should be abstracted
	  for (unsigned int isample=0; isample<samples_.size(); isample++)  histos_[samples_[isample]]->SetMaximum(maxScaleFactor_*findOverallMax(hdata));
	}
      }
    }

    if (!quiet_ && !renormalizeBins_) {
      cout<<"Integral of data, EW, total SM: "<<hdata->Integral()<<" ; "<<totalewk->Integral()<<" ; "<<totalsm->Integral()<<endl;
      cout<<"Chi^2 Test results: "<<hdata->Chi2Test(totalsm,"UW P")<<endl;
      cout<<"KS Test results: "<<hdata->KolmogorovTest(totalsm,"N")<<endl;;
    }
    if (doRatio_) {
      thecanvas->cd(2);
      ratio->Divide(hdata,totalsm);
      ratio->SetMinimum(ratioMin);
      ratio->SetMaximum(ratioMax);
      ratio->GetYaxis()->SetNdivisions(200 + int(ratioMax-ratioMin)+1);    //set ticks ; to be seen if this really works
      ratio->GetYaxis()->SetLabelSize(0.2); //make y label bigger
      ratio->Draw();
      thecanvas->GetPad(2)->SetTopMargin(0.1);
    }
  }

  thecanvas->cd(mainPadIndex);
  if (doleg_)  leg->Draw();

  //  if (doSubtraction_) savename+="-MCSub";
  TString savename = filename;
  if (logy_) savename += "-logY";
  //  savename += scaleAppendToFilename;

  if (!dostack_ && !normalized_)      savename += "-drawPlain";
  else if (!dostack_ && normalized_)  savename += "-drawNorm";
  else savename += "-drawStack";

  //amazingly, \includegraphics cannot handle an extra dot in the filename. so avoid it.
  thecanvas->SaveAs(savename+".eps"); //for me
  //  thecanvas->Print(savename+".C");    //for formal purposes
  thecanvas->SaveAs(savename+".pdf"); //for pdftex
  thecanvas->SaveAs(savename+".png"); //for twiki

}

void drawPlots(const TString var, const int nbins, const float* varbins, const TString xtitle, const TString ytitle, TString filename="") {
  //provide a more natural modification to the argument list....
  drawPlots( var, nbins, 0, 1, xtitle,ytitle, filename,  varbins);
}


//could add xtitle and ytitle
void drawR(const TString vary, const float cutVal, const int nbins, const float low, const float high, const TString& savename) {
  const TString ytitle="N pass / N fail";

  const TString var = "MET"; //hardcoded for now
  TString cstring1 = vary, cstring2=vary;
  cstring1 += " >= ";
  cstring2 += " < ";
  cstring1 += cutVal;
  cstring2 += cutVal;

  loadSamples();

  gROOT->SetStyle("CMS");

  renewCanvas("ratio");

  //in first incarnation, make a separate r(MET) plot for each sample in the list

//   TH1D* qcdPass = new TH1D("qcdPass","",nbins,low,high);
//   TH1D* qcdFail = new TH1D("qcdFail","",nbins,low,high);
//   TH1D* qcdRatio = new TH1D("qcdRatio","",nbins,low,high);

//   qcdPass->Sumw2();
//   qcdFail->Sumw2();
//   qcdRatio->Sumw2();

  resetHistos(); //delete existing histograms

  renewLegend();


  // === begin correlation hack ====
  gROOT->cd();
  TH2D totalsm2d_50("totalsm2d_50","",50,50,100,50,0,TMath::Pi());
  TH2D totalsm2d_SB("totalsm2d_SB","",50,100,150,50,0,TMath::Pi());
  totalsm2d_50.Sumw2();
  totalsm2d_SB.Sumw2();
  TH2D data2d_50("data2d_50","",50,50,100,50,0,TMath::Pi());
  TH2D data2d_SB("data2d_SB","",50,100,150,50,0,TMath::Pi());
  data2d_50.Sumw2();
  data2d_SB.Sumw2();
  // === end correlation hack ===

  TH1D  totalsm_pass("totalsm_pass","",nbins,low,high);
  TH1D  totalsm_fail("totalsm_fail","",nbins,low,high);
  totalsm_pass.Sumw2(); 
  totalsm_fail.Sumw2(); 
  if (totalsm!=0) delete totalsm;
  totalsm =  new TH1D("totalsm","",nbins,low,high);
  totalsm->Sumw2();

  totalsm->SetMarkerColor(sampleColor_["TotalSM"]);
  totalsm->SetLineColor(sampleColor_["TotalSM"]);
  totalsm->SetLineWidth(2);
  totalsm->SetMarkerStyle(0);
  totalsm->SetYTitle(ytitle);

  TString drawopt="hist e";
  float max=-1e9; TString firsthist="";
  for (unsigned int isample=0; isample<samples_.size(); isample++) {

    if (!quiet_) cout <<samples_[isample]<<endl;
    TTree* tree = (TTree*) files_[samples_[isample]]->Get("reducedTree");

    gROOT->cd();

    //need Pass, Fail, and Ratio for each sample
    TString hnameP = var; hnameP += "_"; hnameP += samples_[isample];
    hnameP += "_Pass";
    histos_[hnameP] = new TH1D(hnameP,"",nbins,low,high);
    histos_[hnameP]->Sumw2();

    TString hnameF = var; hnameF += "_"; hnameF += samples_[isample];
    hnameF += "_Fail";
    histos_[hnameF] = new TH1D(hnameF,"",nbins,low,high);
    histos_[hnameF]->Sumw2();

    TString hnameR = var; hnameR += "_"; hnameR += samples_[isample];
    hnameR += "_Ratio";
    histos_[hnameR] = new TH1D(hnameR,"",nbins,low,high);
    histos_[hnameR]->Sumw2();

    //Fill histos
    tree->Project(hnameP,var,getCutString(cstring1).Data());
    tree->Project(hnameF,var,getCutString(cstring2).Data());

    if (addOverflow_)  addOverflowBin( histos_[hnameP] );
    if (addOverflow_)  addOverflowBin( histos_[hnameF] );

    //compute ratio
    histos_[hnameR]->Divide(histos_[hnameP], histos_[hnameF]);

    if (!samples_[isample].Contains("LM")) {
      totalsm_pass.Add(histos_[hnameP]);
      totalsm_fail.Add(histos_[hnameF]);

      TH2D this2d_SB("this2d_SB","",50,100,150,50,0,TMath::Pi());
      tree->Project("this2d_SB","minDeltaPhi:MET",getCutString().Data());
      TH2D this2d_50("this2d_50","",50,50,100,50,0,TMath::Pi());
      tree->Project("this2d_50","minDeltaPhi:MET",getCutString().Data());

      totalsm2d_SB.Add(&this2d_SB);
      totalsm2d_50.Add(&this2d_50);
    }

    //   cout<<"content of bin 2: "<<histos_[hnameP]->GetBinContent(2)<<" / "<< histos_[hnameF]->GetBinContent(2)<<" = "<<histos_[hnameR]->GetBinContent(2)<<endl;

    //now format the histograms
    if (!quiet_) cout<<"setting color to: "<<sampleColor_[samples_[isample]]<<endl;
    histos_[hnameR]->SetLineColor(sampleColor_[samples_[isample]]);
    histos_[hnameR]->SetMarkerStyle(sampleMarkerStyle_[samples_[isample]]);
    histos_[hnameR]->SetMarkerColor(sampleColor_[samples_[isample]]);
    histos_[hnameR]->SetYTitle(ytitle);

    //ad hoc additions
    histos_[hnameR]->SetLineWidth(2);

    //draw
    thecanvas->cd(1);
    if (hnameR.Contains("QCD")) { //HACK draw only qcd
      histos_[hnameR]->Draw(drawopt);
      if (!drawopt.Contains("same")) drawopt+=" same";
      
      if (firsthist="") firsthist = hnameR;
      if (histos_[hnameR]->GetMaximum() > max) max = histos_[hnameR]->GetMaximum();
      leg->AddEntry(histos_[hnameR], sampleLabel_[samples_[isample]]);
    }

  }

  histos_[firsthist]->SetMaximum( max*maxScaleFactor_);
  hinteractive =  histos_[firsthist];

  totalsm->Divide(&totalsm_pass,&totalsm_fail);
  if (drawTotalSM_) {
    totalsm->Draw("hist e same");
    //    leg->Clear();
    leg->AddEntry(totalsm,sampleLabel_["TotalSM"]);
  }

  if (dodata_) {
    gROOT->cd();
    if (!quiet_)   cout<<"Drawing data!"<<endl;
    if (hdata != 0) delete hdata;

    TString hname = var; hname += "_"; hname += "data";
    hdata = new TH1D(hname,"",nbins,low,high);
    hdata->Sumw2();

    TString hnameP = var; hnameP += "_"; hnameP += "dataPass";
    histos_[hnameP] = new TH1D(hnameP,"",nbins,low,high);
    histos_[hnameP]->Sumw2();

    TString hnameF = var; hnameF += "_"; hnameF += "dataFail";
    histos_[hnameF] = new TH1D(hnameF,"",nbins,low,high);
    histos_[hnameF]->Sumw2();

    TTree* dtree = (TTree*) fdata->Get("reducedTree");
    gROOT->cd();
    dtree->Project(hnameP,var,getCutString(cstring1).Data());
    dtree->Project(hnameF,var,getCutString(cstring2).Data());
    if (addOverflow_)  addOverflowBin( histos_[hnameP] );
    if (addOverflow_)  addOverflowBin( histos_[hnameF] );
    //compute ratio
    hdata->Divide(histos_[hnameP], histos_[hnameF]);

    dtree->Project("data2d_SB","minDeltaPhi:MET",getCutString().Data());
    dtree->Project("data2d_50","minDeltaPhi:MET",getCutString().Data());

    //    hdata->UseCurrentStyle(); //maybe not needed anymore
    hdata->SetMarkerColor(kBlack);
    hdata->SetLineWidth(2);
    hdata->SetMarkerStyle(kFullCircle);
    hdata->SetMarkerSize(1);

    thecanvas->cd(1);
    hdata->Draw("SAME");
    leg->AddEntry(hdata,"Data");

    if (hdata->GetMaximum() > max)  {
      histos_[firsthist]->SetMaximum( maxScaleFactor_*hdata->GetMaximum());
      totalsm->SetMaximum(maxScaleFactor_*hdata->GetMaximum());
    }
    else if (doCustomPlotMax_) {
      histos_[firsthist]->SetMaximum( customPlotMax_);
      totalsm->SetMaximum(customPlotMax_);
    }
    if (doCustomPlotMin_) {
      histos_[firsthist]->SetMinimum( customPlotMin_);
      totalsm->SetMinimum(customPlotMin_);
    }

    //    cratio->cd();
    thecanvas->cd(2);
    if (ratio!=0) delete ratio;
    ratio = new TH1D("ratio","data/(SM MC)",nbins,low,high);
    ratio->Sumw2();
    ratio->Divide(hdata,totalsm); 
    ratio->SetMinimum(ratioMin);
    ratio->SetMaximum(ratioMax);
    ratio->Draw();
    cout<<"KS Test results (shape only): "<<hdata->KolmogorovTest(totalsm)<<endl;;
  }
  thecanvas->cd(1);
  leg->Draw();

  thecanvas->SaveAs("mindpPassOverFail-"+savename+".eps");
  thecanvas->SaveAs("mindpPassOverFail-"+savename+".pdf");
  thecanvas->SaveAs("mindpPassOverFail-"+savename+".png");

//   TCanvas* c2d=new TCanvas("c2d","2d",800,800);
//   c2d->Divide(2,2);
//   c2d->cd(1);
//   totalsm2d_50.DrawCopy("colz");
//   c2d->cd(2);
//   totalsm2d_SB.DrawCopy("colz");
//   c2d->cd(3);
//   data2d_50.DrawCopy("colz");
//   c2d->cd(4);
//   data2d_SB.DrawCopy("colz");
  cout<<"Total SM MC correlation [50<MET<100]  = "<<totalsm2d_50.GetCorrelationFactor()<<endl;
  cout<<"Total SM MC correlation [100<MET<150] = "<<totalsm2d_SB.GetCorrelationFactor()<<endl;
  cout<<"Data correlation [50<MET<100]         = "<<data2d_50.GetCorrelationFactor()<<endl;
  cout<<"Data correlation [100<MET<150]        = "<<data2d_SB.GetCorrelationFactor()<<endl;

}

void drawrplots() {
  /*
.L drawReducedTrees.C++
  */

  setStackMode(false);
  doData(true);

  //no met, no mindeltaphi
  setPlotMinimum(1e-2); //setPlotMaximum(1);

 drawTotalSM_=true;
 leg_y1=0.6;
 setLogY(true);

 selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";
 drawR("minDeltaPhi",0.3,20,0,200,"ge1b");

 selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";
 drawR("minDeltaPhi",0.3,20,0,200,"eq1b");

 selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";
 drawR("minDeltaPhi",0.3,20,0,200,"ge2b");

 //coarse
 selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";
 drawR("minDeltaPhi",0.3,4,0,200,"ge1b-4bins");

 selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";
 drawR("minDeltaPhi",0.3,4,0,200,"eq1b-4bins");

 selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";
 drawR("minDeltaPhi",0.3,4,0,200,"ge2b-4bins");
//  const double epsilon=-1;
//  double cbias= hinteractive->GetBinContent(4)/hinteractive->GetBinContent(3);
//  double cbias_syst = (1-cbias)/2;
//  double cbias_up = cbias+cbias_syst;
//  double cbias_down = cbias-cbias_syst;

//  double r_up = hinteractive->GetBinContent(3) * cbias_up;
//  double r_down = hinteractive->GetBinContent(3) * cbias_down;
//  TLine lup(hinteractive->GetBinCenter(4)+epsilon, r_down, hinteractive->GetBinCenter(4)+epsilon, r_up);
//  lup.SetLineColor(kMagenta);
//  lup.SetLineWidth(2);
//  lup.Draw("SAME");

//not so useful
 selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000 && bestTopMass>300";
 drawR("minDeltaPhi",0.3,4,0,200,"ge2b-m3j300");

 //more useful? maybe
 selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000 && bestTopMass>450";
 drawR("minDeltaPhi",0.3,4,0,200,"ge1b-m3j450");

 selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000 && bestTopMass>450 && MET>=100&&MET<150";
 drawR("minDeltaPhi",0.3,5,100,150,"ge1b-SBonly-m3j450");


 //njets == 3
// selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && cutCleaning==1 && njets==3";
// setLogY(true);
//  drawR("minDeltaPhi",0.3,50,0,250);

//  //4 jets
// selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && cutCleaning==1 && njets==4";
// setLogY(true);
//  drawR("minDeltaPhi",0.3,50,0,250);

//  //5 jets
// selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && cutCleaning==1 && njets==5";
// setLogY(true);
//  drawR("minDeltaPhi",0.3,50,0,250);

}

void drawForPreAp_smonly() {
  /*
comment out LM13

.L drawReducedTrees.C++
  */
  setStackMode(false,false);
  doData(true);
  doRatioPlot(false);
  //  setPadDimensions(700,500);
  setLogY(true);

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  drawTotalSM_=true;
  drawMarkers_=false;

  var="MET"; xtitle="E_{T}^{miss} [GeV]";

  //log scale, no stack, full MET range
  nbins=35; low=0; high=350;
  setPlotMinimum(0.1);
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_ge1b");
  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_eq1b");
  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_ge2b");

  //now linear scale with custom y range
  resetPlotMinimum();
  setLogY(false);

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  setPlotMaximum(35);
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_ge1b");
  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  setPlotMaximum(25);
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_eq1b");
  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  setPlotMaximum(11);
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_ge2b");

  //now fail minDeltaPhi for high-ish MET
  resetPlotMinimum();
  setLogY(false);
  nbins=35; low=0; high=350;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==0 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>0";
  setPlotMaximum(125);
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_IDP_ge1b");
  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==0 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>0";
  setPlotMaximum(105);
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_IDP_eq1b");
  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==0 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>0";
  setPlotMaximum(23);
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_IDP_ge2b");

  setLogY(true);
  setPlotMinimum(0.1);
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==0 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>0";
  setPlotMaximum(150);
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_IDP_ge1b");
  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==0 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>0";
  setPlotMaximum(130);
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_IDP_eq1b");
  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==0 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>0";
  setPlotMaximum(50);
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_IDP_ge2b");

}

void drawForPreAp_susy() {

  /*
uncomment LM13

.L drawReducedTrees.C++
  */
  setStackMode(false,false);
  doData(true);
  doRatioPlot(false);
  //  setPadDimensions(700,500);
  setLogY(true);

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  drawTotalSM_=true;
  drawTotalSMSusy_=true;
  drawSusyOnly_=true;
  drawMarkers_=false;

  nbins=60; low=0; high=600;
  setPlotMinimum(0.1);
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_withSusy_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_withSusy_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_withSusy_ge2b");

  setLogY(false);
  resetPlotMinimum();

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  setPlotMaximum(11);
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_withSusy_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  setPlotMaximum(7);
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_withSusy_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  setPlotMaximum(5);
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_withSusy_ge2b");



}

void drawForANandPAS () {
  /*
.L drawReducedTrees.C++
  */
  setStackMode(true);
  doData(true);
  doRatioPlot(false);
  resetPlotMinimum();
  setPadDimensions(700,500);
  setLogY(false);

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(false);

  var="bestTopMass"; xtitle="best 3-jet mass (GeV)";
  // === plots for AN and PAS (aspect ratio adjusted)
  const int nvarbins=5;
  const float varbins[]={0.,160.,180.,260.,400.,800.};
  //SB region
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>=100 && MET<150";
  renormalizeBins_=false;
  setPlotMaximum(25);
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "mcdata_bestM3j_met_100_150_ge1btag_vb");
  renormalizeBins_=true;
  resetPlotMaximum();
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "mcdata_bestM3j_met_100_150_ge1btag_vbrn");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>=100 && MET<150";
  renormalizeBins_=false;
  setPlotMaximum(20);
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "mcdata_bestM3j_met_100_150_eq1btag_vb");
  renormalizeBins_=true;
  resetPlotMaximum();
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "mcdata_bestM3j_met_100_150_eq1btag_vbrn");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>=100 && MET<150";
  renormalizeBins_=false;
  //  setPlotMaximum(8);
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "mcdata_bestM3j_met_100_150_ge2btag_vb");
  renormalizeBins_=true;
  resetPlotMaximum();
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "mcdata_bestM3j_met_100_150_ge2btag_vbrn");

  //-------------- signal region ---------------
  nbins=10; low=0; high=800;

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && cutMET==1";
  renormalizeBins_=false;
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "mcdata_bestM3j_met_50_inf_ge1btag_vb");
  drawPlots(var,nbins,low,high,xtitle,"Events", "mcdata_bestM3j_met_50_inf_ge1btag");
  renormalizeBins_=true;
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "mcdata_bestM3j_met_50_inf_ge1btag_vbrn");


  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && cutMET==1";
  renormalizeBins_=false;
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "mcdata_bestM3j_met_50_inf_eq1btag_vb");
  drawPlots(var,nbins,low,high,xtitle,"Events", "mcdata_bestM3j_met_50_inf_eq1btag");
  renormalizeBins_=true;
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "mcdata_bestM3j_met_50_inf_eq1btag_vbrn");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && cutMET==1";
  renormalizeBins_=false;
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "mcdata_bestM3j_met_50_inf_ge2btag_vb");
  drawPlots(var,nbins,low,high,xtitle,"Events", "mcdata_bestM3j_met_50_inf_ge2btag");
  renormalizeBins_=true;
  drawPlots(var,nvarbins,varbins,xtitle,"Events", "mcdata_bestM3j_met_50_inf_ge2btag_vbrn");

  // -- SB region, compare shapes --
  setStackMode(false);
  renormalizeBins_=false;
  nbins=10; low=0; high=800;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>=100 && MET<150";
  drawPlots(var,nbins,low,high,xtitle,"Arbitrary Units", "m3j_shapesInMC_SB_ge1btag");
  drawPlots(var,nvarbins,varbins,xtitle,"Arbitrary Units", "m3j_shapesInMC_SB_ge1btag_vb");


  resetPadDimensions();

}


void drawSomething() {
  /* for reference, here is all cuts (>=0 b)
selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && cutCleaning==1";

//here is all cuts but with ECAL cleaning removed
selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";

.L drawReducedTrees.C++

  */

  //here we imitate exactly the plots we made in plotWithScale() of the old code

  setStackMode(true);
  doData(true);
  doRatioPlot(true);

  int nbins;
  float low,high;
  TString var,xtitle;

  //  setQCDErrorMode(true);
  // PROBLEMS -- CMS style ruins the hash marks. central value on the TGraphErrors is wrong!
  //might be easiest to accumulate a total MC uncertainty and plot that!

  //the strength and weakness of this tree-based method is that I need to apply the cuts now!
  //most straightforward way is to manually set the cut string before each plot


  // ==== MET plots ====
  setLogY(true);   setPlotMinimum(1e-1);
  nbins=25; low= 0; high=260;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  ratioMin = 0; ratioMax = 2;

  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET");
  //  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_JER0"); //in case we're doing the non-JERbias version

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_ge2b");

  setLogY(true);   setPlotMinimum(1e-1);
  nbins=50; low= 0; high=500;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  ratioMin = 0.5; ratioMax = 1.5;
  selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_METw");

  // == MET plot in limited range (SB and up) ==
  setLogY(false);   setPlotMinimum(0);
  doRatioPlot(false);
  nbins=10; low= 100; high=200;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_SB_ge1b");

  // == MET plot in fail minDeltaPhi
  nbins=10; low= 100; high=200;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==0 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_SBfail_ge1b");

  doRatioPlot(true);
  // === MET plot in ILV/T2 region ===
  setLogY(false);   setPlotMinimum(0);
  nbins=30; low= 0; high=300;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  ratioMin = 0; ratioMax = 2;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1) || (nElectrons==1 && nMuons==0)) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_ILV_ge1b");

  nbins=20; low= 0; high=300;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && (nElectrons==1 && nMuons==0) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_IEV_ge1b");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && (nElectrons==0 && nMuons==1) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_IMV_ge1b");

  // == replacements for Don's ILV plots (but without fancy ttbar categories)
  //these are in the MET>150 region
  doRatioPlot(false);
  setLogY(false);   setPlotMinimum(0);
  nbins=6; low=1; high=7;
  var="njets"; xtitle="Jet multiplicity";
  //ratioMin = 0; ratioMax = 2;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && (nElectrons==0 && nMuons==1) && cutDeltaPhi==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "njet_imv_1btag");
  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && (nElectrons==0 && nMuons==1) && cutDeltaPhi==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "njet_imv_eq1btag");
  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && (nElectrons==0 && nMuons==1) && cutDeltaPhi==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "njet_imv_2btag");

  doRatioPlot(true);
  setLogY(false);   setPlotMinimum(0);
  nbins=10; low= 0; high=200;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  ratioMin = 0; ratioMax = 2;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1  && (nElectrons==0 && nMuons==1) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_imv_1btag");
  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1  && (nElectrons==0 && nMuons==1) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_imv_eq1btag");
  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1  && (nElectrons==0 && nMuons==1) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_imv_2btag");

  doRatioPlot(false);
  setLogY(false);   setPlotMinimum(0);
  nbins=6; low=1; high=7;
  var="njets"; xtitle="Jet multiplicity";
  //ratioMin = 0; ratioMax = 2;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && (nElectrons==1 && nMuons==0) && cutDeltaPhi==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "njet_iev_1btag");
  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && (nElectrons==1 && nMuons==0) && cutDeltaPhi==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "njet_iev_eq1btag");
  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && (nElectrons==1 && nMuons==0) && cutDeltaPhi==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "njet_iev_2btag");

  doRatioPlot(true);
  setLogY(false);   setPlotMinimum(0);
  nbins=10; low= 0; high=200;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  ratioMin = 0; ratioMax = 2;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1  && (nElectrons==1 && nMuons==0) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_iev_1btag");
  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1  && (nElectrons==1 && nMuons==0) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_iev_eq1btag");
  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1  && (nElectrons==1 && nMuons==0) && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "met_iev_2btag");


  // ==== MHT plot === (no MET cut)
  setLogY(true);  setPlotMinimum(1e-1);
  nbins=50; low= 0; high=500;
  ratioMin = 0.5; ratioMax = 1.5;
  var="MHT"; xtitle="H_{T}^{miss} [GeV]";

  selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MHT");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MHT_ge1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MHT_ge2b");

  //to compare to RA2
  setLogY(true);  setPlotMinimum(1);
  nbins=16; low= 0; high=160;
  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MHTveryloose");


  // ==== minDeltaPhi plots ====
  setLogY(false);
  resetPlotMinimum();
  nbins=10; low=0; high = TMath::Pi() + 0.001;
  var="minDeltaPhi"; xtitle="min(#Delta#phi[ jets 1..3, E_{T}^{miss} ] )";

  selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HminDeltaPhiMETj");
  //drawPlots(var,nbins,low,high,xtitle,"Events", "HminDeltaPhiMETj_JER0");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HminDeltaPhiMETj_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HminDeltaPhiMETj_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HminDeltaPhiMETj_ge2b");

  // ==== HT plots ==== //in this case this is not an N-1 plot, because we've already cut on HT
  ratioMin = 0; ratioMax = 3;
  resetPlotMinimum();
  setLogY(false);
  nbins=10; low=300; high = 1000;
  var="HT"; xtitle="H_{T} (GeV)";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_HT_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_HT_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_HT_ge2b");

  // ==== n jets ====
  nbins=8; low=0; high = 8;
  resetPlotMinimum();
  setLogY(false);
  var="njets"; xtitle="Jet multiplicity";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "Hnjets_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "Hnjets_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 &&  MET>=100 && MET<150&& cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "Hnjets_ge2b");

  // ==== n b jets ====
  ratioMin = 0; ratioMax = 2;
  nbins=4; low=0; high = 4;
  resetPlotMinimum();
  setLogY(false);
  var="nbjets"; xtitle="Number of b-tagged jets";
  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && cut3Jets==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "Hnbjets");

  // ==== n Electrons ====
  nbins=3; low=0; high= 3;
  resetPlotMinimum();
  setLogY(false);
  var="nElectrons"; xtitle="Number of electrons";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HnElectrons_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HnElectrons_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HnElectrons_ge2b");

  // ==== n Muons ====
  nbins=3; low=0; high= 3;
  resetPlotMinimum();
  setLogY(false);
  var="nMuons"; xtitle="Number of muons";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HnMuons_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HnMuons_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HnMuons_ge2b");

  // ==== lead jet pT ====
  nbins=10; low=50; high= 450;
  resetPlotMinimum();
  setLogY(false);
  var="jetpt1"; xtitle="p_{T} of lead jet (GeV)";

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "Hjetpt1_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "Hjetpt1_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "Hjetpt1_ge2b");

  // ======= Owen's variables =====
  ratioMin = 0; ratioMax = 3;
  nbins=10; low=50; high= 550;
  resetPlotMinimum();
  setLogY(false);
  var="bestTopMass"; xtitle="3-jet mass (GeV)";

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>=100 && MET<150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "bestTopMass_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>=100 && MET<150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "bestTopMass_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>=100 && MET<150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "bestTopMass_ge2b");

  //signal region
  nbins=20; low=0; high= 800;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "bestTopMass_ge1b_SIG");
  //SB region
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>=100 && MET<150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "bestTopMass_ge1b_wide");
  //LSB region
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET<50";
  drawPlots(var,nbins,low,high,xtitle,"Events", "bestTopMass_ge1b_LSB");

  // ======= study cleaning  =======
  //  bool passBadPFMuon, passInconsistentMuon, passEcalCleaning;
  //no MET cut -- look at BadPFMuon events
  setLogY(true);  resetPlotMinimum();
  nbins=20; low= 0; high=260;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  ratioMin = 0; ratioMax = 2;

  selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passBadPFMuon==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "MET_BadPFMuon");

  selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "MET_InconsistentMuon");

  selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passEcalCleaning==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "MET_FailBE");

  //plot HT (signal region), no b tagging, with and without cleaning
  resetPlotMinimum();
  setLogY(false);
  nbins=10; low=300; high = 1000;
  var="HT"; xtitle="H_{T} (GeV)";
  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passEcalCleaning==0";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HT_FailBE");


   // ==== lead b jet pT ====
   nbins=10; low=30; high= 200;
   resetPlotMinimum();
   setLogY(false);
   var="bjetpt1"; xtitle="p_{T} of lead b jet (GeV)";

   selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
   drawPlots(var,nbins,low,high,xtitle,"Events", "Hbjetpt1_ge1b");

//   selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
//   drawPlots(var,nbins,low,high,xtitle,"Events", "Hjetpt1_eq1b");

//   selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
//   drawPlots(var,nbins,low,high,xtitle,"Events", "Hjetpt1_ge2b");
   
// ==== number of primary vertices ====
   nbins =9; low=1; high=9;
   resetPlotMinimum();
   setLogY(false);
   var="nGoodPV"; xtitle="# of good PVs";
   selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "nGoodPV_SB_ge1b");
   selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>=100 && MET<150 && cutDeltaPhi==0 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "nGoodPV_SBfail_ge1b");




}

void drawMinDeltaPhiMETslices(){

  doOverflowAddition(true);
  setStackMode(true);
  doData(true);
  doRatioPlot(true);
  

  int nbins;
  float low,high;
  TString var,xtitle;
  
  setLogY(false);
  resetPlotMinimum();
  nbins=10; low=0; high = TMath::Pi() + 0.001;
  var="minDeltaPhi"; xtitle="min(#Delta#phi[ jets 1..3, E_{T}^{miss} ] )";
  
  TCut baseSelection ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";
  TCut theseCuts = ""; 

  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";
  for (int ibtag = 0; ibtag<3; ibtag++) { //do this an ugly way for now
    TCut theBTaggingCut = ge1b; TString btagstring = "ge1b";
    if (ibtag==0) { //nothing to do
    }
    else if (ibtag==1) {
      theBTaggingCut = eq1b; 
      btagstring = "eq1b";
    }
    else if (ibtag==2) {
      theBTaggingCut = ge2b; 
      btagstring = "ge2b";
    }
    else assert(0);
    
    //MET < 50
    theseCuts = "MET<50.";
    selection_ = baseSelection && theBTaggingCut && theseCuts;
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+"_MinDeltaPhi_MET_0_50");
    //50 < MET < 100 
    theseCuts = "MET<100. && MET>=50.";
    selection_ = baseSelection && theBTaggingCut && theseCuts;
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+"_MinDeltaPhi_MET_50_100");
    //100 < MET < 150
    theseCuts = "MET<150. && MET>=100.";
    selection_ = baseSelection && theBTaggingCut && theseCuts;
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+"_MinDeltaPhi_MET_100_150");
    //MET < 150
    theseCuts = "MET>=150.";
    selection_ = baseSelection && theBTaggingCut && theseCuts;
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+"_MinDeltaPhi_MET_150_inf");
    
  }
 
}


void drawVJets() {

  doOverflowAddition(true);
  setStackMode(true);
  doData(true);
  doRatioPlot(true);
  setLogY(false);
  
  int nbins;
  float low,high;
  TString var,xtitle;

  TCut baseSelection ="cutHT==1 && cutPV==1 && cutTrigger==1 && passInconsistentMuon==1 && passBadPFMuon==1 && cutDeltaPhi==1 && (nElectrons==0 && nMuons==1)" ;
  TCut njetCut = "njets ==2";  

  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";
  const  TCut pretagCut =  "1";
  for (int ibtag = 0; ibtag<4; ibtag++) { //do this an ugly way for now
    TCut theBTaggingCut = ge1b; TString btagstring = "ge1b";
    if (ibtag==0) { //nothing to do
    }
    else if (ibtag==1) {
      theBTaggingCut = eq1b; 
      btagstring = "eq1b";
    }
    else if (ibtag==2) {
      theBTaggingCut = ge2b; 
      btagstring = "ge2b";
    }
    else if(ibtag==3) {
      theBTaggingCut = pretagCut;
      btagstring = "pre";
    }
    else assert(0);
   
    resetPlotMinimum();
    nbins=10; low=0; high = 10;
    var="njets"; xtitle="njets";
    
    selection_ = baseSelection && theBTaggingCut;
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+"_1muon_njets");
  
    resetPlotMinimum();
    nbins=20; low=10; high = 300;
    var="muonpt1"; xtitle="muonpt1";
 
    selection_ = baseSelection && njetCut && theBTaggingCut;
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+"_1muon_muonpt1");
  }

}

void countABCD() {

  double fitResult[]={52.6, 40.5, 14.1};
  double fitResultErr[]={15.1, 15.9, 8.0};

  loadSamples();

  const  TCut baseSelection = "cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1"; //no MET, no minDeltaPhi, no cleaning, no b tag
  const  TCut passCleaning ="passInconsistentMuon == 1 && passBadPFMuon==1 && weight<1000"; //apply cleaning but without ECAL dead cells
  const  TCut passMinDeltaPhi = "cutDeltaPhi==1";
  const  TCut failMinDeltaPhi = "cutDeltaPhi==0";
  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";

  for (int ibtag = 0; ibtag<3; ibtag++) { //do this an ugly way for now
    TCut theBTaggingCut = ge1b; TString btagstring = "ge1b";
    if (ibtag==0) { //nothing to do
    }
    else if (ibtag==1) {
      theBTaggingCut = eq1b; 
      btagstring = "eq1b";
    }
    else if (ibtag==2) {
      theBTaggingCut = ge2b; 
      btagstring = "ge2b";
    }
    else assert(0);

    TTree* tree= (TTree*) fdata->Get("reducedTree");
    TH1D dummyhist("dummyhist","",1,0,1e9); //kludge!
    dummyhist.Sumw2(); //not really needed for data

    cout<<" ---- "<<btagstring<<" ---- "<<endl;

    TCut METC = "MET>=150 && MET<5000";
    TCut theSIGSelection = baseSelection && passCleaning && passMinDeltaPhi && theBTaggingCut && METC;
    selection_ = theSIGSelection.GetTitle();
    tree->Project("dummyhist","HT",getCutString().Data());
    cout<<"N_SIG = "<<dummyhist.Integral()<<endl;

    dummyhist.Reset();
    if (dummyhist.Integral() != 0) assert(0);

    TCut METD = "MET>=150 && MET<5000";
    TCut theDSelection = baseSelection && passCleaning && failMinDeltaPhi && theBTaggingCut && METD;
    selection_ = theDSelection.GetTitle();
    tree->Project("dummyhist","HT",getCutString().Data());
    double nd = dummyhist.Integral();
    cout<<"N_D = "<<nd<<endl;

    TCut META = "MET>=100 && MET<150";
    TCut theASelection = baseSelection && passCleaning && failMinDeltaPhi && theBTaggingCut && META;
    selection_ = theASelection.GetTitle();
    tree->Project("dummyhist","HT",getCutString().Data());
    double na = dummyhist.Integral();
    cout<<"N_A = "<<na<<endl;

    double R = nd/na;
    double Rerr = R*sqrt(1/nd + 1/na);
    cout<<"R   = "<<R<<" +/- "<<Rerr<<endl;

    double est=R*fitResult[ibtag];
    double estErr = est*sqrt( pow(Rerr/R,2) + pow(fitResultErr[ibtag]/fitResult[ibtag],2));
    cout<<"estimate = "<<est<<" +/- "<<estErr<<endl;

  }
}

void countILV() {
  setQuiet(true);
  loadSamples();
  const  TCut baseSelection = "cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutDeltaPhi==1 && cutMET==1"; //no cleaning, no b tag
  const  TCut passCleaning ="passInconsistentMuon == 1 && passBadPFMuon==1 && weight<1000"; //apply cleaning but without ECAL dead cells

  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";

  const TCut muVeto = "nMuons==0";
  const TCut eVeto = "nElectrons==0";

  const TCut IMV = "nMuons==1";
  const TCut IEV = "nElectrons==1";

  //want to compare SM MC to data
  TCut IMVselection = baseSelection && passCleaning && eVeto && IMV;
  TCut IEVselection = baseSelection && passCleaning && muVeto && IEV;

  TCut IMV_ge1b = IMVselection && ge1b;
  TCut IMV_eq1b = IMVselection && eq1b;
  TCut IMV_ge2b = IMVselection && ge2b;

  TCut IEV_ge1b = IEVselection && ge1b;
  TCut IEV_eq1b = IEVselection && eq1b;
  TCut IEV_ge2b = IEVselection && ge2b;

  doOverflowAddition(false);

  selection_=IMV_ge1b.GetTitle();
  drawPlots("HT",1,0,10000,"","","deleteme");
  cout<<"[IMV ge1b] Integral of data, total SM: "<<hdata->GetBinContent(1)<<" ; "<<totalsm->GetBinContent(1)<< " +/- "<<totalsm->GetBinError(1)<<endl;

  selection_=IMV_eq1b.GetTitle();
  drawPlots("HT",1,0,10000,"","","deleteme");
  cout<<"[IMV eq1b] Integral of data, total SM: "<<hdata->GetBinContent(1)<<" ; "<<totalsm->GetBinContent(1)<< " +/- "<<totalsm->GetBinError(1)<<endl;

  selection_=IMV_ge2b.GetTitle();
  drawPlots("HT",1,0,10000,"","","deleteme");
  cout<<"[IMV ge2b] Integral of data, total SM: "<<hdata->GetBinContent(1)<<" ; "<<totalsm->GetBinContent(1)<< " +/- "<<totalsm->GetBinError(1)<<endl;

  selection_=IEV_ge1b.GetTitle();
  drawPlots("HT",1,0,10000,"","","deleteme");
  cout<<"[IEV ge1b] Integral of data, total SM: "<<hdata->GetBinContent(1)<<" ; "<<totalsm->GetBinContent(1)<< " +/- "<<totalsm->GetBinError(1)<<endl;

  selection_=IEV_eq1b.GetTitle();
  drawPlots("HT",1,0,10000,"","","deleteme");
  cout<<"[IEV eq1b] Integral of data, total SM: "<<hdata->GetBinContent(1)<<" ; "<<totalsm->GetBinContent(1)<< " +/- "<<totalsm->GetBinError(1)<<endl;

  selection_=IEV_ge2b.GetTitle();
  drawPlots("HT",1,0,10000,"","","deleteme");
  cout<<"[IEV ge2b] Integral of data, total SM: "<<hdata->GetBinContent(1)<<" ; "<<totalsm->GetBinContent(1)<< " +/- "<<totalsm->GetBinError(1)<<endl;

}

void drawOwen() {

  /*
.L drawReducedTrees.C++
  */
  
  doOverflowAddition(false);

  const  int nbins = 35;
  const  float min=0;
  const  float max=800;

  bool vb = true;//variable binning ON or OFF
  const int nvarbins=5;
  const float varbins[]={0.,160.,180.,260.,400.,800.};
  //use like this drawSimple("bestTopMass",nvarbins,varbins,histfilename, "bestM3j_met_150_10000_"+btagstring+"tag_qcd","PythiaPUQCDFlat");


  const  TCut baseSelection = "cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1"; //no MET, no minDeltaPhi, no cleaning, no b tag
  // "cutCleaning==1"; //apply all cleaning 
  const  TCut passCleaning ="passInconsistentMuon == 1 && passBadPFMuon==1 && weight<1000"; //apply cleaning but without ECAL dead cells.  Notice the weight cut!
  const  TCut passMinDeltaPhi = "cutDeltaPhi==1";
  const  TCut failMinDeltaPhi = "cutDeltaPhi==0";
  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";

  //
  TString histfilename=  "bestM3j-bins.root";
  TString textfilename=  "drawOwen.output";

  ofstream ofile(textfilename.Data());

  for (int ibtag = 0; ibtag<3; ibtag++) { //do this an ugly way for now
    TCut theBTaggingCut = ge1b; TString btagstring = "ge1b";
    if (ibtag==0) { //nothing to do
    }
    else if (ibtag==1) {
      theBTaggingCut = eq1b; 
      btagstring = "eq1b";
    }
    else if (ibtag==2) {
      theBTaggingCut = ge2b; 
      btagstring = "ge2b";
    }
    else assert(0);
    ofile<<" == "<<btagstring<<" == "<<endl;
 
   //these regions define the PDFs (templates)
    const TCut LSBMET = "MET>=0 && MET<50";
    TCut theLSBSelection = baseSelection && passCleaning && passMinDeltaPhi && theBTaggingCut && LSBMET;
    selection_ = theLSBSelection.GetTitle();
    if(vb){
     drawSimple("bestTopMass",nvarbins,varbins,histfilename, "bestM3j_met_0_50_"+btagstring+"tag_data","data");
     drawPlots("bestTopMass",nvarbins,varbins,"","","deleteme");
    }
    else{
      drawSimple("bestTopMass",nbins,min,max,histfilename, "bestM3j_met_0_50_"+btagstring+"tag_data","data");
      drawPlots("bestTopMass",nbins,min,max,"","","deleteme");
    }
    TFile fh(histfilename,"UPDATE");  
    totalsm->SetName("bestM3j_met_0_50_"+btagstring+"tag_allsm");
    totalsm->Write();          
    fh.Close();   
        
    TCut T2MET = "MET>=50 && MET<5000";
    const  TCut baseT2Selection = "cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && ((nElectrons==0 && nMuons==1) || (nElectrons==1 && nMuons==0))"; //no MET, no minDeltaPhi, no cleaning, no b tag
    TCut theT2Selection = baseT2Selection && passCleaning && passMinDeltaPhi && theBTaggingCut && T2MET;
    selection_ = theT2Selection.GetTitle();
    if (vb) {
      drawSimple("bestTopMass",nvarbins,varbins,histfilename, "bestM3j_met_50_5000_t2_"+btagstring+"tag_data","data");
      drawPlots("bestTopMass",nvarbins,varbins,"","","deleteme");
    }
    else{
      drawSimple("bestTopMass",nbins,min,max,histfilename, "bestM3j_met_50_5000_t2_"+btagstring+"tag_data","data");
      drawPlots("bestTopMass",nbins,min,max,"","","deleteme");
    }
    TFile fh2(histfilename,"UPDATE");  
    totalsm->SetName("bestM3j_met_50_5000_t2_"+btagstring+"tag_allsm");
    totalsm->Write();          
    fh2.Close(); 

    //need T2 in data,ttbar in the SR
    const TCut SRMET = "MET>=150 && MET<5000";
    theT2Selection = baseT2Selection && passCleaning && passMinDeltaPhi && theBTaggingCut && SRMET;
    selection_ = theT2Selection.GetTitle();
    if (vb) {
      drawSimple("bestTopMass",nvarbins,varbins,histfilename, "bestM3j_met_150_5000_t2_"+btagstring+"tag_ttbar","TTbarJets");
      drawSimple("bestTopMass",nvarbins,varbins,histfilename, "bestM3j_met_150_5000_t2_"+btagstring+"tag_data","data");
      drawPlots("bestTopMass",nvarbins,varbins,"","","deleteme");
    }
    else{
      drawSimple("bestTopMass",nbins,min,max,histfilename, "bestM3j_met_150_5000_t2_"+btagstring+"tag_ttbar","TTbarJets");
      drawSimple("bestTopMass",nbins,min,max,histfilename, "bestM3j_met_150_5000_t2_"+btagstring+"tag_data","data");
      drawPlots("bestTopMass",nbins,min,max,"","","deleteme");
    }
    TFile fh22(histfilename,"UPDATE");  
    totalnonttbar->SetName("bestM3j_met_150_5000_t2_"+btagstring+"tag_nonttbar");
    totalewk->SetName("bestM3j_met_150_5000_t2_"+btagstring+"tag_allewk");
    totalnonttbar->Write();          
    totalewk->Write();          
    fh22.Close(); 
    
    // == signal region ==
    TCut theSRSelection = baseSelection && passCleaning && passMinDeltaPhi && theBTaggingCut && SRMET;
    selection_ = theSRSelection.GetTitle();
    if (vb) {
      //plot all samples
      for (unsigned int isample=0; isample<samples_.size(); isample++) {
	TString oname=sampleOwenName_[samples_[isample]];
	drawSimple("bestTopMass",nvarbins,varbins,histfilename, "bestM3j_met_150_5000_"+btagstring+"tag_"+oname,samples_[isample]);
      }
    }
    else {
      //plot all samples
      for (unsigned int isample=0; isample<samples_.size(); isample++) {
	TString oname=sampleOwenName_[samples_[isample]];
	drawSimple("bestTopMass",nbins,min,max,histfilename, "bestM3j_met_150_5000_"+btagstring+"tag_"+oname,samples_[isample]);
      }
    }

    //need invdphi in data in the SR
    TCut invdpSelection = baseSelection && passCleaning && failMinDeltaPhi && theBTaggingCut && SRMET;
    TString nameOfIDPhist = "bestM3j_met_150_5000_invdphi_";
    nameOfIDPhist+=btagstring; nameOfIDPhist+="tag_";
    selection_ = invdpSelection.GetTitle();
    if (vb) {
      drawSimple("bestTopMass",nvarbins,varbins,histfilename, nameOfIDPhist+"qcd" ,"PythiaPUQCD");
      drawSimple("bestTopMass",nvarbins,varbins,histfilename, nameOfIDPhist+"data" ,"data");
      drawPlots("bestTopMass",nvarbins,varbins,"","","deleteme");
    }
    else {
      drawSimple("bestTopMass",nbins,min,max,histfilename, nameOfIDPhist+"qcd" ,"PythiaPUQCD");
      drawSimple("bestTopMass",nbins,min,max,histfilename, nameOfIDPhist+"data" ,"data");
      drawPlots("bestTopMass",nbins,min,max,"","","deleteme");
    }
    TFile fh3(histfilename,"UPDATE");
    totalnonqcd->SetName(nameOfIDPhist+"nonqcd");
    totalnonqcd->Write();
    fh3.Close();

    //now for a flexible MET region
    for (int metCutLow = 80; metCutLow <=100; metCutLow+=10) {
      for (int metCutHigh = 150; metCutHigh <=150; metCutHigh+=5) {
	TString metCutString; metCutString.Form("MET >= %d && MET < %d",metCutLow,metCutHigh);
	ofile<<metCutString<<endl;
	TCut METselection(metCutString.Data());

	theT2Selection = baseT2Selection && passCleaning && passMinDeltaPhi && theBTaggingCut && METselection;
	TString nameOfT2Hist;
	nameOfT2Hist.Form( "bestM3j_met_%d_%d_t2_%stag_",metCutLow,metCutHigh,btagstring.Data());

	TCut theSelection = baseSelection && passCleaning && passMinDeltaPhi && theBTaggingCut && METselection;
	selection_ = theSelection.GetTitle();
	TString nameOfHist;
	nameOfHist.Form( "bestM3j_met_%d_%d_%stag_",metCutLow,metCutHigh,btagstring.Data());
	
	//*** nominal selection ***//
	if (vb) {
	  //plot all samples
	  for (unsigned int isample=0; isample<samples_.size(); isample++) {
	    TString oname=sampleOwenName_[samples_[isample]];
	    ofile<<oname<<" = "<<drawSimple("bestTopMass",nvarbins,varbins,histfilename, nameOfHist+oname,samples_[isample])<<endl;
	  }
	  //plot data
	  drawSimple("bestTopMass",nvarbins,varbins,histfilename, nameOfHist+"data","data");
	  //fills plots that are combinations of various samples (to be accessed via global pointers)
	  drawPlots("bestTopMass",nvarbins,varbins,"","","deleteme");
	}
	else {
	  //plot all samples
	  for (unsigned int isample=0; isample<samples_.size(); isample++) {
	    TString oname=sampleOwenName_[samples_[isample]];
	    ofile<<oname<<" = "<<drawSimple("bestTopMass",nbins,min,max,histfilename, nameOfHist+oname,samples_[isample])<<endl;
	  }
	  //plot data
	  drawSimple("bestTopMass",nbins,min,max,histfilename, nameOfHist+"data","data");
	  //fills plots that are combinations of various samples (to be accessed via global pointers)
	  drawPlots("bestTopMass",nbins,min,max,"","","deleteme");
	}
	TFile fh4(histfilename,"UPDATE");
	totalnonttbar->SetName(nameOfHist+"nonttbar");
	totalnonqcd->SetName(nameOfHist+"nonqcd");
	totalewk->SetName(nameOfHist+"allewk");
	totalsm->SetName(nameOfHist+"allsm");
	totalnonttbar->Write();
	totalnonqcd->Write();
	totalewk->Write();
	totalsm->Write();
	fh4.Close();
	
	//*** T2 selection ***//
	selection_ = theT2Selection.GetTitle(); //change to t2 selection
	if (vb) {
	  drawSimple("bestTopMass",nvarbins,varbins,histfilename, nameOfT2Hist+"data","data");
	  drawSimple("bestTopMass",nvarbins,varbins,histfilename, nameOfT2Hist+"ttbar","TTbarJets");
	  drawPlots("bestTopMass",nvarbins,varbins,"","","deleteme");
	}
	else {
	  drawSimple("bestTopMass",nbins,min,max,histfilename, nameOfT2Hist+"data","data");
	  drawSimple("bestTopMass",nbins,min,max,histfilename, nameOfT2Hist+"ttbar","TTbarJets");
	  drawPlots("bestTopMass",nbins,min,max,"","","deleteme");
	}
	TFile fh44(histfilename,"UPDATE");  
	totalnonttbar->SetName(nameOfT2Hist+"nonttbar");
	totalewk->SetName(nameOfT2Hist+"allewk");
	totalnonttbar->Write();          
	totalewk->Write();          
	fh44.Close(); 

	//*** fail mdp selection ***//
      	theSelection = baseSelection && passCleaning && failMinDeltaPhi && theBTaggingCut && METselection;
	selection_ = theSelection.GetTitle();
	nameOfHist.Form( "bestM3j_met_%d_%d_invdphi_%stag_",metCutLow,metCutHigh,btagstring.Data());
	if(vb){
	  drawSimple("bestTopMass",nvarbins,varbins,histfilename, nameOfHist+"data","data");
	  drawSimple("bestTopMass",nvarbins,varbins,histfilename, nameOfHist+"qcd","PythiaPUQCD");
	  drawPlots("bestTopMass",nvarbins,varbins,"","","deleteme");
	}
	else{
	  drawSimple("bestTopMass",nbins,min,max,histfilename, nameOfHist+"data","data");
	  drawSimple("bestTopMass",nbins,min,max,histfilename, nameOfHist+"qcd","PythiaPUQCD");
	  drawPlots("bestTopMass",nbins,min,max,"","","deleteme");
	}
	TFile fh5(histfilename,"UPDATE");
	totalnonqcd->SetName(nameOfHist+"nonqcd");
	totalsm->SetName(nameOfHist+"allsm");
	totalnonqcd->Write();
	totalsm->Write();
	fh5.Close();
      }
    }
     
  }
  ofile.close();

}
