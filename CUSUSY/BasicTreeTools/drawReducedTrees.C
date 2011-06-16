/*
====== this is the nominal code for drawing RA2b data/MC comparison plots =======
-- drawPlots() -- main plotting routine to draw pretty stacks of MC with data on top
-- drawSimple() -- very simple function to draw exactly one sample and put it in a file
-- drawR() -- draws r(MET). a bit more kludgely, but works.
-- drawSignificance -- draws S/sqrt(B) as a function of the variable

The functions that do the heavy lifting, along with various utility functions,
are in drawReducedTrees.h.

Samples to plot are currently hard-coded in loadSamples().
New functions addSample() and removeSample() add some run-time flexibility.

Various routines here call the above functions. Some of them (e.g. drawSomething() ) I usually
use in a quasi-interactive mode.

-- drawOwen() -- uses drawPlots() and drawSimple() to make a file with the histograms needed for
toys and fits to data.

Details:
This code replaces drawCutflowPlots.C (which had replaced drawBasicPlots.C). 
The input files are the reducedTrees.

Depending on the exact setup in loadSamples(), the QCD and Single-top samples must be added together with 'hadd' 
in order to get one file per sample.

Potential improvements:
 -- there are becoming way too many configuration options. i've stopped adding setter functions for them out of laziness
 -- the calls to renormBins() currently have a kludge that hard-codes the reference bin to 2
 -- The samples to be plotted are currently defined at compile time.
Could make it so that all samples are always loaded, but there is an
independent list that controls which samples are plotted. This would
allow the user to switch the plotted samples on the fly.
[any mechanism like this would also have to allow the user to control
the order in which the samples are stack. This is currently defined
by the order of the push_backs in loadSamples()]

---- update: this is now partially fixed. samples_ can be manipulated after compile time.
More functions for manipulation (e.g. to allow adding a sample in the middle of the list) still
need to be written

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

// can be checked out of UserCode/joshmt
#include "MiscUtil.cxx"

#include <fstream>

#include <iostream>
#include <map>
#include <set>


TString inputPath = "/cu2/kreis/reducedTrees/V00-01-02/";
//TString cutdesc = "Baseline0_PF_JERbias_pfMEThigh_PFLep0e0mu_minDP_MuonEcalCleaning";
TString cutdesc = "TCHET";
//TString cutdesc = "Baseline0_PF_pfMEThigh_PFLepRA20e0mu_minDP_MuonEcalCleaning";
//TString cutdesc = "Baseline0_PF_pfMEThigh_PFLep0e0mu_minDP_MuonEcalCleaning";

#include "drawReducedTrees.h"

void morerstuff() {
  /*
.L drawReducedTrees.C++
  */
  setStackMode(false);
  doData(false);

  useFlavorHistoryWeights_=false;

  loadSamples(); //needs to come first! this is why this should be turned into a class, with this guy in the ctor
  clearSamples();
  addSample("PythiaPUQCD");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";
  setLogY(true);
  drawR("minDeltaPhi",0.3,13,70,230,"ge1b");
  setLogY(false);
  drawR("minDeltaPhi",0.3,13,70,230,"ge1b");

  //  TF1 *fa = new TF1("fa","[0]*exp(-1*[1]*x)",-3,3); 

}


void biasStudy() {
  /*
.L drawReducedTrees.C++
  */
  useFlavorHistoryWeights_=false;

  setStackMode(false);
  doData(true);

  //no met, no mindeltaphi
  setPlotMinimum(1e-2); //setPlotMaximum(1);

 drawTotalSM_=true;
 leg_y1=0.6;
 setLogY(true);
 
 setQuiet(true);

 double cbias,err;

 selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000 &&njets==3";
 drawR("minDeltaPhi",0.3,4,0,200,"ge1b-4bins");
 cbias= hinteractive->GetBinContent(4)/hinteractive->GetBinContent(3);
 err = jmt::errAoverB(hinteractive->GetBinContent(4),hinteractive->GetBinError(4),hinteractive->GetBinContent(3),hinteractive->GetBinError(3));
 cout<<"no extra cut cbias = "<<cbias <<" +/- "<<err<<endl;

 selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000 && (MET/HT<0.2) &&njets==3";
 drawR("minDeltaPhi",0.3,4,0,200,"ge1b-METoverHTcut-4bins");
 cbias= hinteractive->GetBinContent(4)/hinteractive->GetBinContent(3);
 err = jmt::errAoverB(hinteractive->GetBinContent(4),hinteractive->GetBinError(4),hinteractive->GetBinContent(3),hinteractive->GetBinError(3));
 cout<<"MET/HT<0.2   cbias = "<<cbias <<" +/- "<<err<<endl;

 selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1  && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000 && (MET/HT>=0.2) &&njets==3";
 drawR("minDeltaPhi",0.3,4,0,200,"ge1b-METoverHTcut-4bins");
 cbias= hinteractive->GetBinContent(4)/hinteractive->GetBinContent(3);
 err = jmt::errAoverB(hinteractive->GetBinContent(4),hinteractive->GetBinError(4),hinteractive->GetBinContent(3),hinteractive->GetBinError(3));
 cout<<"MET/HT>=0.2  cbias = "<<cbias <<" +/- "<<err<<endl;


}

void drawrplots() {
  /*
.L drawReducedTrees.C++
  */

  useFlavorHistoryWeights_=false;

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


void drawSoftJets() {
/*
.L drawReducedTrees.C++
*/

//this seems completely uninteresting!

  setStackMode(true);
  doData(true);
  int nbins;
  float low,high;
  TString var,xtitle;
  doRatioPlot(false);
  setLogY(false);

  nbins=10; low=0; high=10;
  var="nLooseJets20_30"; xtitle="Soft Jet multiplicity";
  //ratioMin = 0; ratioMax = 2;
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 &&cutMET==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "nsoftjets_ge1b");
  

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

  // TO DO add weight<1000 cut everywhere!


  // ==== MET plots ====
  setLogY(true);   setPlotMinimum(1e-1);
  nbins=25; low= 0; high=260;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  ratioMin = 0; ratioMax = 2;

  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET");
  //  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_JER0"); //in case we're doing the non-JERbias version

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
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

  //check bad PV events
   nbins =6; low=0; high=5;
   resetPlotMinimum();
   setLogY(false);
   var="nGoodPV"; xtitle="# of good PVs";
   selection_ ="nbjets>=1 && cutHT==1 && cutPV==0 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&MET>50";
   drawPlots(var,nbins,low,high,xtitle,"Events", "nGoodPV_PVfail_ge1b");
   selection_ ="cutHT==1 && cutPV==0 && cutTrigger==1 && cutDeltaPhi==1";
   drawPlots(var,nbins,low,high,xtitle,"Events", "nGoodPV_PVfail_veryloose");


   //exploration
   setLogY(false);
   resetPlotMinimum();
 
   nbins=20; low=-TMath::Pi(); high = TMath::Pi();
   var="METphi-jetphi1"; xtitle="(METphi-jetphi1)";
   selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>100 && passInconsistentMuon==1 && passBadPFMuon==1";
   drawPlots(var,nbins,low,high,xtitle,"Events", "HdeltaPhiMETjet1_nomindp");

   nbins=20; low=-TMath::Pi(); high = TMath::Pi();
   var="METphi-jetphi2"; xtitle="METphi-jetphi2";
   selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>100 && passInconsistentMuon==1 && passBadPFMuon==1 && (abs(METphi-jetphi1)>1)";
   drawPlots(var,nbins,low,high,xtitle,"Events", "HdeltaPhiMETjet2_metphi1cut");

   var="METphi-jetphi2"; xtitle="METphi-jetphi2";
   selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>150 && passInconsistentMuon==1 && passBadPFMuon==1 && (abs(METphi-jetphi1)>1)";
   drawPlots(var,nbins,low,high,xtitle,"Events", "HdeltaPhiMETjet2_metphi1cut_sig");

   var="METphi-jetphi3"; xtitle="METphi-jetphi3";
   selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && MET>150 && passInconsistentMuon==1 && passBadPFMuon==1 && (abs(METphi-jetphi1)>1) && (abs(METphi-jetphi2)>0.3)";
   drawPlots(var,nbins,low,high,xtitle,"Events", "HdeltaPhiMETjet3_metphi12cut_sig");

  // ==== MET/HT ====
   //addSample("LM13");
  doRatioPlot(false);
  setLogY(true);
  nbins=40;
  low=0;
  high=1;
  doOverflowAddition(true);
  var="MET/HT"; xtitle="MET/HT";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";
  drawPlots(var,nbins,low,high,xtitle,"Events", "metOverHt_ge1b_noMetCut");

  setLogY(false);
  resetPlotMinimum();
  nbins=10;

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 && MET>100 &&MET<150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "metOverHt_ge1b_SB");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&MET>=150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "metOverHt_ge1b_SIG");

  // ==== MET/sqrt(HT) ====
   //addSample("LM13");
  doRatioPlot(false);
  nbins=20;  low=0;  high=20;
  doOverflowAddition(true);
  var="MET/sqrt(HT)"; xtitle="MET/#sqrt{HT}";

   setLogY(false);
  resetPlotMinimum();

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 && MET>100 &&MET<150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "metOverRootHt_ge1b_SB");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&MET>=150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "metOverRootHt_ge1b_SIG");


  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&MET>=100";
  drawPlots(var,nbins,low,high,xtitle,"Events", "metOverRootHt_ge1b_metOver100");

}

void studySignificance() {
/*
.L drawReducedTrees.C++
  */

  setStackMode(true);
  doData(true);

  int nbins;
  float low,high;
  TString var,xtitle;
  doOverflowAddition(true);
  doRatioPlot(false);


  setLogY(false);
  resetPlotMinimum();

  //this compares the effectiveness of MET and MET/sqrt(HT); MET is better
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000";

  nbins=20;  low=0;  high=20;
  var="MET/sqrt(HT)"; xtitle="MET/#sqrt{HT}";
  drawSignificance(var,nbins,low,high,"sigfile");

  nbins=30;  low=0;  high=300;
  var="MET"; xtitle="MET";
  drawSignificance(var,nbins,low,high,"sigfile");


  //after MET cut, what can we do?
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&cutMET==1";
  nbins=20;  low=0;  high=20;
  var="MET/sqrt(HT)"; xtitle="MET/#sqrt{HT}";
  drawSignificance(var,nbins,low,high,"sigfile"); //nothing good here

  //reoptimize minDeltaPhi
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&cutMET==1";
  nbins=31;  low=0;  high=TMath::Pi();
  var="minDeltaPhi"; xtitle="minDeltaPhi";
  drawSignificance(var,nbins,low,high,"sigfile");

  var="minDeltaPhiAll30"; xtitle="minDeltaPhi All 30"; //not as good!
  drawSignificance(var,nbins,low,high,"sigfile");


  //remove 3Jet cut
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&cutMET==1";
  nbins=7; low=0; high=7;
  var="njets"; xtitle="jet multiplicity";
  drawSignificance(var,nbins,low,high,"sigfile");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&cutMET==1";
  nbins=40; low=200; high=1000;
  var="HT"; xtitle="HT";
  drawSignificance(var,nbins,low,high,"sigfile");

  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&cutMET==1";
  nbins=4; low=0; high=4;
  var="nbjets"; xtitle="n b tags";
  drawSignificance(var,nbins,low,high,"sigfile");

  //nominal cuts
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 &&weight<1000 &&cutMET==1";
  nbins=40;  low=0;  high=400;
  var="jetpt1"; xtitle="lead jet pt";
  drawSignificance(var,nbins,low,high,"sigfile");

  nbins=40;  low=0;  high=400;
  var="bjetpt1"; xtitle="lead b jet pt";
  drawSignificance(var,nbins,low,high,"sigfile");

}

void drawMinDeltaPhiMETslices(){
  useFlavorHistoryWeights_=false;
  doOverflowAddition(true);
  doData(true);
  doRatioPlot(true);
  setQuiet(false);
  ratioMin=0; ratioMax=3;
  
  bool doVB = false;
  bool doSlices = true;
  if(!doSlices){
    leg_y1=0.6; //and also remember to use correct color scheme
    owenColor_=true;
  }  

  int nbins;
  float low,high;
  TString var,xtitle;

  nbins=10; low=0; high = TMath::Pi() + 0.001;
  var="minDeltaPhi"; xtitle="min(#Delta#phi[ jets 1..3, E_{T}^{miss} ] )";
  const int nvarbins=2;
  const float varbins[]={0.,0.3,high};
  
  TCut baseSelection =""; 
  TCut theseCuts = ""; 
  TString extraName = "";

  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";
  const  TCut pretag =  "1";
  const  TCut antitag = "nbjets == 0";
  for (int ibtag = 0; ibtag<5; ibtag++) { 
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
    else if (ibtag==3) {
      theBTaggingCut = pretag;
      btagstring = "pre";
    }
    else if (ibtag==4) {
      theBTaggingCut = antitag;
      btagstring = "antib";
    }
    else assert(0);
    

    //+++++++++++NOMINAL SELECTION+++++++++++//
    baseSelection = "cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passInconsistentMuon==1 && passBadPFMuon==1";//nominal
    extraName = "_nominal";
     
    if(doSlices){
      //MinDeltaPhi in MET slices
      setStackMode(true); setLogY(false); resetPlotMinimum();
      //MET < 50
      theseCuts = "MET<50.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_0_50");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_0_50_vb");
      //50 < MET < 100 
      theseCuts = "MET<100. && MET>=50.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_50_100");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_50_100_vb");
      //100 < MET < 150
      theseCuts = "MET<150. && MET>=100.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_100_150");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_100_150_vb");
      //MET < 150
      theseCuts = "MET>=150.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_150_inf"); 
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_150_inf_vb");
    }
    else{
      //Ratio plot
      selection_ = baseSelection && theBTaggingCut;
      setStackMode(false); setPlotMinimum(1e-2); drawTotalSM_=true; setLogY(true); 
      drawR("minDeltaPhi",0.3,20,0,200,btagstring+extraName);
      drawR("minDeltaPhi",0.3,4,0,200,btagstring+extraName+"_4bins");
    }
    
    /*
    //[+++++++++++++njets==2++++++++]//
    baseSelection = "cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && passInconsistentMuon==1 && passBadPFMuon==1 && njets==2"; //remove cut3Jets and add njets
    extraName = "_njetsEQ2";

    if(doSlices){
      //MinDeltaPhi in MET slices
      setStackMode(true); setLogY(false); resetPlotMinimum();
      //MET < 50
      theseCuts = "MET<50.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_0_50");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_0_50_vb");
      //50 < MET < 100 
      theseCuts = "MET<100. && MET>=50.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_50_100");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_50_100_vb");
      //100 < MET < 150
      theseCuts = "MET<150. && MET>=100.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_100_150");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_100_150_vb");
      //MET < 150
      theseCuts = "MET>=150.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_150_inf");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_150_inf_vb");
    }
    else{
      //Ratio plot
      selection_ = baseSelection && theBTaggingCut;
      setStackMode(false); setPlotMinimum(1e-2); drawTotalSM_=true; setLogY(true); 
      drawR("minDeltaPhi",0.3,20,0,200,btagstring+extraName);
      drawR("minDeltaPhi",0.3,4,0,200,btagstring+extraName+"_4bins");
    }
    */
    
    /*
    //[+++++++++++m3j>350++++++++++++]//
    baseSelection = "cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passInconsistentMuon==1 && passBadPFMuon==1 && bestTopMass>350";//nominal plus m3j cut
    extraName = "_m3jGT350";
    
    if(doSlices){
      //MinDeltaPhi in MET slices
      setStackMode(true); setLogY(false); resetPlotMinimum();
      //MET < 50
      theseCuts = "MET<50.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_0_50");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_0_50_vb");
      //50 < MET < 100 
      theseCuts = "MET<100. && MET>=50.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_50_100");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_50_100_vb");
      //100 < MET < 150
      theseCuts = "MET<150. && MET>=100.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_100_150");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_100_150_vb");
      //MET < 150
      theseCuts = "MET>=150.";
      selection_ = baseSelection && theBTaggingCut && theseCuts;
      drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_150_inf");
      if(doVB) drawPlots(var,nvarbins,varbins,xtitle,"Events", btagstring+extraName+"_MinDeltaPhi_MET_150_inf_vb");
    }
    else{
      //Ratio plot
      selection_ = baseSelection && theBTaggingCut;
      setStackMode(false); setPlotMinimum(1e-2); drawTotalSM_=true; setLogY(true); 
      drawR("minDeltaPhi",0.3,20,0,200,btagstring+extraName);
      drawR("minDeltaPhi",0.3,4,0,200,btagstring+extraName+"_4bins");
    }
    */

  }//b-tag cut loop 
}


void drawMETPlots() {
  useFlavorHistoryWeights_=false;
  setStackMode(true);
  doData(true);
  doRatioPlot(true);
  setQuiet(false);

  int nbins;
  float low,high;
  TString var,xtitle;

  // ==== MET plots ====
  setLogY(true);   setPlotMinimum(1e-1);
  nbins=25; low= 0; high=500;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  ratioMin = 0; ratioMax = 3;

  /*
  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_ge2b");

  selection_ ="nbjets==0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET_antib");
  */

  setLogY(false);
  setPlotMaximum(100);
  selection_ ="nbjets==0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>=150.";
  resetPlotMaximum();
  drawPlots(var,10,150,high,xtitle,"Events", "H_MET_antib");
}


void drawVJets() {

  useFlavorHistoryWeights_=true;

  doOverflowAddition(true);
  setStackMode(true);
  doData(true);
  doRatioPlot(true);
  setLogY(false);
  
  setQuiet(true);

  int nbins;
  float low,high;
  TString var,xtitle;

  TCut baseSelection ="cutHT==1 && cutPV==1 && cutTrigger==1 && passInconsistentMuon==1 && passBadPFMuon==1 && cutDeltaPhi==1 && (nElectrons==0 && nMuons==1)";// && bestTopMass>350" ;
  TString extraName = ""; //"_m3jGT350";
  TCut njetCut = "njets ==2";  

  const TCut METcut = "MET>40";
  const TCut METcut2 = "MET>60";
  const TCut METcut3 = "MET>100";
  const TCut MTcut = "MT_Wlep>30";

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
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_1muon_njets");


    selection_ = baseSelection && njetCut && theBTaggingCut;  

    resetPlotMinimum();
    nbins=10; low=0; high = 200;
    var="MET"; xtitle="MET (GeV)"; 
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_1muon2jet_met");
   
    resetPlotMinimum();
    nbins=20; low=10; high = 300;
    var="muonpt1"; xtitle="muonpt1"; 
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_1muon2jet_muonpt1");
    
    { //give this its own scope
      cout<<"no extra cuts"<<endl;
      double ndata= hdata->Integral();
      double nnonW=0,nW=0;
      double nnonW_err=0,nW_err=0;
      for ( std::map<TString, TH1D*>::iterator ihist = histos_.begin(); ihist!=histos_.end(); ++ihist) {
	if (ihist->second != 0) {
	  if (TString(ihist->second->GetName()).Contains("WJets")) {
	    nW = ihist->second->Integral();
	    nW_err = jmt::errOnIntegral(ihist->second);
	  }
	  else {
	    nnonW += ihist->second->Integral();
	    nnonW_err += pow(jmt::errOnIntegral(ihist->second),2);
	  }
	}
      }
      cout<<"(data - nonW)/W = ("<<ndata<<" - "<<nnonW<<") / "<<nW<<endl;
      double err = jmt::errAoverB( ndata-nnonW, sqrt(ndata+nnonW_err), nW,nW_err);
      cout<<"                = "<<(ndata-nnonW)/nW<<" +/- "<<err<<endl;
    }

    selection_ = baseSelection && njetCut && theBTaggingCut && METcut;  
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_1muon2jet_muonpt1_MET40");
    { //give this its own scope
      cout<<"MET>40"<<endl;
      double ndata= hdata->Integral();
      double nnonW=0,nW=0;
      double nnonW_err=0,nW_err=0;
      for ( std::map<TString, TH1D*>::iterator ihist = histos_.begin(); ihist!=histos_.end(); ++ihist) {
	if (ihist->second != 0) {
	  if (TString(ihist->second->GetName()).Contains("WJets")) {
	    nW = ihist->second->Integral();
	    nW_err = jmt::errOnIntegral(ihist->second);
	  }
	  else {
	    nnonW += ihist->second->Integral();
	    nnonW_err += pow(jmt::errOnIntegral(ihist->second),2);
	  }
	}
      }
      cout<<"(data - nonW)/W = ("<<ndata<<" - "<<nnonW<<") / "<<nW<<endl;
      double err = jmt::errAoverB( ndata-nnonW, sqrt(ndata+nnonW_err), nW,nW_err);
      cout<<"                = "<<(ndata-nnonW)/nW<<" +/- "<<err<<endl;
    }

    selection_ = baseSelection && njetCut && theBTaggingCut && METcut2;  
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_1muon2jet_muonpt1_MET60");
    { //give this its own scope
      cout<<"MET>60"<<endl;
      double ndata= hdata->Integral();
      double nnonW=0,nW=0;
      double nnonW_err=0,nW_err=0;
      for ( std::map<TString, TH1D*>::iterator ihist = histos_.begin(); ihist!=histos_.end(); ++ihist) {
	if (ihist->second != 0) {
	  if (TString(ihist->second->GetName()).Contains("WJets")) {
	    nW = ihist->second->Integral();
	    nW_err = jmt::errOnIntegral(ihist->second);
	  }
	  else {
	    nnonW += ihist->second->Integral();
	    nnonW_err += pow(jmt::errOnIntegral(ihist->second),2);
	  }
	}
      }
      cout<<"(data - nonW)/W = ("<<ndata<<" - "<<nnonW<<") / "<<nW<<endl;
      double err = jmt::errAoverB( ndata-nnonW, sqrt(ndata+nnonW_err), nW,nW_err);
      cout<<"                = "<<(ndata-nnonW)/nW<<" +/- "<<err<<endl;
    }

    selection_ = baseSelection && njetCut && theBTaggingCut && METcut3;  
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_1muon2jet_muonpt1_MET100");
    { //give this its own scope
      cout<<"MET>100"<<endl;
      double ndata= hdata->Integral();
      double nnonW=0,nW=0;
      double nnonW_err=0,nW_err=0;
      for ( std::map<TString, TH1D*>::iterator ihist = histos_.begin(); ihist!=histos_.end(); ++ihist) {
	if (ihist->second != 0) {
	  if (TString(ihist->second->GetName()).Contains("WJets")) {
	    nW = ihist->second->Integral();
	    nW_err = jmt::errOnIntegral(ihist->second);
	  }
	  else {
	    nnonW += ihist->second->Integral();
	    nnonW_err += pow(jmt::errOnIntegral(ihist->second),2);
	  }
	}
      }
      cout<<"(data - nonW)/W = ("<<ndata<<" - "<<nnonW<<") / "<<nW<<endl;
      double err = jmt::errAoverB( ndata-nnonW, sqrt(ndata+nnonW_err), nW,nW_err);
      cout<<"                = "<<(ndata-nnonW)/nW<<" +/- "<<err<<endl;
    }


    selection_ = baseSelection && njetCut && theBTaggingCut && MTcut;  
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_1muon2jet_muonpt1_MT15");
    { //give this its own scope
      cout<<"MT>30"<<endl;
      double ndata= hdata->Integral();
      double nnonW=0,nW=0;
      double nnonW_err=0,nW_err=0;
      for ( std::map<TString, TH1D*>::iterator ihist = histos_.begin(); ihist!=histos_.end(); ++ihist) {
	if (ihist->second != 0) {
	  if (TString(ihist->second->GetName()).Contains("WJets")) {
	    nW = ihist->second->Integral();
	    nW_err = jmt::errOnIntegral(ihist->second);
	  }
	  else {
	    nnonW += ihist->second->Integral();
	    nnonW_err += pow(jmt::errOnIntegral(ihist->second),2);
	  }
	}
      }
      cout<<"(data - nonW)/W = ("<<ndata<<" - "<<nnonW<<") / "<<nW<<endl;
      double err = jmt::errAoverB( ndata-nnonW, sqrt(ndata+nnonW_err), nW,nW_err);
      cout<<"                = "<<(ndata-nnonW)/nW<<" +/- "<<err<<endl;
      cout<<"W / SM = "<<nW/(nW+nnonW)<<endl;
    }


    selection_ = baseSelection && njetCut && theBTaggingCut;  
    nbins=20; low=0; high=300;
    //nbins=25; low=-5; high=10;
    var="MT_Wlep"; xtitle="leptonic W M_{T}";
    drawPlots(var,nbins,low,high,xtitle,"Events", btagstring+extraName+"_1muon2jet_MT");

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
    tree->Project("dummyhist","HT",getCutString(1.,"",selection_, "",0).Data());
    cout<<"N_SIG = "<<dummyhist.Integral()<<endl;

    dummyhist.Reset();
    if (dummyhist.Integral() != 0) assert(0);

    TCut METD = "MET>=150 && MET<5000";
    TCut theDSelection = baseSelection && passCleaning && failMinDeltaPhi && theBTaggingCut && METD;
    selection_ = theDSelection.GetTitle();
    tree->Project("dummyhist","HT",getCutString(1.,"",selection_,"",0).Data());
    double nd = dummyhist.Integral();
    cout<<"N_D = "<<nd<<endl;

    TCut META = "MET>=100 && MET<150";
    TCut theASelection = baseSelection && passCleaning && failMinDeltaPhi && theBTaggingCut && META;
    selection_ = theASelection.GetTitle();
    tree->Project("dummyhist","HT",getCutString(1.,"",selection_,"",0).Data());
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

void flvHistReweighting() {

  //we want to make a way to apply the k-factors 

  //can use the extraSelection part of the getCutString()

  dodata_=false;
  const TString samplename = "WJets"; //if this is changed to data, lumiScale_ should not be applied
  setQuiet(true);
  loadSamples();
  
  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";

  TTree* tree=0;
  if (samplename=="data") {
    tree = (TTree*) fdata->Get("reducedTree");
  }
  else {
    for (unsigned int isample=0; isample<samples_.size(); isample++) {
      if ( samples_[isample].Contains(samplename)) {
	tree = (TTree*) files_[samples_[isample]]->Get("reducedTree");
      }
    }
  }
  if (tree==0) {cout<<"Something went wrong finding your sample!"<<endl; return;}

  gROOT->cd();

  TH1D dummyhist("dummyhist","",1,0,1e9); //the typical kludge to count events; 1 bin histogram with large range
  dummyhist.Sumw2();

  //nominal case
  tree->Project("dummyhist","HT",getCutString(lumiScale_,"",selection_,"",0).Data());
  float n_nom = dummyhist.GetBinContent(1);
  float nerr_nom = dummyhist.GetBinError(1);
  cout<<" nominal     = "<<n_nom<<" +/- "<<nerr_nom<<endl;
  dummyhist.Reset();

  //with flv hist
  tree->Project("dummyhist","HT",getCutString(lumiScale_,"flavorHistoryWeight",selection_,"",0).Data());
  float  n_fh = dummyhist.GetBinContent(1);
  float  nerr_fh = dummyhist.GetBinError(1);
  cout<<" full reweight    = "<<n_fh<<" +/- "<<nerr_fh<<endl;

  const int ncat=11;
  const float kfactors[]  = {2,2,1,1,2,1,2,1,2,1,1};
  const float kerr_up[]   = {1,1,1,1,1,1,1,1,1,1,0};
  const float kerr_down[] = {1,1,0.5,0.5,1,0.5,1,0.5,1,0.5,0};

  float total_up=0;
  float total_down=0;

  //loop over the variations and up/down
  for (int iupdown = 0; iupdown<=1; iupdown++ ) {
    for (int ivariations = 0; ivariations<ncat; ivariations++) {
      dummyhist.Reset();
      TString manualreweight="";
      for (int ii=0; ii<ncat; ii++) {
	float k = kfactors[ii];
	//vary only one at a time
	if (ivariations == ii) {
	  if (iupdown==0) k += kerr_up[ii];
	  else            k -= kerr_down[ii];
	}
	TString thisterm;
	thisterm.Form( "((flavorHistory==%d) * %f)", ii+1, k);
	manualreweight += thisterm;
	if (ii+1 != ncat) manualreweight += " + ";
      }
      tree->Project("dummyhist","HT",getCutString(lumiScale_,manualreweight,selection_,"",0).Data()); //this will be normalized *wrong* because it doesn't have the N/sum(k_i) factor
      float  n = dummyhist.GetBinContent(1);
      float  nerr = dummyhist.GetBinError(1);
      dummyhist.Reset();
      assert(0);//please double check that the following few lines of code have the weights and cuts handled correctly
      tree->Project("dummyhist","HT","1"); //weight of 1
      double unweighted = dummyhist.GetBinContent(1);
      dummyhist.Reset();
      tree->Project("dummyhist","HT",manualreweight); //weighted with just the k factors
      double reweighted = dummyhist.GetBinContent(1);

      n *= unweighted/reweighted;
      nerr *= unweighted/reweighted;

      TString uddesc= iupdown==0 ? " up " : "down";
      cout<<"reweight "<< uddesc<<" [" << ivariations+1<<"] = "<<n<<" +/- "<<nerr<<" ; "<<n-n_fh<<endl;
      if (iupdown==0) total_up += pow(n-n_fh,2);
      else total_down += pow(n-n_fh,2);
    }
  }

  cout<<" up  syst = "<<sqrt(total_up)<<endl;
  cout<<"down syst = "<<sqrt(total_down)<<endl;

}

void pdfUncertainties() {
  dodata_=false;
  const TString samplename = "WJetsZ2"; //if this is changed to data, lumiScale_ should not be applied
  setQuiet(true);
  loadSamples();

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && weight<1000";

  TTree* tree=0;
  if (samplename=="data") {
    tree = (TTree*) fdata->Get("reducedTree");
  }
  else {
    for (unsigned int isample=0; isample<samples_.size(); isample++) {
      if ( samples_[isample] == samplename) {
	tree = (TTree*) files_[samples_[isample]]->Get("reducedTree");
      }
    }
  }
  if (tree==0) {cout<<"Something went wrong finding your sample!"<<endl; return;}

  gROOT->cd();
  TString optfh= useFlavorHistoryWeights_ && samplename.Contains("WJets") ? "flavorHistoryWeight" : "";

  TH1D dummyhist("dummyhist","",1,0,1e9); //the typical kludge to count events; 1 bin histogram with large range
  dummyhist.Sumw2();
  //nominal case
  tree->Project("dummyhist","HT",getCutString(lumiScale_,optfh,selection_,"",0).Data());
  float n = dummyhist.GetBinContent(1);
  float nerr = dummyhist.GetBinError(1);

  TH1D HpdfUnc("HpdfUnc","",100,n - 3*nerr, n+3*nerr);
  for (int i=1; i<=44; i++) { //hard code that there are 44+1 pdf weights
    dummyhist.Reset(); //maybe not needed
    tree->Project("dummyhist","HT",getCutString(lumiScale_,optfh,selection_,"",i).Data());
    HpdfUnc.Fill(dummyhist.GetBinContent(1));
  }

  cout<<" nominal = "<<n<<" +/- "<<nerr<<endl;
  cout<<" w/pdf   = "<<HpdfUnc.GetMean()<<" +/- "<<HpdfUnc.GetRMS()<<endl;

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
  bool oldSaveSetting = savePlots_;
  savePlots_=false;

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

  savePlots_=oldSaveSetting;
}

void drawOwen(TCut extracut="") {

  /*
.L drawReducedTrees.C++
  */
  
  //can't handle a tighter MET cut this way, because it affects the method differently than, say, a tighter HT cut
  assert( !TString(extracut.GetTitle()).Contains("MET"));

  savePlots_=false; //don't save eps,png,pdf files

  loadSamples(false); //make sure to load single top as 3 pieces
  //this only works if loadSamples() hasn't been called yet in this root session!

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
  TString histfilename;
  histfilename.Form("bestM3j-%s.%s.root", vb ? "vb" : "bins", TString(extracut.GetTitle())=="" ? "baseline" : jmt::fortranize(extracut.GetTitle()).Data());
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
    TCut theLSBSelection = baseSelection && passCleaning && passMinDeltaPhi && theBTaggingCut && LSBMET && extracut;
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
    TCut theT2Selection = baseT2Selection && passCleaning && passMinDeltaPhi && theBTaggingCut && T2MET && extracut;
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
    theT2Selection = baseT2Selection && passCleaning && passMinDeltaPhi && theBTaggingCut && SRMET && extracut;
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
    TCut theSRSelection = baseSelection && passCleaning && passMinDeltaPhi && theBTaggingCut && SRMET && extracut;
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
    TCut invdpSelection = baseSelection && passCleaning && failMinDeltaPhi && theBTaggingCut && SRMET && extracut;
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

	theT2Selection = baseT2Selection && passCleaning && passMinDeltaPhi && theBTaggingCut && METselection && extracut;
	TString nameOfT2Hist;
	nameOfT2Hist.Form( "bestM3j_met_%d_%d_t2_%stag_",metCutLow,metCutHigh,btagstring.Data());

	TCut theSelection = baseSelection && passCleaning && passMinDeltaPhi && theBTaggingCut && METselection && extracut;
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
      	theSelection = baseSelection && passCleaning && failMinDeltaPhi && theBTaggingCut && METselection && extracut;
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
