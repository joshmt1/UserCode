/*
--- ROOT version ---
source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.18/x86_64-slc5-gcc46-opt/root/bin/thisroot.csh

--- setup --

You must have symlink to or copy of MiscUtil.cxx in the working directory.
This is available at:
https://github.com/joshmt1/UserCode/blob/master/MiscUtil.cxx

For interactive mode, you must also make a symlink as follows:
ln -s drawReducedTrees.h drawDelphes.h


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
#include "TH3.h"
#include "TLine.h"
#include "TCut.h"
#include "TLatex.h"
#include "TChain.h"

//from Luke
#include "TSelectorMultiDraw.h"

//container for mSugra cross-sections
#include "CrossSectionTable.h"

// can be checked out of UserCode/joshmt
#include "MiscUtil.cxx"

#include <fstream>

#include <iostream>
#include <iomanip>
#include <map>
#include <set>

TString inputPath = "";
TString dataInputPath =  "";

double lumiScale_ = 3000e3;//300 fb-1

#include "drawReducedTrees.h"

void initSamples(TString option="all skimmed") {

  option.ReplaceAll("all","tt bj nm1 nm2 nm3 stoc");
  option.ReplaceAll("signal","nm1 nm2 nm3 stoc");

  //v10 is the first version with New samples (including new b-tagging)
  //v09 is the latest and greatest (old samples)
  //pre-approval was with v07

  TString pathbase ="/cu3/joshmt/Upgrade/";
  if (option.Contains("PhaseI_50")) pathbase += "PhaseI_50PileUp/";
  else pathbase += "PhaseII_Configuration4v2_140PileUp/";

  if (option.Contains("skimmed"))
    inputPath = pathbase+"v10/skimmed/"; //delcmscornell
  else if (option.Contains("htskim"))
    inputPath = pathbase+"v10/ht1000/"; //delcmscornell
  else if (option=="comp7to10")
    inputPath = pathbase+"comp7to10/";
  else 
    inputPath = pathbase+"v10/"; //delcmscornell

  dodata_=false; //no data for upgrade studies!
  cmEnergy_=14; //TeV
  reducedTreeName_ = "simpleTree";


  addToSamplesAll("tt-4p-0-600-v1510_14TEV");
  addToSamplesAll("tt-4p-600-1100-v1510_14TEV");
  addToSamplesAll("tt-4p-1100-1700-v1510_14TEV");
  addToSamplesAll("tt-4p-1700-2500-v1510_14TEV");
  addToSamplesAll("tt-4p-2500-100000-v1510_14TEV");

  addToSamplesAll("tj-4p-0-500-v1510_14TEV");
  addToSamplesAll("tj-4p-500-1000-v1510_14TEV");
  addToSamplesAll("tj-4p-1000-1600-v1510_14TEV");
  addToSamplesAll("tj-4p-1600-2400-v1510_14TEV");
  addToSamplesAll("tj-4p-2400-100000-v1510_14TEV");

  addToSamplesAll("tB-4p-0-500-v1510_14TEV");
  addToSamplesAll("tB-4p-500-900-v1510_14TEV");
  addToSamplesAll("tB-4p-900-1500-v1510_14TEV");
  addToSamplesAll("tB-4p-1500-2200-v1510_14TEV");
  addToSamplesAll("tB-4p-2200-100000-v1510_14TEV");
 
  addToSamplesAll("LL-4p-0-100-v1510_14TEV");
  addToSamplesAll("LL-4p-100-200-v1510_14TEV");
  addToSamplesAll("LL-4p-200-500-v1510_14TEV");
  addToSamplesAll("LL-4p-500-900-v1510_14TEV");
  addToSamplesAll("LL-4p-900-1400-v1510_14TEV");
  addToSamplesAll("LL-4p-1400-100000-v1510_14TEV");

  addToSamplesAll("ttB-4p-0-900-v1510_14TEV");
  addToSamplesAll("ttB-4p-1600-2500-v1510_14TEV");
  addToSamplesAll("ttB-4p-2500-100000-v1510_14TEV");
  addToSamplesAll("ttB-4p-900-1600-v1510_14TEV");
  addToSamplesAll("LLB-4p-0-400-v1510_14TEV");
  addToSamplesAll("LLB-4p-400-900-v1510_14TEV");
  addToSamplesAll("LLB-4p-900-100000-v1510_14TEV");
  addToSamplesAll("Bjj-vbf-4p-0-700-v1510_14TEV");
  addToSamplesAll("Bjj-vbf-4p-1400-2300-v1510_14TEV");
  //  addToSamplesAll("Bjj-vbf-4p-2300-3400_14TEV");
  addToSamplesAll("Bjj-vbf-4p-2300-3400-v1510_14TEV");
  addToSamplesAll("Bjj-vbf-4p-700-1400-v1510_14TEV");
  addToSamplesAll("BB-4p-0-300-v1510_14TEV");
  addToSamplesAll("BB-4p-1300-2100-v1510_14TEV");
  addToSamplesAll("BB-4p-2100-100000-v1510_14TEV");
  addToSamplesAll("BB-4p-300-700-v1510_14TEV");
  addToSamplesAll("BB-4p-700-1300-v1510_14TEV");
  addToSamplesAll("BBB-4p-0-600-v1510_14TEV");
  addToSamplesAll("BBB-4p-1300-100000-v1510_14TEV");
  addToSamplesAll("BBB-4p-600-1300-v1510_14TEV");


  addToSamplesAll("Bj-4p-0-300-v1510_14TEV");
  addToSamplesAll("Bj-4p-1100-1800-v1510_14TEV");
  addToSamplesAll("Bj-4p-1800-2700-v1510_14TEV");
  addToSamplesAll("Bj-4p-2700-3700-v1510_14TEV");
  addToSamplesAll("Bj-4p-300-600-v1510_14TEV");
  addToSamplesAll("Bj-4p-3700-100000-v1510_14TEV");
  addToSamplesAll("Bj-4p-600-1100-v1510_14TEV");

  addToSamplesAll("B-4p-0-1-v1510_14TEV");

  TString signal1 = "naturalModel1";//old name: "susyhit_Scenario1_v02";
  TString signal2 = "naturalModel2";
  TString signal3 = "naturalModel3";
  TString signal4 = "stoc";
  addToSamplesAll(signal1);
  addToSamplesAll(signal2);
  addToSamplesAll(signal3);
  addToSamplesAll(signal4);


  if (option=="comp7to10") {

    addToSamplesAll("tt-4p-0-600-v1510_14TEV_v07");
    addToSamplesAll("tt-4p-600-1100-v1510_14TEV_v07");
    addToSamplesAll("tt-4p-1100-1700-v1510_14TEV_v07");
    addToSamplesAll("tt-4p-1700-2500-v1510_14TEV_v07");
    addToSamplesAll("tt-4p-2500-100000-v1510_14TEV_v07");

    addToSamplesAll("tt-4p-0-600-v1510_14TEV_v10");
    addToSamplesAll("tt-4p-600-1100-v1510_14TEV_v10");
    addToSamplesAll("tt-4p-1100-1700-v1510_14TEV_v10");
    addToSamplesAll("tt-4p-1700-2500-v1510_14TEV_v10");
    addToSamplesAll("tt-4p-2500-100000-v1510_14TEV_v10");

  }

  loadSamples("delphes");
  usePUweight_=false;
  useTrigEff_=false;

  isPreliminary_=false;
  currentConfig_=configDescriptions_.getDefault();
 clearSamples();
  resetChains(); //important in the case that we call initHiggsSamples() more than once

  if (option.Contains("combinesm")) {
    TString smname = "tt-4p-0-600-v1510_14TEV";
    addSample(smname,kAzure-3,"Total SM");
    chainSamples(smname,"tt-4p-600-1100-v1510_14TEV");
    chainSamples(smname,"tt-4p-1100-1700-v1510_14TEV");
    chainSamples(smname,"tt-4p-1700-2500-v1510_14TEV");
    chainSamples(smname,"tt-4p-2500-100000-v1510_14TEV");
    
    TString tj_name = "tj-4p-0-500-v1510_14TEV"; 
    chainSamples(smname,tj_name);
    chainSamples(smname,"tj-4p-500-1000-v1510_14TEV");
    chainSamples(smname,"tj-4p-1000-1600-v1510_14TEV");
    chainSamples(smname,"tj-4p-1600-2400-v1510_14TEV");
    chainSamples(smname,"tj-4p-2400-100000-v1510_14TEV");

    TString tb_name = "tB-4p-0-500-v1510_14TEV"; 
    chainSamples(smname,tb_name);
    chainSamples(smname,"tB-4p-500-900-v1510_14TEV");
    chainSamples(smname,"tB-4p-900-1500-v1510_14TEV");
    chainSamples(smname,"tB-4p-1500-2200-v1510_14TEV");
    chainSamples(smname,"tB-4p-2200-100000-v1510_14TEV");
    
    TString ll_name = "LL-4p-0-100-v1510_14TEV";
    chainSamples(smname,ll_name);
    chainSamples(smname,"LL-4p-100-200-v1510_14TEV");
    chainSamples(smname,"LL-4p-200-500-v1510_14TEV");
    chainSamples(smname,"LL-4p-500-900-v1510_14TEV");
    chainSamples(smname,"LL-4p-900-1400-v1510_14TEV");
    chainSamples(smname,"LL-4p-1400-100000-v1510_14TEV");

    TString bj_name = "Bj-4p-0-300-v1510_14TEV";
    chainSamples(smname,bj_name);
    chainSamples(smname,"Bj-4p-1100-1800-v1510_14TEV");
    chainSamples(smname,"Bj-4p-1800-2700-v1510_14TEV");
    chainSamples(smname,"Bj-4p-2700-3700-v1510_14TEV");
    chainSamples(smname,"Bj-4p-300-600-v1510_14TEV");
    chainSamples(smname,"Bj-4p-3700-100000-v1510_14TEV");
    chainSamples(smname,"Bj-4p-600-1100-v1510_14TEV");

    chainSamples(smname,"B-4p-0-1-v1510_14TEV");

    chainSamples(smname,"BB-4p-0-300-v1510_14TEV");
    chainSamples(smname,"BB-4p-1300-2100-v1510_14TEV");
    chainSamples(smname,"BB-4p-2100-100000-v1510_14TEV");
    chainSamples(smname,"BB-4p-300-700-v1510_14TEV");
    chainSamples(smname,"BB-4p-700-1300-v1510_14TEV");
    chainSamples(smname,"BBB-4p-0-600-v1510_14TEV");
    chainSamples(smname,"BBB-4p-1300-100000-v1510_14TEV");
    chainSamples(smname,"BBB-4p-600-1300-v1510_14TEV");
    chainSamples(smname,"Bjj-vbf-4p-0-700-v1510_14TEV");
    chainSamples(smname,"Bjj-vbf-4p-1400-2300-v1510_14TEV");
    //    chainSamples(smname,"Bjj-vbf-4p-2300-3400_14TEV");
    chainSamples(smname,"Bjj-vbf-4p-2300-3400-v1510_14TEV");
    chainSamples(smname,"Bjj-vbf-4p-700-1400-v1510_14TEV");
    chainSamples(smname,"LLB-4p-0-400-v1510_14TEV");
    chainSamples(smname,"LLB-4p-400-900-v1510_14TEV");
    chainSamples(smname,"LLB-4p-900-100000-v1510_14TEV");
    chainSamples(smname,"ttB-4p-0-900-v1510_14TEV");
    chainSamples(smname,"ttB-4p-1600-2500-v1510_14TEV");
    chainSamples(smname,"ttB-4p-2500-100000-v1510_14TEV");
    chainSamples(smname,"ttB-4p-900-1600-v1510_14TEV");

  }
  else if (option.Contains("combinefs")) {
   TString smname = "tt-4p-0-600-v1510_14TEV";
    addSample(smname,kAzure-3,"Total SM");
    chainSamples(smname,"tt-4p-600-1100-v1510_14TEV");
    chainSamples(smname,"tt-4p-1100-1700-v1510_14TEV");
    chainSamples(smname,"tt-4p-1700-2500-v1510_14TEV");
    chainSamples(smname,"tt-4p-2500-100000-v1510_14TEV");
    
    TString tj_name = "tj-4p-0-500-v1510_14TEV"; 
    chainSamples(smname,tj_name);
    chainSamples(smname,"tj-4p-500-1000-v1510_14TEV");
    chainSamples(smname,"tj-4p-1000-1600-v1510_14TEV");
    chainSamples(smname,"tj-4p-1600-2400-v1510_14TEV");
    chainSamples(smname,"tj-4p-2400-100000-v1510_14TEV");

    TString tb_name = "tB-4p-0-500-v1510_14TEV"; 
    chainSamples(smname,tb_name);
    chainSamples(smname,"tB-4p-500-900-v1510_14TEV");
    chainSamples(smname,"tB-4p-900-1500-v1510_14TEV");
    chainSamples(smname,"tB-4p-1500-2200-v1510_14TEV");
    chainSamples(smname,"tB-4p-2200-100000-v1510_14TEV");
    /*
    chainSamples(smname,"ttB-4p-0-900-v1510_14TEV");
    chainSamples(smname,"ttB-4p-1600-2500-v1510_14TEV");
    chainSamples(smname,"ttB-4p-2500-100000-v1510_14TEV");
    chainSamples(smname,"ttB-4p-900-1600-v1510_14TEV");
    */
    if (option.Contains("combinefssig")) { //chain signal
      cout<<"Adding signal to SM chain!"<<endl;
      chainSamples(smname,signal1);
    }
  }
  else if (option.Contains("combinenonfs")) {
 
    TString smname = "LL-4p-0-100-v1510_14TEV";
    addSample(smname,kGreen,"Total non-FS");
    chainSamples(smname,"LL-4p-100-200-v1510_14TEV");
    chainSamples(smname,"LL-4p-200-500-v1510_14TEV");
    chainSamples(smname,"LL-4p-500-900-v1510_14TEV");
    chainSamples(smname,"LL-4p-900-1400-v1510_14TEV");
    chainSamples(smname,"LL-4p-1400-100000-v1510_14TEV");

    TString bj_name = "Bj-4p-0-300-v1510_14TEV";
    chainSamples(smname,bj_name);
    chainSamples(smname,"Bj-4p-1100-1800-v1510_14TEV");
    chainSamples(smname,"Bj-4p-1800-2700-v1510_14TEV");
    chainSamples(smname,"Bj-4p-2700-3700-v1510_14TEV");
    chainSamples(smname,"Bj-4p-300-600-v1510_14TEV");
    chainSamples(smname,"Bj-4p-3700-100000-v1510_14TEV");
    chainSamples(smname,"Bj-4p-600-1100-v1510_14TEV");

    chainSamples(smname,"B-4p-0-1-v1510_14TEV");

    chainSamples(smname,"BB-4p-0-300-v1510_14TEV");
    chainSamples(smname,"BB-4p-1300-2100-v1510_14TEV");
    chainSamples(smname,"BB-4p-2100-100000-v1510_14TEV");
    chainSamples(smname,"BB-4p-300-700-v1510_14TEV");
    chainSamples(smname,"BB-4p-700-1300-v1510_14TEV");
    chainSamples(smname,"BBB-4p-0-600-v1510_14TEV");
    chainSamples(smname,"BBB-4p-1300-100000-v1510_14TEV");
    chainSamples(smname,"BBB-4p-600-1300-v1510_14TEV");
    chainSamples(smname,"Bjj-vbf-4p-0-700-v1510_14TEV");
    chainSamples(smname,"Bjj-vbf-4p-1400-2300-v1510_14TEV");
    //    chainSamples(smname,"Bjj-vbf-4p-2300-3400_14TEV");
    chainSamples(smname,"Bjj-vbf-4p-2300-3400-v1510_14TEV");
    chainSamples(smname,"Bjj-vbf-4p-700-1400-v1510_14TEV");
    chainSamples(smname,"LLB-4p-0-400-v1510_14TEV");
    chainSamples(smname,"LLB-4p-400-900-v1510_14TEV");
    chainSamples(smname,"LLB-4p-900-100000-v1510_14TEV");

    //    chainSamples(smname,"ttB-4p-0-900-v1510_14TEV");
    //    chainSamples(smname,"ttB-4p-1600-2500-v1510_14TEV");
    //    chainSamples(smname,"ttB-4p-2500-100000-v1510_14TEV");
    //    chainSamples(smname,"ttB-4p-900-1600-v1510_14TEV");

  }
  else if (option=="comp7to10") {
    TString tt_name = "tt-4p-0-600-v1510_14TEV_v07";
    addSample(tt_name,kRed,"v07 tt");
    chainSamples(tt_name,"tt-4p-600-1100-v1510_14TEV_v07");
    chainSamples(tt_name,"tt-4p-1100-1700-v1510_14TEV_v07");
    chainSamples(tt_name,"tt-4p-1700-2500-v1510_14TEV_v07");
    chainSamples(tt_name,"tt-4p-2500-100000-v1510_14TEV_v07");
    
    tt_name = "tt-4p-0-600-v1510_14TEV_v10";
    addSample(tt_name,kBlue,"v10 tt");
    chainSamples(tt_name,"tt-4p-600-1100-v1510_14TEV_v10");
    chainSamples(tt_name,"tt-4p-1100-1700-v1510_14TEV_v10");
    chainSamples(tt_name,"tt-4p-1700-2500-v1510_14TEV_v10");
    chainSamples(tt_name,"tt-4p-2500-100000-v1510_14TEV_v10");
   
  }
  else {
 
    if (option.Contains("bj")) {

      TString ot_name="BBB-4p-0-600-v1510_14TEV";
      addSample(ot_name,kOrange+2,"Other SM");
      //      chainSamples(bb_name,"BBB-4p-0-600-v1510_14TEV");
      chainSamples(ot_name,"BBB-4p-1300-100000-v1510_14TEV");
      chainSamples(ot_name,"BBB-4p-600-1300-v1510_14TEV");

      chainSamples(ot_name,"LL-4p-0-100-v1510_14TEV");
      chainSamples(ot_name,"LL-4p-100-200-v1510_14TEV");
      chainSamples(ot_name,"LL-4p-200-500-v1510_14TEV");
      chainSamples(ot_name,"LL-4p-500-900-v1510_14TEV");
      chainSamples(ot_name,"LL-4p-900-1400-v1510_14TEV");
      chainSamples(ot_name,"LL-4p-1400-100000-v1510_14TEV");
      
      chainSamples(ot_name,"LLB-4p-0-400-v1510_14TEV");
      chainSamples(ot_name,"LLB-4p-400-900-v1510_14TEV");
      chainSamples(ot_name,"LLB-4p-900-100000-v1510_14TEV");

      chainSamples(ot_name,"ttB-4p-0-900-v1510_14TEV");
      chainSamples(ot_name,"ttB-4p-1600-2500-v1510_14TEV");
      chainSamples(ot_name,"ttB-4p-2500-100000-v1510_14TEV");
      chainSamples(ot_name,"ttB-4p-900-1600-v1510_14TEV");
    }

    if (option.Contains("tt")) {

      TString tb_name = "tB-4p-0-500-v1510_14TEV"; 
      addSample(tb_name,kCyan-6,"Single top");
      //      chainSamples(tt_name,tb_name); //combine all top samples into one
      chainSamples(tb_name,"tB-4p-500-900-v1510_14TEV");
      chainSamples(tb_name,"tB-4p-900-1500-v1510_14TEV");
      chainSamples(tb_name,"tB-4p-1500-2200-v1510_14TEV");
      chainSamples(tb_name,"tB-4p-2200-100000-v1510_14TEV");

      chainSamples(tb_name,"tj-4p-0-500-v1510_14TEV");
      chainSamples(tb_name,"tj-4p-500-1000-v1510_14TEV");
      chainSamples(tb_name,"tj-4p-1000-1600-v1510_14TEV");
      chainSamples(tb_name,"tj-4p-1600-2400-v1510_14TEV");
      chainSamples(tb_name,"tj-4p-2400-100000-v1510_14TEV");
  
    }
        if (option.Contains("bj")) {

      TString bb_name = "BB-4p-0-300-v1510_14TEV";
      addSample(bb_name,kPink+3,"VV");
      //      chainSamples(bj_name,bb_name); //combine these with V+jets
      chainSamples(bb_name,"BB-4p-1300-2100-v1510_14TEV");
      chainSamples(bb_name,"BB-4p-2100-100000-v1510_14TEV");
      chainSamples(bb_name,"BB-4p-300-700-v1510_14TEV");
      chainSamples(bb_name,"BB-4p-700-1300-v1510_14TEV");
    }


    if (option.Contains("bj")) {

      TString bj_name = "Bj-4p-0-300-v1510_14TEV";
      addSample(bj_name,kViolet+5,"V+jets");
      chainSamples(bj_name,"Bj-4p-1100-1800-v1510_14TEV");
      chainSamples(bj_name,"Bj-4p-1800-2700-v1510_14TEV");
      chainSamples(bj_name,"Bj-4p-2700-3700-v1510_14TEV");
      chainSamples(bj_name,"Bj-4p-300-600-v1510_14TEV");
      chainSamples(bj_name,"Bj-4p-3700-100000-v1510_14TEV");
      chainSamples(bj_name,"Bj-4p-600-1100-v1510_14TEV");
      chainSamples(bj_name,"Bjj-vbf-4p-0-700-v1510_14TEV");
      chainSamples(bj_name,"Bjj-vbf-4p-1400-2300-v1510_14TEV");
      //      chainSamples(bj_name,"Bjj-vbf-4p-2300-3400_14TEV");
      chainSamples(bj_name,"Bjj-vbf-4p-2300-3400-v1510_14TEV");
      chainSamples(bj_name,"Bjj-vbf-4p-700-1400-v1510_14TEV");
      chainSamples(bj_name,"B-4p-0-1-v1510_14TEV"); //      addSample("B-4p-0-1-v1510_14TEV",kOrange+3,"B");

    }
   if (option.Contains("tt") ) {
      TString tt_name = "tt-4p-0-600-v1510_14TEV";
      addSample(tt_name,kAzure-4,"t#bar{t}");
      chainSamples(tt_name,"tt-4p-600-1100-v1510_14TEV");
      chainSamples(tt_name,"tt-4p-1100-1700-v1510_14TEV");
      chainSamples(tt_name,"tt-4p-1700-2500-v1510_14TEV");
      chainSamples(tt_name,"tt-4p-2500-100000-v1510_14TEV");
    
    }
  }
  //add signal last so that it sits on top of the stack
  //add signal last so that it sits on top of the stack
  if (option.Contains("nm3"))   addSample(signal3,6,"NM3");  
  if (option.Contains("nm2"))   addSample(signal2,4,"NM2");  
  if (option.Contains("nm1")) addSample(signal1,2,"NM1");

  if (option.Contains("stoc"))  addSample(signal4,kPink+8,"STOC");  


}

