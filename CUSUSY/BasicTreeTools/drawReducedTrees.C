/*
after replacing drawBasicPlot.C with drawCutflowPlots.C,
I will now work on replacing drawCutflowPlots.C with this code.

This code will not use the plots created by the Nminus1 code,
but rather it will work from the reducedTrees
*/

#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TH1.h"
#include "TLine.h"

//not clear yet whether I want to use HistHolder.h or not

#include <iostream>
#include <map>

//holds a list of the sample names (using same code as file name)
std::vector<TString> samples_;
std::map<TString, TFile*> files_;
std::map<TString, TH1D*> histos_;
std::map<TString, UInt_t> sampleColor_;
std::map<TString, TString> sampleLabel_;
std::map<TString, UInt_t> sampleMarkerStyle_;

TString inputPath = "/cu2/joshmt/";
TString cutdesc = "Baseline0_PF_pfMEThigh_PFLep0e0mu_minDP_MuonCleaning";
TString selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && cutCleaning==1";

bool logy_=false;
bool dostack_=false;
bool doleg_=true;
bool dodata_=true;
bool addOverflow_=true;
//bool doSubtraction_=false;

bool doVerticalLine_=false;
double verticalLinePosition_=0;

bool doCustomPlotMax_=false;
double customPlotMax_=0;

bool doCustomPlotMin_=false;
double customPlotMin_=0;


TCanvas* cnorm=0;
TLegend* leg=0;
THStack* thestack=0;

bool loaded_=false; //bookkeeping

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
}

void setStackMode(bool dostack) {
  dostack_=dostack;
}

void doData(bool dodata) {
  dodata_=dodata;
}

void drawLegend(bool doleg) {
  doleg_=doleg;
}


void addOverflowBin(TH1D* theHist) {
  //this code was written for when there was a customizable plot range (post-histo creation)
  //it could be made a lot simpler now

  int lastVisibleBin = theHist->GetNbinsX();
  //  cout<<theHist<<"  "<<lastVisibleBin<<"\t";

  //in case there is no custom range, the code should just add the overflow bin to the last bin of the histo
  double lastBinContent = theHist->GetBinContent(lastVisibleBin);
  double lastBinError = pow(theHist->GetBinError(lastVisibleBin),2); //square in prep for addition

  cout<<"Overflow addition: "<<lastBinContent<<" +/- "<<sqrt(lastBinError)<<" --> ";

  //now loop over the bins that aren't being shown at the moment (including the overflow bin)
  for (int ibin = lastVisibleBin+1; ibin <= 1 + theHist->GetNbinsX() ; ++ibin) {
    lastBinContent += theHist->GetBinContent( ibin);
    lastBinError += pow(theHist->GetBinError( ibin),2);
  }
  lastBinError = sqrt(lastBinError);

  theHist->SetBinContent(lastVisibleBin,lastBinContent);
  theHist->SetBinError(lastVisibleBin,lastBinError);
  cout<<lastBinContent<<" +/- "<<lastBinError<<endl;
}

void drawVerticalLine() {
  if (cnorm==0) return;

  //this is a fine example of ROOT idiocy
  TVirtualPad* thePad = cnorm->GetPad(0);
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

  //this block controls what samples will enter your plot
  //  samples_.push_back("QCD");
  //  samples_.push_back("PythiaQCD");
  samples_.push_back("PythiaPUQCDFlat");
  samples_.push_back("TTbarJets");
  samples_.push_back("SingleTop");
  samples_.push_back("WJets");
  samples_.push_back("ZJets");
  samples_.push_back("Zinvisible");
  samples_.push_back("LM13");

  //these 3 blocks are just a "dictionary"
  //no need to ever comment these out
  sampleColor_["LM13"] = kGray; //borrowed from a different sample
  sampleColor_["QCD"] = kYellow;
  sampleColor_["PythiaQCD"] = kYellow-7;
  sampleColor_["PythiaPUQCDFlat"] = kYellow-5;
  sampleColor_["TTbarJets"]=kRed+1;
  sampleColor_["SingleTop"] = kMagenta;
  sampleColor_["WJets"] = kGreen-3;
  sampleColor_["ZJets"] = kAzure-2;
  sampleColor_["Zinvisible"] = kOrange-3;

  sampleLabel_["LM13"] = "LM13";
  sampleLabel_["QCD"] = "QCD";
  sampleLabel_["PythiaQCD"] = "QCD (Pythia Z2)";
  sampleLabel_["PythiaPUQCDFlat"] = "QCD (Pileup)";
  sampleLabel_["TTbarJets"]="t#bar{t}";
  sampleLabel_["SingleTop"] = "Single-Top";
  sampleLabel_["WJets"] = "W#rightarrowl#nu";
  sampleLabel_["ZJets"] = "Z/#gamma*#rightarrowl^{+}l^{-}";
  sampleLabel_["Zinvisible"] = "Z#rightarrow#nu#nu";

  sampleMarkerStyle_["LM13"] = kFullStar;
  sampleMarkerStyle_["QCD"] = kFullCircle;
  sampleMarkerStyle_["PythiaQCD"] = kOpenCircle;
  sampleMarkerStyle_["PythiaPUQCDFlat"] = kOpenCircle;
  sampleMarkerStyle_["TTbarJets"]= kFullSquare;
  sampleMarkerStyle_["SingleTop"] = kOpenSquare;
  sampleMarkerStyle_["WJets"] = kMultiply;
  sampleMarkerStyle_["ZJets"] = kFullTriangleUp;
  sampleMarkerStyle_["Zinvisible"] = kFullTriangleDown;

  for (unsigned int isample=0; isample<samples_.size(); isample++) {
    TString fname="reducedTree.";
    fname+=cutdesc;
    fname+=".";
    fname+=samples_[isample];
    fname+=".root";
    fname.Prepend(inputPath);
    files_[samples_[isample]] = new TFile(fname);
    if (files_[samples_[isample]]->IsZombie() ) cout<<"file error with "<<samples_[isample]<<endl;
    else     cout<<"Added sample: "<<samples_[isample]<<endl;
  }

  //load data file too
  TString dname="reducedTree.";
  dname+=cutdesc;
  dname+=".data.root";
  dname.Prepend(inputPath);
  if (dostack_ && dodata_) {
    fdata = new TFile(dname);
    if (fdata->IsZombie()) cout<<"Problem with data file! "<<dname<<endl;
  }

}

void drawStack(const TString var, const int nbins, const float low, const float high, const TString xtitle, const TString ytitle, TString filename="") {
  loadSamples();

  if (filename=="") filename=var;

  //  TH1D* thestackH=0;

  gROOT->SetStyle("CMS");

  if (cnorm!=0) delete cnorm;
  cnorm= new TCanvas("cnorm","the canvas",600,550);
  cnorm->cd()->SetRightMargin(0.04);
  if (logy_) cnorm->SetLogy();

  //FIXME -- could add recode these numbers as a config option
  if (leg!=0) delete leg;
  leg = new TLegend(0.696, 0.35, 0.94, 0.92);
  leg->SetBorderSize(0);
  leg->SetLineStyle(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);

  if (dostack_) {
    if (thestack!= 0 ) delete thestack;
    thestack = new THStack("thestack","--");
  }

  //here is the part that is really different from the previous implementation
  //need to make new histograms
  TString opt="hist e";
  for (unsigned int isample=0; isample<samples_.size(); isample++) {
    cout <<samples_[isample]<<endl;

    if (histos_.find(samples_[isample]) != histos_.end() ) {
      if (histos_[samples_[isample]] != 0) {
	delete  histos_[samples_[isample]];
	histos_[samples_[isample]] = 0;
      }
    }

    gROOT->cd();
    //should each histo have a different name? maybe
    TString hname = var; hname += "_"; hname += samples_[isample];
    histos_[samples_[isample]] = new TH1D(hname,"",nbins,low,high);

    histos_[samples_[isample]]->SetXTitle(xtitle);
    histos_[samples_[isample]]->SetYTitle(ytitle);

    //qcd reweighting not implemented yet

    TTree* tree = (TTree*) files_[samples_[isample]]->Get("reducedTree");
    TString weightedcut="weight"; 
    if (selection_!="") {
      weightedcut += "*(";
      weightedcut+=selection_;
      weightedcut+=")";
    }
    //    cout<<selection_<<"\t"<<weightedcut<<endl;
    gROOT->cd();
    tree->Project(hname,var,weightedcut.Data());
    //now the histo is filled
    
    if (addOverflow_)  addOverflowBin( histos_[samples_[isample]] ); //manipulates the TH1D

    //now just do a bunch of histogram formatting
    if (!dostack_) {
      //set line color instead of fill color for this type of plot
      histos_[samples_[isample]]->SetLineColor(sampleColor_[samples_[isample]]);
      histos_[samples_[isample]]->SetMarkerStyle(sampleMarkerStyle_[samples_[isample]]);
      histos_[samples_[isample]]->SetMarkerColor(sampleColor_[samples_[isample]]);
      
      //ad hoc additions
      histos_[samples_[isample]]->SetLineWidth(2);
      histos_[samples_[isample]]->SetMarkerColor(sampleColor_[samples_[isample]]);
    }
    else {
      histos_[samples_[isample]]->SetFillColor(sampleColor_[samples_[isample]]);
      histos_[samples_[isample]]->SetMarkerSize(0);
    }
    leg->AddEntry(histos_[samples_[isample]], sampleLabel_[samples_[isample]]);

    if (dostack_) {
      thestack->Add(histos_[samples_[isample]] );
//       if (doSubtraction_) {
// 	if (thestackH==0) {
// 	  const TH1D* temphist=hh.find(hname,files_[samples_[isample]]);
// 	  thestackH=new TH1D("thestackH","total MC",temphist->GetNbinsX(),temphist->GetXaxis()->GetXmin(),temphist->GetXaxis()->GetXmax());
// 	  thestackH->Sumw2();
// 	}
// 	thestackH->Add(hh.find(hname,files_[samples_[isample]]));
//       }
    }
    else {
      histos_[samples_[isample]]->Draw(opt);
      if (!opt.Contains("same")) opt+=" same";
    }
  }

  if (!dostack_) {
    //not implemented yet
    //    hh.normalize(); //normalize all histos to 1
    //hh.SetMaximum( hh.GetMaximum()*1.05);
    //should implement use of customPlotMax_, i guess
  }
  else {
    thestack->Draw("hist");
    thestack->GetHistogram()->GetXaxis()->SetTitle(xtitle);
    thestack->GetHistogram()->GetYaxis()->SetTitle(ytitle);

    if (doVerticalLine_) drawVerticalLine(); //i want to draw the data last

    if (dodata_) {
      gROOT->cd();
      cout<<"Drawing data!"<<endl;
      if (hdata != 0) delete hdata;
      TString hname = var; hname += "_"; hname += "data";
      hdata = new TH1D(hname,"",nbins,low,high);
      
      TTree* dtree = (TTree*) fdata->Get("reducedTree");
      gROOT->cd();
      dtree->Project(hname,var,selection_.Data());
      //now the histo is filled

      hdata->UseCurrentStyle(); //maybe not needed anymore
      hdata->SetMarkerColor(kBlack);
      hdata->SetLineWidth(2);
      hdata->SetMarkerStyle(kFullCircle);
      hdata->SetMarkerSize(1);
      if (addOverflow_)     addOverflowBin(hdata); // manipulates the histogram!

//       if (doSubtraction_) { //plot data - MC
// 	hdataSubtracted = (TH1D*)hdata->Clone("hdataSubtracted");
// 	hdataSubtracted->Add(thestackH,-1);
// 	cout<<hdataSubtracted->Integral()<<" "<<hdata->Integral()<<endl;

// 	hdataSubtracted->UseCurrentStyle();
// 	hdataSubtracted->SetMarkerColor(kBlue);
// 	hdataSubtracted->SetLineColor(kBlue);
// 	hdataSubtracted->SetLineWidth(2);
// 	hdataSubtracted->SetMarkerStyle(kOpenCircle);
// 	hdataSubtracted->SetMarkerSize(1);
//       }

      hdata->Draw("SAME");
      //      cout<<hdata->GetMaximum()<<"\t"<<thestack->GetMaximum()<<endl;
      if (hdata->GetMaximum() > thestack->GetMaximum()) {
	thestack->SetMaximum( hdata->GetMaximum());
      }
      //      if (doSubtraction_) hdataSubtracted->Draw("SAME");
    }
    if (doCustomPlotMax_) thestack->SetMaximum(customPlotMax_);
    if (doCustomPlotMin_) thestack->SetMinimum(customPlotMin_);
  }

  if (doleg_)  leg->Draw();

  TString savename = filename;
  //  if (doSubtraction_) savename+="-MCSub";
  if (logy_) savename += "-logY";
  //  savename += scaleAppendToFilename;

  if (!dostack_)    savename += "-drawNormalized";
  else savename += "-drawStack";

  //amazingly, \includegraphics cannot handle an extra dot in the filename. so avoid it.
  cnorm->SaveAs(savename+".eps"); //for me
  //  cnorm->Print(savename+".C");    //for formal purposes
  cnorm->SaveAs(savename+".pdf"); //for pdftex
  cnorm->SaveAs(savename+".png"); //for twiki


}


void drawSomething() {

  /* for reference, here is all cuts (>=0 b)
selection_ ="nbjets>=0 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && cutCleaning==1";
  */

  //here we imitate exactly the plots we made in plotWithScale() of the old code

  setStackMode(true);
  doData(true);

  //the strength and weakness of this tree-based method is that I need to apply the cuts now!
  //most straightforward way is to manually set the cut string before each plot

  int nbins,low,high;
  TString var,xtitle;

  // ==== MET plots ====
  setLogY(true);   setPlotMinimum(1e-3);
  nbins=50; low= 0; high=260;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";

  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && cutCleaning==1";
  drawStack(var,nbins,low,high,xtitle,"Events", "H_MET");

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && cutCleaning==1";
  drawStack(var,nbins,low,high,xtitle,"Events", "H_MET_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && cutCleaning==1";
  drawStack(var,nbins,low,high,xtitle,"Events", "H_MET_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && cutCleaning==1";
  drawStack(var,nbins,low,high,xtitle,"Events", "H_MET_ge2b");

  // ==== minDeltaPhi plots ====
  setLogY(false);
  //shit...forgot to put this variable in the tree!
  nbins=25; low=0; high = TMath::Pi() + 0.001;
  var="minDeltaPhi"; xtitle="min(#Delta#phi[ jets 1..3, E_{T}^{miss} ] )";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutCleaning==1";
  drawStack(var,nbins,low,high,xtitle,"Events", "HminDeltaPhiMETj_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutCleaning==1";
  drawStack(var,nbins,low,high,xtitle,"Events", "HminDeltaPhiMETj_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutCleaning==1";
  drawStack(var,nbins,low,high,xtitle,"Events", "HminDeltaPhiMETj_ge2b");

  // ==== HT plots ====
  setLogY(false);
  //shit...forgot to put this variable in the tree!
  nbins=25; low=0; high = TMath::Pi() + 0.001;
  var="minDeltaPhi"; xtitle="min(#Delta#phi[ jets 1..3, E_{T}^{miss} ] )";
  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutCleaning==1";
  drawStack(var,nbins,low,high,xtitle,"Events", "HminDeltaPhiMETj_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutCleaning==1";
  drawStack(var,nbins,low,high,xtitle,"Events", "HminDeltaPhiMETj_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutCleaning==1";
  drawStack(var,nbins,low,high,xtitle,"Events", "HminDeltaPhiMETj_ge2b");

}
