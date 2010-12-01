//arguably stupid to redo the work that went into drawBasicPlot.C
//but I'm doing it anyway (lots of copy/paste)
//also:
//as much as I hate it, it might be nice to make this macro compilable

//QCD and SingleTop need to be added with hadd

#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "/afs/cern.ch/user/j/joshmt/root/util/HistHolder.h"
//i am using histHolder, but reimplementing a lot of stuff in PlotUtil,etc

#include <iostream>
#include <map>

//holds a list of the sample names (using same code as file name)
std::vector<TString> samples_;
std::map<TString, TFile*> files_;
std::map<TString, UInt_t> sampleColor_;
std::map<TString, TString> sampleLabel_;
std::map<TString, UInt_t> sampleMarkerStyle_;

bool loaded_=false; //bookkeeping

bool customrange_active=false;
double customrange_low=0;
double customrange_high=0;

void setPlotRange(double low,double high) {
  customrange_low=low;
  customrange_high=high;
  customrange_active=true;
}

void resetPlotRange() {
  customrange_active=false;
}

void loadSamples() {
  loaded_=true;
  const TString cutdesc = "Baseline0_PF_pfMEThigh_PFLep_minDP";

  samples_.push_back("QCD");
  samples_.push_back("TTbarJets");
  samples_.push_back("SingleTop");
  samples_.push_back("WJets");
  samples_.push_back("ZJets");
  samples_.push_back("Zinvisible");
  samples_.push_back("LM13");

  sampleColor_["LM13"] = kGray; //borrowed from a different sample
  sampleColor_["QCD"] = kYellow;
  sampleColor_["TTbarJets"]=kRed+1;
  sampleColor_["SingleTop"] = kMagenta;
  sampleColor_["WJets"] = kGreen-3;
  sampleColor_["ZJets"] = kAzure-2;
  sampleColor_["Zinvisible"] = kOrange-3;

  sampleLabel_["LM13"] = "LM13";
  sampleLabel_["QCD"] = "QCD";
  sampleLabel_["TTbarJets"]="t#bar{t}";
  sampleLabel_["SingleTop"] = "Single-Top";
  sampleLabel_["WJets"] = "W#rightarrowl#nu";
  sampleLabel_["ZJets"] = "Z#rightarrowl^{+}l^{-}";
  sampleLabel_["Zinvisible"] = "Z#rightarrow#nu#nu";

  sampleMarkerStyle_["LM13"] = kFullStar;
  sampleMarkerStyle_["QCD"] = kFullCircle;
  sampleMarkerStyle_["TTbarJets"]= kFullSquare;
  sampleMarkerStyle_["SingleTop"] = kOpenSquare;
  sampleMarkerStyle_["WJets"] = kMultiply;
  sampleMarkerStyle_["ZJets"] = kFullTriangleUp;
  sampleMarkerStyle_["Zinvisible"] = kFullTriangleDown;

  for (unsigned int isample=0; isample<samples_.size(); isample++) {
    TString fname="cutflowPlots.";
    fname+=cutdesc;
    fname+=".";
    fname+=samples_[isample];
    fname+=".root";
    files_[samples_[isample]] = new TFile(fname);
    if (files_[samples_[isample]]->IsZombie() ) cout<<"file error with "<<samples_[isample]<<endl;
    else     cout<<"Added sample: "<<samples_[isample]<<endl;
  }

}

TCanvas* cnorm;
TLegend* leg;
void drawNormalized(const TString hname, const TString xtitle, const TString ytitle)
{

  gROOT->SetStyle("CMS");

  if (!loaded_)  loadSamples();

  cnorm= new TCanvas("cnorm","normalized",600,400);
  cnorm->cd()->SetRightMargin(0.04);

  leg = new TLegend(0.696, 0.35, 0.94, 0.92);
  leg->SetBorderSize(0);
  //leg->SetFillColor(0);
  //leg->SetLineColor(1);
  leg->SetLineStyle(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);

  //first load the desired histogram from each file
  HistHolder hh;
  TString opt="hist e";
  for (unsigned int isample=0; isample<samples_.size(); isample++) {
    hh.load(hname,files_[samples_[isample]]);
    hh.find(hname,files_[samples_[isample]])->SetXTitle(xtitle);
    hh.find(hname,files_[samples_[isample]])->SetYTitle(ytitle);

    if (customrange_active) {
      hh.find(hname,files_[samples_[isample]])->GetXaxis()->SetRangeUser(customrange_low,customrange_high);
    }

    //set line color instead of fill color for this type of plot
    hh.find(hname,files_[samples_[isample]])->SetLineColor(sampleColor_[samples_[isample]]);
    hh.find(hname,files_[samples_[isample]])->SetMarkerStyle(sampleMarkerStyle_[samples_[isample]]);
    hh.find(hname,files_[samples_[isample]])->SetMarkerColor(sampleColor_[samples_[isample]]);

    //ad hoc additions
    hh.find(hname,files_[samples_[isample]])->SetLineWidth(2);
    hh.find(hname,files_[samples_[isample]])->SetMarkerColor(sampleColor_[samples_[isample]]);

    leg->AddEntry( hh.find(hname,files_[samples_[isample]]), sampleLabel_[samples_[isample]]);
    hh.find(hname,files_[samples_[isample]])->Draw(opt);
    if (!opt.Contains("same")) opt+=" same";
  }
  hh.normalize(); //normalize all histos to 1
  hh.SetMaximum( hh.GetMaximum()*1.05);

  leg->Draw();

  //amazingly, \includegraphics cannot handle an extra dot in the filename. so avoid it.
  cnorm->SaveAs(hname+"-drawNormalized.eps"); //for me
  //  cnorm->Print(hname+".drawNormalized.C");    //for formal purposes
  cnorm->SaveAs(hname+"-drawNormalized.pdf"); //for pdftex
  cnorm->SaveAs(hname+"-drawNormalized.png"); //for twiki

}


void plotStuff()
{
  //this first one is not so informative because of HT=0 events
  drawNormalized("H_HT_cutInclusive","HT [GeV]","Arbitrary units");
  drawNormalized("H_NJets_cutHT","# of jets","Arbitrary units");

  setPlotRange(0,5);
  drawNormalized("H_NElectrons_cut3Jets","# of electrons","Arbitrary units");
  drawNormalized("H_NMuons_cut3Jets","# of muons","Arbitrary units");

  resetPlotRange();
  drawNormalized("H_MET_cutMuVeto","E^{miss}_{T} [GeV]","Arbitrary units");
  setPlotRange(0,5);
  drawNormalized("H_NBJets_cutMET","# of b-tagged jets","Arbitrary units");

}
