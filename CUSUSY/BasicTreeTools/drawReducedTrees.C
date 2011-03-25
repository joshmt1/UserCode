/*
after replacing drawBasicPlot.C with drawCutflowPlots.C,
I will now work on replacing drawCutflowPlots.C with this code.

This code will not use the plots created by the Nminus1 code,
but rather it will work from the reducedTrees

known bugs:
non-stacked drawing is not yet implemented
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
std::map<TString, TString> sampleLabel_;
std::map<TString, UInt_t> sampleMarkerStyle_;

TString inputPath = "/cu2/joshmt/";
//TString cutdesc = "Baseline0_PF_JERbias_pfMEThigh_PFLep0e0mu_minDP_MuonEcalCleaning";
TString cutdesc = "Baseline0_PF_JERbias_pfMEThigh_PFLepRA20e0mu_minDP_MuonEcalCleaning";
//TString cutdesc = "Baseline0_PF_pfMEThigh_PFLepRA20e0mu_minDP_MuonEcalCleaning";
//TString cutdesc = "Baseline0_PF_pfMEThigh_PFLep0e0mu_minDP_MuonEcalCleaning";
TString selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && cutCleaning==1";

float leg_x1 = 0.696, leg_x2=0.94, leg_y1=0.5, leg_y2=0.92;

bool doRatio_=false;
bool logy_=false;
bool dostack_=true;
bool doleg_=true;
bool dodata_=true;
bool addOverflow_=true;
//bool doSubtraction_=false;
bool drawQCDErrors_=false;

bool doVerticalLine_=false;
double verticalLinePosition_=0;

bool doCustomPlotMax_=false;
double customPlotMax_=0;

bool doCustomPlotMin_=false;
double customPlotMin_=0;


TCanvas* thecanvas=0;
//TCanvas* cratio=0;
TLegend* leg=0;
THStack* thestack=0;
TH1D* totalsm=0;
TH1D* totalewk=0;
TH1D* totalqcdttbar=0;
TH1D* ratio=0; float ratioMin=0; float ratioMax=2;
TGraphErrors* qcderrors=0;
bool loaded_=false; //bookkeeping

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

int mainpadWidth = 600; int mainpadHeight=550;
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
    cout<< thecanvas->GetPad(1)->GetXlowNDC() <<"\t"
	<< thecanvas->GetPad(1)->GetWNDC() <<"\t"
	<< thecanvas->GetPad(1)->GetYlowNDC() <<"\t"
	<< thecanvas->GetPad(1)->GetHNDC() <<endl;
    if (logy_) thecanvas->GetPad(1)->SetLogy();
  }
  else { if (logy_) thecanvas->SetLogy(); }


  int cdarg = opt.Contains("ratio") ? 1 : 0;
  thecanvas->cd(cdarg);

}


void resetHistos() {
  for ( std::map<TString, TH1D*>::iterator i = histos_.begin(); i!=histos_.end(); ++i) {
    if (i->second != 0) {
      delete  i->second;
      i->second= 0;
    }
  }
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
  cout<<weightedcut<<endl;
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
  cout<<lastBinContent<<" +/- "<<lastBinError<<endl;
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

  //this block controls what samples will enter your plot
  //careful -- QCD must have 'QCD' in its name somewhere.
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
  sampleColor_["PythiaQCD"] = kYellow;
  sampleColor_["PythiaPUQCDFlat"] = kYellow;
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
  if (dname.Contains("JERbias")) dname.ReplaceAll("JERbias_",""); //JERbias not relevant for data
  if ( dodata_) {
    fdata = new TFile(dname);
    if (fdata->IsZombie()) cout<<"Problem with data file! "<<dname<<endl;
  }

}

float drawSimple(const TString var, const int nbins, const float low, const float high, const TString filename, 
		 const TString histname , const TString samplename) {

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
	cout <<samples_[isample]<<endl;
	tree = (TTree*) files_[samples_[isample]]->Get("reducedTree");
      }
    }
  }
  if (tree==0) {cout<<"Something went wrong finding your sample!"<<endl; return 0;}
  gROOT->cd();
  
  //when owen explicitly uses a histo type, it is a TH1F
  TH1F hh(histname,histname,nbins,low,high);
  hh.Sumw2();
   
  tree->Project(histname,var,getCutString().Data());
  float theIntegral = hh.Integral(0,nbins+1);

  //at this point i've got a histogram. what more could i want?
  TFile fout(filename,"UPDATE");
  hh.Write();
  fout.Close();
  return theIntegral;
}

void drawPlots(const TString var, const int nbins, const float low, const float high, const TString xtitle, const TString ytitle, TString filename="") {
  loadSamples();

  if (filename=="") filename=var;

  //  TH1D* thestackH=0;

  gROOT->SetStyle("CMS");
  //gStyle->SetHatchesLineWidth(1);

  TString canvasOpt = doRatio_ ? "ratio" : "";
  const int mainPadIndex = doRatio_ ? 1 : 0;
  renewCanvas(canvasOpt);

  thecanvas->cd(mainPadIndex);

  //FIXME -- could add recode these numbers as a config option
  if (leg!=0) delete leg;
  leg = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2);
  leg->SetBorderSize(0);
  leg->SetLineStyle(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);

  if (dostack_) {
    if (thestack!= 0 ) delete thestack;
    thestack = new THStack("thestack","--");
    if (totalsm!=0) delete totalsm;
    totalsm = new TH1D("totalsm","",nbins,low,high);
    totalsm->Sumw2();
    if (totalewk!=0) delete totalewk;
    totalewk = new TH1D("totalewk","",nbins,low,high);
    totalewk->Sumw2();
    if (totalqcdttbar!=0) delete totalqcdttbar;
    totalqcdttbar = new TH1D("totalqcdttbar","",nbins,low,high);
    totalqcdttbar->Sumw2();
    if (doRatio_) {
      if (ratio!=0) delete ratio;
      ratio = new TH1D("ratio","data/(SM MC)",nbins,low,high);
      ratio->Sumw2();
    }
  }
  //here is the part that is really different from the previous implementation
  //need to make new histograms
  resetHistos(); //delete existing histograms
  TString opt="hist e";
  for (unsigned int isample=0; isample<samples_.size(); isample++) {
    cout <<samples_[isample]<<endl;

    gROOT->cd();
    //should each histo have a different name? maybe
    TString hname = var; hname += "_"; hname += samples_[isample];
    histos_[samples_[isample]] = new TH1D(hname,"",nbins,low,high);
    histos_[samples_[isample]]->Sumw2();

    histos_[samples_[isample]]->SetXTitle(xtitle);
    histos_[samples_[isample]]->SetYTitle(ytitle);

    //qcd reweighting not implemented yet

    TTree* tree = (TTree*) files_[samples_[isample]]->Get("reducedTree");
    gROOT->cd();
    tree->Project(hname,var,getCutString().Data());
    //now the histo is filled
    
    if (addOverflow_)  addOverflowBin( histos_[samples_[isample]] ); //manipulates the TH1D

    //if we're going to draw QCD errors, create a TGraphErrors from the QCD histogram
    if (drawQCDErrors_ && samples_[isample].Contains("QCD")) {
      if (qcderrors!=0) delete qcderrors;
      qcderrors = new TGraphErrors(histos_[samples_[isample]]);
      qcderrors->SetFillStyle(3353);
      qcderrors->SetFillColor(1);
    }
    if (!samples_[isample].Contains("LM")) {
      totalsm->Add(histos_[samples_[isample]]);
        cout << "totalsm: " << samples_[isample] << endl;
    }
    if (!samples_[isample].Contains("LM") && !samples_[isample].Contains("QCD") && !samples_[isample].Contains("TTbar")) {
      totalewk->Add(histos_[samples_[isample]]);
      cout << "totalewk: " << samples_[isample] << endl;
    }
    if (samples_[isample].Contains("QCD") || samples_[isample].Contains("TTbar")){
      totalqcdttbar->Add(histos_[samples_[isample]]);
      cout << "totalqcdttbar: " << samples_[isample] << endl;
    }

    //now just do a bunch of histogram formatting
    if (!dostack_) {
      //set line color instead of fill color for this type of plot
      histos_[samples_[isample]]->SetLineColor(sampleColor_[samples_[isample]]);
      histos_[samples_[isample]]->SetMarkerStyle(sampleMarkerStyle_[samples_[isample]]);
      histos_[samples_[isample]]->SetMarkerColor(sampleColor_[samples_[isample]]);
      
      //ad hoc additions
      histos_[samples_[isample]]->SetLineWidth(2);
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

    if (drawQCDErrors_) qcderrors->Draw("2 same");

    if (dodata_) {
      gROOT->cd();
      cout<<"Drawing data!"<<endl;
      if (hdata != 0) delete hdata;
      TString hname = var; hname += "_"; hname += "data";
      hdata = new TH1D(hname,"",nbins,low,high);
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

      cout<<"Integral of data, total SM: "<<hdata->Integral()<<" ; "<<totalsm->Integral()<<endl;
      cout<<"Chi^2 Test results: "<<hdata->Chi2Test(totalsm,"UW P")<<endl;
      cout<<"KS Test results: "<<hdata->KolmogorovTest(totalsm,"N")<<endl;;
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
    if (doCustomPlotMax_) thestack->SetMaximum(customPlotMax_);
    if (doCustomPlotMin_) thestack->SetMinimum(customPlotMin_);
  }

  thecanvas->cd(mainPadIndex);
  if (doleg_)  leg->Draw();

  //  if (doSubtraction_) savename+="-MCSub";
  TString savename = filename;
  if (logy_) savename += "-logY";
  //  savename += scaleAppendToFilename;

  if (!dostack_)    savename += "-drawPlain";
  else savename += "-drawStack";

  //amazingly, \includegraphics cannot handle an extra dot in the filename. so avoid it.
  thecanvas->SaveAs(savename+".eps"); //for me
  //  thecanvas->Print(savename+".C");    //for formal purposes
  thecanvas->SaveAs(savename+".pdf"); //for pdftex
  thecanvas->SaveAs(savename+".png"); //for twiki

}

//could add xtitle and ytitle
void drawR(const TString vary, const float cutVal, const int nbins, const float low, const float high) {

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
  TString drawopt="";
  float max=-1e9; TString firsthist="";
  for (unsigned int isample=0; isample<samples_.size(); isample++) {

    cout <<samples_[isample]<<endl;
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

    //   cout<<"content of bin 2: "<<histos_[hnameP]->GetBinContent(2)<<" / "<< histos_[hnameF]->GetBinContent(2)<<" = "<<histos_[hnameR]->GetBinContent(2)<<endl;

    //now format the histograms
    cout<<"setting color to: "<<sampleColor_[samples_[isample]]<<endl;
    histos_[hnameR]->SetLineColor(sampleColor_[samples_[isample]]);
    histos_[hnameR]->SetMarkerStyle(sampleMarkerStyle_[samples_[isample]]);
    histos_[hnameR]->SetMarkerColor(sampleColor_[samples_[isample]]);

    //ad hoc additions
    histos_[hnameR]->SetLineWidth(2);

    //draw
    thecanvas->cd(1);
    histos_[hnameR]->Draw(drawopt);
    if (!drawopt.Contains("same")) drawopt+=" same";

    if (firsthist="") firsthist = hnameR;
    if (histos_[hnameR]->GetMaximum() > max) max = histos_[hnameR]->GetMaximum();
  }
  histos_[firsthist]->SetMaximum( max*1.05);
  hinteractive =  histos_[firsthist];

  if (dodata_) {
    gROOT->cd();
    cout<<"Drawing data!"<<endl;
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

    //    hdata->UseCurrentStyle(); //maybe not needed anymore
    hdata->SetMarkerColor(kBlack);
    hdata->SetLineWidth(2);
    hdata->SetMarkerStyle(kFullCircle);
    hdata->SetMarkerSize(1);

    thecanvas->cd(1);
    hdata->Draw("SAME");
    if (hdata->GetMaximum() > max)  histos_[firsthist]->SetMaximum( hdata->GetMaximum());

    //    cratio->cd();
    thecanvas->cd(2);
    if (ratio!=0) delete ratio;
    ratio = new TH1D("ratio","data/(SM MC)",nbins,low,high);
    ratio->Sumw2();
    ratio->Divide(hdata,histos_[firsthist]);
    ratio->SetMinimum(ratioMin);
    ratio->SetMaximum(ratioMax);
    ratio->Draw();
  }


}

void drawrplots() {
  /*
.L drawReducedTrees.C++
  */

  setStackMode(false);
  doData(false);

  //no met, no mindeltaphi
selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutCleaning==1";
setLogY(true);
 drawR("minDeltaPhi",0.3,50,0,250);

 //njets == 3
selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && cutCleaning==1 && njets==3";
setLogY(true);
 drawR("minDeltaPhi",0.3,50,0,250);

 //4 jets
selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && cutCleaning==1 && njets==4";
setLogY(true);
 drawR("minDeltaPhi",0.3,50,0,250);

 //5 jets
selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cutEleVeto==1 && cutMuVeto==1 && cutCleaning==1 && njets==5";
setLogY(true);
 drawR("minDeltaPhi",0.3,50,0,250);

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

  //  setQCDErrorMode(true);
  // PROBLEMS -- CMS style ruins the hash marks. central value on the TGraphErrors is wrong!
  //might be easiest to accumulate a total MC uncertainty and plot that!

  //the strength and weakness of this tree-based method is that I need to apply the cuts now!
  //most straightforward way is to manually set the cut string before each plot

  int nbins;
  float low,high;
  TString var,xtitle;

  // ==== MET plots ====
  setLogY(true);   setPlotMinimum(1e-2);
  nbins=25; low= 0; high=260;
  var="MET"; xtitle="E_{T}^{miss} [GeV]";
  ratioMin = 0; ratioMax = 2;

  selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "H_MET");

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

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HminDeltaPhiMETj_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HminDeltaPhiMETj_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && passInconsistentMuon==1 && passBadPFMuon==1";
  drawPlots(var,nbins,low,high,xtitle,"Events", "HminDeltaPhiMETj_ge2b");

  // ==== HT plots ==== //in this case this is not an N-1 plot, because we've already cut on HT
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
  nbins=10; low=50; high= 550;
  resetPlotMinimum();
  setLogY(false);
  var="bestTopMass"; xtitle="Best top mass (GeV)";

  selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>=100 && MET<150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "bestTopMass_ge1b");

  selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>=100 && MET<150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "bestTopMass_eq1b");

  selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1 && MET>=100 && MET<150";
  drawPlots(var,nbins,low,high,xtitle,"Events", "bestTopMass_ge2b");

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


  // not yet in ntuple
//   // ==== lead b jet pT ====
//   nbins=20; low=30; high= 200;
//   resetPlotMinimum();
//   setLogY(false);
//   var="bjetpt1"; xtitle="p_{T} of lead jet (GeV)";

//   selection_ ="nbjets>=1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
//   drawPlots(var,nbins,low,high,xtitle,"Events", "Hjetpt1_ge1b");

//   selection_ ="nbjets==1 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
//   drawPlots(var,nbins,low,high,xtitle,"Events", "Hjetpt1_eq1b");

//   selection_ ="nbjets>=2 && cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && passInconsistentMuon==1 && passBadPFMuon==1";
//   drawPlots(var,nbins,low,high,xtitle,"Events", "Hjetpt1_ge2b");

}


void drawOwen(bool doMinDPPass) {

  /*
.L drawReducedTrees.C++
  */
  
  doOverflowAddition(false);

  const  int nbins = 35;
  const  float min=0;
  const  float max=800;

  const  TCut baseSelection = "cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1"; //no MET, no minDeltaPhi, no cleaning, no b tag
  // "cutCleaning==1"; //apply all cleaning 
  const  TCut passCleaning ="passInconsistentMuon == 1 && passBadPFMuon==1"; //apply cleaning but without ECAL dead cells
  const  TCut passMinDeltaPhi = "cutDeltaPhi==1";
  const  TCut failMinDeltaPhi = "cutDeltaPhi==0";
  const  TCut ge1b =  "nbjets >= 1";
  const  TCut ge2b =  "nbjets >= 2";
  const  TCut eq1b =  "nbjets == 1";

  //
  TString histfilename= doMinDPPass ? "bestM3j-bins.pass.root" : "bestM3j-bins.fail.root";
  TString textfilename= doMinDPPass ? "drawOwen.pass.output" : "drawOwen.fail.output";
  TCut theMinDeltaPhiCut = doMinDPPass ? passMinDeltaPhi : failMinDeltaPhi;

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
    TCut LSBMET = "MET>=0 && MET<50";
    TCut theLSBSelection = baseSelection && passCleaning && theMinDeltaPhiCut && theBTaggingCut && LSBMET;
    selection_ = theLSBSelection.GetTitle();
    drawSimple("bestTopMass",nbins,min,max,histfilename, "bestM3j_met_0_50_"+btagstring+"tag_data","data");
    drawSimple("bestTopMass",nbins,min,max,histfilename, "bestM3j_met_0_50_"+btagstring+"tag_qcd","PythiaPUQCDFlat");
    drawPlots("bestTopMass",nbins,min,max,"","","deleteme");
    TFile fh(histfilename,"UPDATE");  
    totalsm->SetName("bestM3j_met_0_50_"+btagstring+"tag_totalsm");
    totalewk->SetName("bestM3j_met_0_50_"+btagstring+"tag_totalewk");
    totalqcdttbar->SetName("bestM3j_met_0_50_"+btagstring+"tag_totalqcdttbar");
    totalsm->Write();          
    totalewk->Write();
    totalqcdttbar->Write();
    fh.Close();   
    
    TCut SIGMET = "MET>=150 && MET<10000";
    TCut theSIGSelection = baseSelection && passCleaning && theMinDeltaPhiCut && theBTaggingCut && SIGMET;
    selection_ = theSIGSelection.GetTitle();
    drawSimple("bestTopMass",nbins,min,max,histfilename, "bestM3j_met_150_10000_"+btagstring+"tag_qcd","PythiaPUQCDFlat");
    drawSimple("bestTopMass",nbins,min,max,histfilename, "bestM3j_met_150_10000_"+btagstring+"tag_ttbar","TTbarJets");
    drawSimple("bestTopMass",nbins,min,max,histfilename, "bestM3j_met_150_10000_"+btagstring+"tag_data","data");
    

    //this one is going to be handled by the flexible MET version....
    //   TCut SBMET = "MET>=70 && MET<150";
    //   TCut theSBSelection = baseSelection && passCleaning && theMinDeltaPhiCut && theBTaggingCut && SBMET;
    //   selection_ = theSBSelection.GetTitle();
    //   drawSimple("bestTopMass",nbins,min,max,histfilename, "bestM3j_met_70_150_ge1btag_ttbar","TTbarJets");
    
    //now for a flexible MET region
    //  int metCutLow = 60;
    //  int metCutHigh = 140;
    for (int metCutLow = 70; metCutLow <=100; metCutLow+=30) {
      for (int metCutHigh = 150; metCutHigh <=150; metCutHigh+=5) {
	TString metCutString; metCutString.Form("MET >= %d && MET < %d",metCutLow,metCutHigh);
	ofile<<metCutString<<endl;
	TCut METselection(metCutString.Data());
	TCut theSelection = baseSelection && passCleaning && theMinDeltaPhiCut && theBTaggingCut && METselection;
	selection_ = theSelection.GetTitle();
	TString nameOfHist;
	nameOfHist.Form( "bestM3j_met_%d_%d_%stag_",metCutLow,metCutHigh,btagstring.Data());
	ofile<<"ttbar = "<<drawSimple("bestTopMass",nbins,min,max,histfilename, nameOfHist+"ttbar","TTbarJets")<<endl;
	ofile<<"qcd   = "<<drawSimple("bestTopMass",nbins,min,max,histfilename, nameOfHist+"qcd","PythiaPUQCDFlat")<<endl;
	drawSimple("bestTopMass",nbins,min,max,histfilename, nameOfHist+"data","data");
	drawPlots("bestTopMass",nbins,min,max,"","","deleteme");
	TFile fh(histfilename,"UPDATE");
	totalsm->SetName(nameOfHist+"totalsm");
	totalewk->SetName(nameOfHist+"totalewk");
	totalqcdttbar->SetName(nameOfHist+"totalqcdttbar");
	totalsm->Write();
	totalewk->Write();
	totalqcdttbar->Write();
	fh.Close();
      }
    }
    
  }
  ofile.close();

  //trick to get the total SM histo filled. make sure SingleTop is commented out in the master list
  /*
  drawPlots("bestTopMass",nbins,min,max,"","","deleteme");
  TFile fh(histfilename,"UPDATE");
  totalsm->SetName("bestM3j_met_80_150_ge1btag_allsm");
  totalsm->Write();
  fh.Close();

  drawSimple("bestTopMass",nbins,min,max,histfilename, "bestM3j_met_80_150_ge1btag_data","data"); //for the real fit
  drawSimple("bestTopMass",nbins,min,max,histfilename, "bestM3j_met_80_150_ge1btag_qcd","PythiaPUQCDFlat");
  drawSimple("bestTopMass",nbins,min,max,histfilename, "bestM3j_met_80_150_ge1btag_ttbar","TTbarJets");

  drawSimple("bestTopMass",nbins,min,max,histfilename, "bestM3j_met_80_150_ge1btag_wjets","WJets");
  drawSimple("bestTopMass",nbins,min,max,histfilename, "bestM3j_met_80_150_ge1btag_zjets","ZJets");

  drawSimple("bestTopMass",nbins,min,max,histfilename, "bestM3j_met_80_150_ge1btag_zinvis","Zinvisible");
  */
}
