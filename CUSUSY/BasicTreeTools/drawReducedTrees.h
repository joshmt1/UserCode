// -*- C++ -*-

//useful for playing around with plots in interactive ROOT
TH1D* hinteractive=0;

//holds a list of the *active* samples (to be plotted)
std::vector<TString> samples_;
//hold a list of all samples
std::set<TString> samplesAll_;
//these maps use the sample names as keys
std::map<TString, TFile*> files_;
std::map<TString, TH1D*> histos_;
std::map<TString, UInt_t> sampleColor_;
std::map<TString, TString> sampleOwenName_;
std::map<TString, TString> sampleLabel_;
std::map<TString, UInt_t> sampleMarkerStyle_;
TFile* fdata=0;
TH1D* hdata=0;

//default selection
TString selection_ ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && cutMET==1 && cutDeltaPhi==1 && cutCleaning==1";

float leg_x1 = 0.696, leg_x2=0.94, leg_y1=0.5, leg_y2=0.92;
double lumiScale_ = 1.;

bool quiet_=false;
bool doRatio_=false;
bool logy_=false;
bool dostack_=true;
bool doleg_=true;
bool dodata_=true;
bool addOverflow_=true;
//bool doSubtraction_=false;
bool drawQCDErrors_=false;
bool renormalizeBins_=false;//no setter function
bool owenColor_ = false;

bool normalized_=false;

bool useFlavorHistoryWeights_=false;//no setter function
float flavorHistoryScaling_=-1;

bool savePlots_ = true; //no setter function
bool drawTotalSM_=false; //no setter function
bool drawTotalSMSusy_=false;//no setter function
bool drawSusyOnly_=false;//no setter function
bool drawMarkers_=true;//no setter function

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

// == set configuration options ==
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

void setLumiScale(double lumiscale){
  lumiScale_ = lumiscale;
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

void renewLegend() {

  if (leg!=0) delete leg;
  leg = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2);
  leg->SetBorderSize(0);
  leg->SetLineStyle(0);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);

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

void fillFlavorHistoryScaling() {

  TTree* tree=0;
  for (unsigned int isample=0; isample<samples_.size(); isample++) {
    if ( samples_[isample].Contains("WJets")) { //will pick out WJets or WJetsZ2
      tree = (TTree*) files_[samples_[isample]]->Get("reducedTree");
    }
  }
  if (tree==0) {cout<<"Did not find a WJets sample!"<<endl; return;}

  gROOT->cd();
  TH1D dummyU("dummyU","",1,0,1e9);
  TH1D dummyk("dummyk","",1,0,1e9);
  assert(0);//please check that the following cuts have the weights and cuts handled correctly
  tree->Draw("HT>>dummyU","1","goff");
  tree->Draw("HT>>dummyk","flavorHistoryWeight","goff");
  flavorHistoryScaling_ = dummyU.Integral() / dummyk.Integral();
  if (!quiet_) cout<<"flavor history scaling factor = "<<flavorHistoryScaling_<<endl;

}

TString getCutString(double lumiscale= 1., TString extraWeight="", TString thisSelection="", TString extraSelection="", int pdfWeightIndex=0) {
  TString weightedcut="weight"; 
  
  weightedcut += "*(";
  weightedcut +=lumiscale;
  weightedcut+=")";
  
  if (extraWeight=="flavorHistoryWeight") {
    if (flavorHistoryScaling_ <0) {
      fillFlavorHistoryScaling();
    }
    extraWeight.Form("flavorHistoryWeight*%f",flavorHistoryScaling_);
  }
  if (extraWeight!="") {
    weightedcut += "*(";
    weightedcut +=extraWeight;
    weightedcut+=")";
  }
  if (pdfWeightIndex != 0) {
    TString pdfString;
    pdfString.Form("*pdfWeights[%d]",pdfWeightIndex);
    weightedcut += pdfString;
  }
  if (thisSelection!="") {
    weightedcut += "*(";
    weightedcut+=thisSelection;
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

//add a sample to be plotted to the *end* of the list
void addSample(const TString & newsample) {

  //see if it is already there
  for (std::vector<TString>::iterator it = samples_.begin(); it!=samples_.end(); ++it) {
    if ( *it == newsample) {
      cout<<newsample<<" is already on the list!"<<endl;
      return;
    }
  }
  //if it isn't there, go ahead and add it

  if ( samplesAll_.find(newsample) != samplesAll_.end() ) {
    samples_.push_back(newsample);
  }
  else {
    cout<<"Could not find sample with name "<<newsample<<endl;
  }
}

void clearSamples() {
  samples_.clear();
}

void removeSample(const TString & sample) {

  for (std::vector<TString>::iterator it = samples_.begin(); it!=samples_.end(); ++it) {
    if ( *it == sample) {
      if (!quiet_) cout<<sample<<" removed from plotting list!"<<endl;
      samples_.erase(it);
      return;
    }
  }

  //if we get down here then we didn't find the sample
  cout<<sample<<" could not be found on the plotting list, so I did not remove it!"<<endl;
}

void loadSamples(bool joinSingleTop=true) {
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
  if (joinSingleTop) samples_.push_back("SingleTop");
  else {
    samples_.push_back("SingleTop-sChannel");
    samples_.push_back("SingleTop-tChannel");
    samples_.push_back("SingleTop-tWChannel");
  }
  samples_.push_back("WJets");

  samples_.push_back("ZJets");
  samples_.push_back("Zinvisible");
  samples_.push_back("LM13");

  //samplesAll_ should have *every* available sample
  //samplesAll_.insert("QCD");
  //samplesAll_.insert("PythiaQCD");
  samplesAll_.insert("PythiaPUQCD");
  //samplesAll_.insert("PythiaPUQCDFlat");
  samplesAll_.insert("TTbarJets");
  samplesAll_.insert("WJets");
  samplesAll_.insert("ZJets");
  //samplesAll_.insert("WJetsZ2");
  //samplesAll_.insert("ZJetsZ2");
  samplesAll_.insert("Zinvisible");
  samplesAll_.insert("SingleTop");
  samplesAll_.insert("SingleTop-sChannel");
  samplesAll_.insert("SingleTop-tChannel");
  samplesAll_.insert("SingleTop-tWChannel");
  samplesAll_.insert("LM13");
  samplesAll_.insert("LM9");

  //these blocks are just a "dictionary"
  //no need to ever comment these out
  if (!owenColor_) {
    sampleColor_["LM13"] = kRed-9;//kGray;
    sampleColor_["QCD"] = kYellow;
    sampleColor_["PythiaQCD"] = kYellow;
    sampleColor_["PythiaPUQCD"] = kYellow;
    sampleColor_["PythiaPUQCDFlat"] = kYellow;
    sampleColor_["TTbarJets"]=kRed+1;
    sampleColor_["SingleTop"] = kMagenta;
    sampleColor_["WJets"] = kGreen-3;
    sampleColor_["WJetsZ2"] = kGreen-3;
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
    sampleColor_["WJetsZ2"] = kOrange;
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
  sampleLabel_["PythiaQCD"] = "QCD (Z2)";
  sampleLabel_["PythiaPUQCDFlat"] = "QCD (Z2+PU)"; 
  sampleLabel_["PythiaPUQCD"] = "QCD (Z2+PU)";
  sampleLabel_["TTbarJets"]="t#bar{t}";
  sampleLabel_["SingleTop"] = "Single-Top";
  sampleLabel_["WJets"] = "W#rightarrowl#nu";
  sampleLabel_["WJetsZ2"] = "W#rightarrowl#nu (Z2)";
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
  sampleMarkerStyle_["WJetsZ2"] = kMultiply;
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
  sampleOwenName_["WJetsZ2"] = "wjets";
  sampleOwenName_["ZJets"] = "zjets";
  sampleOwenName_["Zinvisible"] = "zinvis";
  sampleOwenName_["SingleTop-sChannel"] = "singletops";
  sampleOwenName_["SingleTop-tChannel"] = "singletopt";
  sampleOwenName_["SingleTop-tWChannel"] = "singletoptw";
  sampleOwenName_["TotalSM"] = "totalsm";
  sampleOwenName_["Total"] = "total";  

  for (std::set<TString>::iterator isample=samplesAll_.begin(); isample!=samplesAll_.end(); ++isample) {
    TString fname="reducedTree.";
    fname+=cutdesc;
    fname+=".";
    fname += *isample;
    fname+=".root";
    fname.Prepend(inputPath);
    files_[*isample] = new TFile(fname);
    if (files_[*isample]->IsZombie() ) {cout<<"file error with "<<*isample<<endl; files_[*isample]=0;}
    else { if (!quiet_)    cout<<"Added sample: "<<*isample<<endl;}
  }

  //load data file too
  TString dname="reducedTree.";
  dname+=cutdesc;
  dname+=".data.root";
  dname.Prepend(inputPath);
  if (dname.Contains("JERbias")) {
    dname.ReplaceAll("JERbias_",""); //JERbias not relevant for data
    dname.ReplaceAll("JERbias6_",""); //JERbias not relevant for data
  }
  if ( dodata_) {
    fdata = new TFile(dname);
    if (fdata->IsZombie()) cout<<"Problem with data file! "<<dname<<endl; 
  }

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
   
  TString optfh= useFlavorHistoryWeights_ && samplename.Contains("WJets") ? "flavorHistoryWeight" : "";
  tree->Project(histname,var,getCutString(lumiScale_,optfh,selection_,"",0).Data());
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
    TString hname = jmt::fortranize(var); hname += "_"; hname += samples_[isample];
    histos_[samples_[isample]] = (varbins==0) ? new TH1D(hname,"",nbins,low,high) : new TH1D(hname,"",nbins,varbins);
    histos_[samples_[isample]]->Sumw2();

    //qcd reweighting not implemented yet

    TTree* tree = (TTree*) files_[samples_[isample]]->Get("reducedTree");
    gROOT->cd();
    TString weightopt= useFlavorHistoryWeights_ && samples_[isample].Contains("WJets") ? "flavorHistoryWeight" : "";
    tree->Project(hname,var,getCutString(lumiScale_,weightopt,selection_,"",0).Data());
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
    TString hname = jmt::fortranize(var); hname += "_"; hname += "data";
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

    if (!quiet_)    cout<<"Data underflow: " <<hdata->GetBinContent(0)<<endl;//BEN
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

    if(!quiet_ && dostack_ && nbins<11){//BEN - 11 is an arbitrary number that isn't too big so we don't print out too much stuff.
      for(int i=1; i<=nbins; i++){
	cout << "data: " << hdata->GetBinContent(i) << " +- " << hdata->GetBinError(i) << ", totalsm: " << totalsm->GetBinContent(i) << " +- " << totalsm->GetBinError(i) << ", ratio: " << ratio->GetBinContent(i) << " +- " << ratio->GetBinError(i) << endl;
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
  if (savePlots_) {
    thecanvas->SaveAs(savename+".eps"); //for me
    //  thecanvas->Print(savename+".C");    //for formal purposes
    thecanvas->SaveAs(savename+".pdf"); //for pdftex
    thecanvas->SaveAs(savename+".png"); //for twiki
  }

}

void drawPlots(const TString var, const int nbins, const float* varbins, const TString xtitle, const TString ytitle, TString filename="") {
  //provide a more natural modification to the argument list....
  drawPlots( var, nbins, 0, 1, xtitle,ytitle, filename,  varbins);
}

void drawSignificance(const TString & var, const int nbins, const float low, const float high, const TString & savename) {

  bool oldSaveSetting = savePlots_;
  savePlots_=false;
  drawPlots(var,nbins,low,high,var,"",savename);

  TH1D* hsusy=0;
  for (std::vector<TString>::iterator it = samples_.begin(); it!=samples_.end(); ++it) {
    if ( (*it).Contains("LM") ) hsusy = histos_[ *it];
  }
  if (hsusy==0) {
    cout<<"Didn't find a signal sample"<<endl;
    return;
  }

  for (int ibin= 1; ibin<=nbins; ibin++) {
    double B = totalsm->Integral(ibin,nbins);
    double S = hsusy->Integral(ibin,nbins);
    if (B>0)    cout<<totalsm->GetBinLowEdge(ibin)<<"\t"<<S/sqrt(B)<<endl;//"\t"<<totalsm->GetBinContent(ibin)<<endl;
    else     cout<<ibin<<" B is zero"<<endl;
  }



  savePlots_=oldSaveSetting;
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

  //terrible hack to decide if bias correction should be calculated
  bool calcBiasCorr = false;
  if(nbins==4 && low>-0.01 && low<0.01 && high>199.99 && high<200.01) calcBiasCorr=true;
  float cb_qcd=0, cb_qcd_err=0, cb_sm=0, cb_sm_err=0, cb_data=0, cb_data_err=0;
  float cp_qcd=0, cp_qcd_err=0, cp_sm=0, cp_sm_err=0, cp_data=0, cp_data_err=0;
  float n_qcd_sb = 0, n_qcd_sb_err = 0, n_qcd_sig = 0, n_qcd_sig_err = 0;
  float n_qcd_a = 0, n_qcd_a_err = 0, n_qcd_d = 0, n_qcd_d_err = 0;

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
    if (useFlavorHistoryWeights_) assert(0); // this needs to be implemented
    tree->Project(hnameP,var,getCutString(lumiScale_,"",selection_,cstring1,0).Data());
    tree->Project(hnameF,var,getCutString(lumiScale_,"",selection_,cstring2,0).Data());
    
    if (addOverflow_)  addOverflowBin( histos_[hnameP] );
    if (addOverflow_)  addOverflowBin( histos_[hnameF] );

    //compute ratio
    histos_[hnameR]->Divide(histos_[hnameP], histos_[hnameF]);

    if (!samples_[isample].Contains("LM")) {
      totalsm_pass.Add(histos_[hnameP]);
      totalsm_fail.Add(histos_[hnameF]);

      //comment out filling of these for now to save time
      TH2D this2d_SB("this2d_SB","",50,100,150,50,0,TMath::Pi());
      //      tree->Project("this2d_SB","minDeltaPhi:MET",getCutString().Data());
      TH2D this2d_50("this2d_50","",50,50,100,50,0,TMath::Pi());
      //      tree->Project("this2d_50","minDeltaPhi:MET",getCutString().Data());

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
      
      if(calcBiasCorr){
	cp_qcd = histos_[hnameR]->GetBinContent(3)/ histos_[hnameR]->GetBinContent(2);
	cp_qcd_err = jmt::errAoverB( histos_[hnameR]->GetBinContent(3), histos_[hnameR]->GetBinError(3), histos_[hnameR]->GetBinContent(2), histos_[hnameR]->GetBinError(2)); 
	cb_qcd = histos_[hnameR]->GetBinContent(4)/ histos_[hnameR]->GetBinContent(3);
	cb_qcd_err = jmt::errAoverB( histos_[hnameR]->GetBinContent(4), histos_[hnameR]->GetBinError(4), histos_[hnameR]->GetBinContent(3), histos_[hnameR]->GetBinError(3)); 
	n_qcd_sb = histos_[hnameP]->GetBinContent(3);
	n_qcd_sb_err = histos_[hnameP]->GetBinError(3);
	n_qcd_sig = histos_[hnameP]->GetBinContent(4);
	n_qcd_sig_err = histos_[hnameP]->GetBinError(4);
     	n_qcd_a = histos_[hnameF]->GetBinContent(3);
	n_qcd_a_err = histos_[hnameF]->GetBinError(3);
	n_qcd_d = histos_[hnameF]->GetBinContent(4);
	n_qcd_d_err = histos_[hnameF]->GetBinError(4);
      }
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
  if(calcBiasCorr){
    cp_sm = totalsm->GetBinContent(3)/totalsm->GetBinContent(2);
    cp_sm_err = jmt::errAoverB(totalsm->GetBinContent(3),totalsm->GetBinError(3),totalsm->GetBinContent(2),totalsm->GetBinError(2));
    cb_sm = totalsm->GetBinContent(4)/totalsm->GetBinContent(3);
    cb_sm_err = jmt::errAoverB(totalsm->GetBinContent(4),totalsm->GetBinError(4),totalsm->GetBinContent(3),totalsm->GetBinError(3));
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
    dtree->Project(hnameP,var,getCutString(1.,"",selection_,cstring1,0).Data());
    dtree->Project(hnameF,var,getCutString(1.,"",selection_,cstring2,0).Data());
    if (addOverflow_)  addOverflowBin( histos_[hnameP] );
    if (addOverflow_)  addOverflowBin( histos_[hnameF] );
    //compute ratio
    hdata->Divide(histos_[hnameP], histos_[hnameF]);

    dtree->Project("data2d_SB","minDeltaPhi:MET",getCutString(1.,"",selection_,"",0).Data());
    dtree->Project("data2d_50","minDeltaPhi:MET",getCutString(1.,"",selection_,"",0).Data());

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

    if(calcBiasCorr){
      cp_data = hdata->GetBinContent(3)/hdata->GetBinContent(2);
      cp_data_err = jmt::errAoverB(hdata->GetBinContent(3),hdata->GetBinError(3),hdata->GetBinContent(2),hdata->GetBinError(2)); 
      cb_data = hdata->GetBinContent(4)/hdata->GetBinContent(3);
      cb_data_err = jmt::errAoverB(hdata->GetBinContent(4),hdata->GetBinError(4),hdata->GetBinContent(3),hdata->GetBinError(3)); 
    }
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
  if(calcBiasCorr){
    cout<<endl;
    cout<<"Pseudo bias correction (using 50<MET<100 and 100<MET<150):" << endl;
    cout<<"QCD MC: "<<cp_qcd<<" +/- "<<cp_qcd_err<<endl;
    cout<<"SM MC: "<<cp_sm<<" +/- "<<cp_sm_err<<endl;
    cout<<"Data: "<<cp_data<<" +/- "<<cp_data_err<<endl;
    cout<<endl;
    cout<<"True bias correction (using 100<MET<150 and 150<MET):" << endl;
    cout<<"QCD MC: "<<cb_qcd<<" \\pm "<<cb_qcd_err<<endl;
    cout<<"SM MC: "<<cb_sm<<" \\pm "<<cb_sm_err<<endl;
    cout<<"Data: "<<cb_data<<" \\pm "<<cb_data_err<<endl;
    cout << endl;
    cout << "QCD Event Counts" << endl;
    cout << "SB: " << n_qcd_sb << " +- " << n_qcd_sb_err << endl;
    cout << "SIG: " << n_qcd_sig << " +- " << n_qcd_sig_err << endl;
    cout << "A: " << n_qcd_a << " +- " << n_qcd_a_err << endl;
    cout << "D: " << n_qcd_d << " +- " << n_qcd_d_err << endl;
    cout << endl;
  }
  cout<<"End of drawR()"<<endl;
}
