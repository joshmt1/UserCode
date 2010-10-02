#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TString.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLine.h"
#include "TText.h"
#include "TF1.h"

#include "TStyle.h"

#include <iostream>

TCanvas * Cabcd;
TCanvas * Cfit;

/*
usage:

.L doABCD.C++
ABCD a("MET","minDeltaPhiMET")
a.doABCD()

*/

//const TString cutDescription_ = "RA2METminDP_NoMET_NoDeltaPhi.ge2b";
const TString cutDescription_ = "RA2METminDP_NoMET.ge0b";

class ABCD {
public :
  ABCD(TString xvar,TString yvar);
  ~ABCD();
  void doABCD();

  TH2D* Habcd; //for user picture
  //for fitting and calculation (extended ABCD)
  TH1D* Hx_L;
  TH1D* Hx_H;
  TH1D* Hx_ratio;
  TF1* expfunc_;
private :
  void loopOverTree(TTree* theTree, const bool countSignalRegion, bool doExtendedEstimate=false);

  bool doFit_;

  TString xvar_;
  double xl_low;
  double xl_high;
  double xh_low;
  double xh_high;

  int fitnum_;

  double xplotmin;
  double xplotmax;

  TString yvar_;
  double yl_low;
  double yl_high;
  double yh_low;
  double yh_high;

  double yplotmin;
  double yplotmax;

  bool flippedY_;

  double extendedEstimate_;
  double extendedEstimateErr_;
  double total[2][2];
  double error[2][2];

};

ABCD::ABCD(TString xvar,TString yvar) :
  Habcd(0), Hx_L(0), Hx_H(0), Hx_ratio(0), expfunc_(0),
  doFit_(true),
  xvar_(xvar),
  xl_low(0),xl_high(0),xh_low(0),xh_high(0),
  fitnum_(10), //number of bins for extended ABCD fit
  xplotmin(0),xplotmax(0),
  yvar_(yvar),
  yl_low(0),yl_high(0),yh_low(0),yh_high(0),
  yplotmin(0),yplotmax(0),
  flippedY_(false),
  extendedEstimate_(0),
  extendedEstimateErr_(0)
{
  gStyle->SetOptFit(1);

  for (int i=0; i<2; i++) {
    for (int j=0; j<2; j++) {
      total[i][j] = error[i][j] = 0;
      //      cout<<total[i][j]<<"\t"<<  error[i][j] <<endl;
    }
  }

  //define the boundaries of the different regions
  if (xvar_=="MET" || xvar_=="MHT") {
    xl_low = 50;  //40
    xl_high = 150;
    
    xh_low = 150; 
    xh_high = 1e6;

    xplotmin=xl_low; //*0.9
    xplotmax=440;
  }
  else {std::cout<<"ERROR -- xvar is not valid!"<<std::endl;}

  if (yvar_ =="minDeltaRbj") {
    yl_low = 0.5;
    yl_high = 2;
    
    yh_low = 2; 
    yh_high = 50;

    yplotmin=yl_low; // *0.9
    yplotmax=5.5;//6;
    flippedY_=true;
    doFit_=false;
  }
  else if (yvar_=="DeltaPhiMPTMET") {
    yl_low = 0;
    yl_high = 1.8;
    yh_low = yl_high; 
    yh_high = TMath::Pi();
    yplotmin=yl_low;
    yplotmax=yh_high;
    flippedY_=true;
    doFit_=false;
  }
  else if (yvar_=="minDeltaPhiMET" || yvar_=="minDeltaPhiMHT") {
    yl_low = 0;
    yl_high = 0.3;
    
    yh_low = 0.3; 
    yh_high = TMath::Pi();
    yplotmin=yl_low;
    yplotmax=yh_high;
    flippedY_=false;
  }
  else {std::cout<<"ERROR -- yvar is not valid!"<<std::endl;}

}

ABCD::~ABCD()
{
}


void ABCD::loopOverTree(TTree* theTree, const bool countSignalRegion, bool doExtendedEstimate) {

  /*
this code got a little ugly with the extended estimate
break into its own function?
  */

  double xval;
  double yval;
  double weight;
  theTree->SetBranchAddress(xvar_.Data(),&xval);
  theTree->SetBranchAddress(yvar_.Data(),&yval);
  theTree->SetBranchAddress("weight",&weight);


  int y0=1;
  int y1=0;
  if (flippedY_) {
    y0=0;
    y1=1;
  }

  for (Long64_t ievent = 0; ievent<theTree->GetEntries() ; ievent++) {
    theTree->LoadTree(ievent);
    theTree->GetEntry(ievent);

    if (!doExtendedEstimate)  Habcd->Fill(xval,yval,weight);

    int xbinid=-1;
    if ((xval >= xl_low) && (xval < xl_high) )      xbinid=0;
    else if ((xval >= xh_low) && (xval < xh_high) ) xbinid=1;

    if (xbinid== -1) continue; //not in one of the bins

    int ybinid=-1;
    if ((yval >= yl_low) && (yval < yl_high) )      ybinid=y0;
    else if ((yval >= yh_low) && (yval < yh_high) ) ybinid=y1;    

    if (ybinid== -1) continue; //not in one of the bins

    //skip the signal region counting for backgrounds
    if (xbinid==1 && ybinid==0 && !countSignalRegion && !doExtendedEstimate) continue;

    if (doExtendedEstimate) {
      //we are only interested in the high MET control region
      if (xbinid==1 && ybinid==1 ) {
	double weightedest=weight * expfunc_->Eval(xval);
	extendedEstimate_ += weightedest;
	extendedEstimateErr_ += weightedest*weightedest;
      }
    }
    else {
      total[xbinid][ybinid] += weight;
      error[xbinid][ybinid] += weight*weight;
      
      //for a fit
      if (doFit_) {
	// in principle we could use an array of TH1D* based on ybinid, instead of these if statements
	//do i use y0 and y1 here, or 0 and 1?
	if (xbinid ==0 && ybinid==1) Hx_L->Fill(xval,weight);
	else if (xbinid ==0 && ybinid==0) Hx_H->Fill(xval,weight);
      }
    }
  }

}



void ABCD::doABCD() {

  Cabcd=new TCanvas("Cabcd");

  if (doFit_) Cfit = new TCanvas("Cfit");
  //make a loop over events and pull events from tree

  //count events in each box for naive ABCD approach

  //can do a fit to the low MET ratio, then a second loop over events
  //to get an extended estimate (Ben/RA2 approach)

  //could even contemplate setting up a 2D fit for uncorrelated variables

  //option to handle QCD tree plus 'background' trees
  //this option seems to work for the simple ABCD. Not tested for extended method (probably works?)

  const TString dir="/cu1/joshmt/ABCDtrees/"; //must include trailing /

  //-----------------------------------------
  const  bool dottbar=false;
  const  bool dosignal=false;

  TString myname=dir+"ABCDtree.";
  myname += cutDescription_;
  //get the QCD tree
  TFile fqcd(myname+".QCD.root");
  TTree* Tqcd = (TTree*) fqcd.Get("ABCDtree");

  TFile fsig(myname+".LM13.root"); //FIXME hard-coded!
  TTree* Tsignal = (TTree*) fsig.Get("ABCDtree");

  TFile fttbar(myname+".TTbarJets.root");
  TTree* Tttbar = (TTree*) fttbar.Get("ABCDtree");

  gROOT->cd();
  //create a diagnostic plot
  int nbinsx=50;
  int nbinsy=50;
  Habcd = new TH2D("Habcd","ABCD regions",nbinsx,xplotmin,xplotmax,nbinsy,yplotmin,yplotmax);
  Habcd->Sumw2();
  Habcd->SetXTitle(xvar_);
  Habcd->SetYTitle(yvar_);

  if (doFit_) {
    Hx_L = new TH1D("Hx_L","low y region",fitnum_,xl_low,xl_high);
    Hx_H = new TH1D("Hx_H","high y region",fitnum_,xl_low,xl_high);
    Hx_ratio = new TH1D("Hx_ratio","ratio of y regions",fitnum_,xl_low,xl_high);
    Hx_L->Sumw2();
    Hx_H->Sumw2();
    Hx_ratio->Sumw2();
  }

  //this calls the loop over events. also fills the fit histograms
  loopOverTree(Tqcd,true);
  if (dottbar)  loopOverTree(Tttbar,false);
  if (dosignal) loopOverTree(Tsignal,false);

  double totalerr =  error[0][1]*total[0][0]*total[0][0]*total[1][1]*total[1][1]/pow(total[0][1],4);
  totalerr += error[0][0]*pow(total[1][1]/total[0][1],2);
  totalerr += error[1][1]*pow(total[0][0]/total[0][1],2);
  totalerr = sqrt(totalerr);

  for (int i=0; i<2; i++) {
    for (int j=0; j<2; j++) {
      error[i][j] = sqrt(error[i][j]);
    }
  }

  cout<<"Actual number found in SR = "<<total[1][0]<<" +/- "<<error[1][0]<<endl;
  cout<<"Estimated number in SR    = "<<total[1][1]*total[0][0]/total[0][1]<<" +/- "<<totalerr<<endl;
  cout<<"----"<<endl;
  cout<<"Actual numbers in low MET regions = "<<total[0][0]<<" +/- "<<error[0][0]<<endl
      <<"                                  = "<<total[0][1]<<" +/- "<<error[0][1]<<endl;
  cout<<"Actual number found in high MET control region = "<<  total[1][1]<<" +/- "<<error[1][1]<<endl;

  if (doFit_) {
    Cfit->cd();
    /*
    Cfit->Divide(2,2);
    Cfit->cd(1);
    Hx_H->Draw();
    Cfit->cd(3);
    Hx_L->Draw();
    */
    Hx_ratio->Divide(Hx_H,Hx_L);
    expfunc_ = new TF1("expfunc","[0]*exp([1]*x)",xl_low,xl_high);
    expfunc_->SetParameter(0,Hx_ratio->GetBinContent(1));
    expfunc_->SetParameter(1,-2e-2); //completely empirical
    //    Cfit->cd(4);
    Hx_ratio->Fit(expfunc_,"R");
    loopOverTree(Tqcd,true,true);
    if (dottbar)  loopOverTree(Tttbar,false,true);
    if (dosignal) loopOverTree(Tsignal,false,true);
    extendedEstimateErr_ = sqrt(extendedEstimateErr_);
    cout<<"Extended SR estimate      = "<<extendedEstimate_<<" +/- "<<extendedEstimateErr_<<endl;
  }

  Cabcd->cd();
  Habcd->DrawCopy("COLZ"); //draw copy in order to avoid losing the picture when th2d goes out of scope
  TLine* vline_ll = new TLine(xl_low,yplotmin,xl_low,yplotmax);
  TLine* vline_lh = new TLine(xl_high,yplotmin,xl_high,yplotmax);
  vline_ll->SetLineWidth(2);
  vline_lh->SetLineWidth(2);
  vline_ll->Draw();
  vline_lh->Draw();
  TLine* hline_ll = new TLine(xplotmin,yl_low,xplotmax,yl_low);
  TLine* hline_lh = new TLine(xplotmin,yl_high,xplotmax,yl_high);
  hline_ll->SetLineWidth(2);
  hline_lh->SetLineWidth(2);
  hline_ll->Draw();
  hline_lh->Draw();

  TLine* vline_hl = new TLine(xh_low,yplotmin,xh_low,yplotmax);
  TLine* vline_hh = new TLine(xh_high,yplotmin,xh_high,yplotmax);
  vline_hl->SetLineWidth(2);
  vline_hh->SetLineWidth(2);
  vline_hl->Draw();
  vline_hh->Draw();
  TLine* hline_hl = new TLine(xplotmin,yh_low,xplotmax,yh_low);
  TLine* hline_hh = new TLine(xplotmin,yh_high,xplotmax,yh_high);
  hline_hl->SetLineWidth(2);
  hline_hh->SetLineWidth(2);
  hline_hl->Draw();
  hline_hh->Draw();

  double sry = flippedY_ ? yl_low : yh_low;
  TText sr(xh_low,sry,"Signal Region");
  sr.DrawClone();


}
