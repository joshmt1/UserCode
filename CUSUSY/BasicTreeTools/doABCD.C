//#include "TFile.h"

TCanvas * Cabcd;
void doABCD() {

  Cabcd=new TCanvas("Cabcd");
  //make a loop over events and pull events from tree

  //count events in each box for naive ABCD approach

  //should perhaps add a way to add some cuts on the other tree variable
  //this way,
  //eventually I can do a study with multiple methods that all cover the same
  //datasets

  //can do a fit to the low MET ratio, then a second loop over events
  //to get an extended estimate (Ben/RA2 approach)

  //could even contemplate setting up a 2D fit for uncorrelated variables

  //add option to handle QCD tree plus 'background' trees

  const TString dir="";

  //define the boundaries of the different regions
  const TString xvar = "MET";
  const double xl_low = 40; 
  const double xl_high = 120;

  const double xh_low = 150; 
  const double xh_high = 1e6;

  
  const TString yvar = "minDeltaRbj";
  const double yl_low = 0.5;
  const double yl_high = 2;

  const double yh_low = 2; 
  const double yh_high = 50;
  
  /*
  const TString yvar = "minDeltaPhi";
  const double yl_low = 0;
  const double yl_high = 0.3;

  const double yh_low = 0.3; 
  const double yh_high = 4;
  */
  //-----------------------------------------

  double xplotmin=0.9*xl_low;
  double xplotmax=400;
  double yplotmin=0.9*yl_low;
  double yplotmax=6;

  bool flippedY=false;
  if ( yvar=="minDeltaRbj") flippedY=true;
  else if (yvar=="minDeltaPhi") flippedY=false;
  else {
    cout<<"Can't find the ybar in my list!"<<endl;
    return;
  }

  double xval;
  double yval;
  double weight;

  //get the QCD tree
  TFile fqcd(dir+"ABCDtree.QCD.root");
  TTree* Tqcd = (TTree*) fqcd.Get("ABCDtree");
  Tqcd->SetBranchAddress(xvar.Data(),&xval);
  Tqcd->SetBranchAddress(yvar.Data(),&yval);
  Tqcd->SetBranchAddress("weight",&weight);

  double total[2][2]={0};
  double error[2][2]={0};

  //create a diagnostic plot
  int nbinsx=50;
  int nbinsy=50;
  TH2D Habcd("Habcd","ABCD regions",nbinsx,xplotmin,xplotmax,nbinsy,yplotmin,yplotmax);
  Habcd.Sumw2();
  Habcd.SetXTitle(xvar);
  Habcd.SetYTitle(yvar);

  int y0=1;
  int y1=0;
  if (flippedY) {
    y0=0;
    y1=1;
  }

  for (Long64_t iqcd = 0; iqcd<Tqcd->GetEntries() ; iqcd++) {
    Tqcd->LoadTree(iqcd);
    Tqcd->GetEntry(iqcd);

    Habcd.Fill(xval,yval,weight);

    int xbinid=-1;
    if ((xval >= xl_low) && (xval < xl_high) )      xbinid=0;
    else if ((xval >= xh_low) && (xval < xh_high) ) xbinid=1;

    if (xbinid== -1) continue; //not in one of the bins

    int ybinid=-1;
    if ((yval >= yl_low) && (yval < yl_high) )      ybinid=y0;
    else if ((yval >= yh_low) && (yval < yh_high) ) ybinid=y1;    

    if (ybinid== -1) continue; //not in one of the bins
    
    total[xbinid][ybinid] += weight;
    error[xbinid][ybinid] += weight*weight;
    
  }
  //first let's ignore the different ABCD and just do the minDeltaRbj configuration
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

  Habcd.DrawCopy("COLZ"); //draw copy in order to avoid losing the picture when th2d goes out of scope
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

  double sry = flippedY ? yl_low : yh_low;
  TText sr(xh_low,sry,"Signal Region");
  sr.DrawClone();

}
