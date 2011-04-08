
double abError(double a, double aE, double b, double bE){

  double Err = 0;
  Err = sqrt( a*a*bE*bE + b*b*aE*aE);

  return Err;
}

void backgdSubError( double Ndata, double Nbackground, double Nbackground_err, double NQCD, double NQCD_err) {

  double Ndata_err = sqrt(Ndata);

  if (NQCD<=0) return;
  double S = (Ndata - Nbackground)/NQCD;

  double S_err = sqrt(Ndata_err*Ndata_err + Nbackground_err*Nbackground_err + NQCD_err*NQCD_err*pow(Nbackground-Ndata,2)/pow(NQCD,2))/NQCD;

  cout<<"S = "<<S<<" +/- "<<S_err<<endl;

}

void backgdSubError( double Ndata, double Nbackground, double Nbackground_err, double NQCD, double NQCD_err, double NQCDC, double NQCDC_err) {

  double Ndata_err = sqrt(Ndata);
  
  if (NQCD<=0) return;
  
  double S = (Ndata - Nbackground)/NQCD;
  double S_err = sqrt(Ndata_err*Ndata_err + Nbackground_err*Nbackground_err + NQCD_err*NQCD_err*pow(Nbackground-Ndata,2)/pow(NQCD,2))/NQCD;

  double SIG = S*NQCDC;
  double SIG_err = abError(S,S_err, NQCDC, NQCDC_err);

  cout<<"S = "<<S<<" +/- "<<S_err<<endl;
  cout<<"SIG = "<<SIG<<" +/- "<<SIG_err<<endl;

}

void nfromtree(double *n, TString file, TString cut) {

  TChain ch("ABCDtree");
  ch.Add(file);

  TH1D hht("hht","ht",1,0,1e6);
  hht.Sumw2();
  TString selection;
  selection.Form("weight*(%s)",cut.Data());
  cout<<selection<<endl;
  ch.Project("hht","HT",selection);
  n[0] = hht.GetBinContent(1);
  n[1] = hht.GetBinError(1);

}

void runit()
{

  /*
  TString qcdfile="/cu1/joshmt/ABCDtrees/11Mar30/ABCDtree.Baseline0_PF_JERbias_pfMEThigh_PFLepRA20e0mu_minDP_MuonCleaning_NoMET_NoDeltaPhi.ge0b.PythiaPUQCDFlat.root";
  //hack -- I have made a directory with symlinks to only the background I want
  TString smfile="SMnonQCD/*.root";
  TString datafile="/cu1/joshmt/ABCDtrees/11Mar30/ABCDtree.Baseline0_PF_pfMEThigh_PFLepRA20e0mu_minDP_MuonCleaning_NoMET_NoDeltaPhi.ge0b.data-0.root";
  */

  TString qcdfile="/cu1/kreis/ABCDtrees/36_Mar26/ABCDtree.Baseline0_PF_JERbias_pfMEThigh_PFLepRA20e0mu_minDP_MuonCleaning_NoMET_NoDeltaPhi.ge0b.PythiaPUQCD.root";
  TString smfile="/cu1/kreis/ABCDtrees/36_Mar26/SMnonQCD/*.root";
  TString datafile="/cu1/kreis/ABCDtrees/36_Mar26/ABCDtree.Baseline0_PF_JERbias_pfMEThigh_PFLepRA20e0mu_minDP_MuonCleaning_NoMET_NoDeltaPhi.ge0b.data-0.root";
 
  //  TString cut = "minDeltaPhiMET>0.3 &&  nbSSVM>=2";
  //  TString cut = "minDeltaPhiMET<0.3 && MET >=150 &&  nbSSVM>=1";
  //  TString cut = "minDeltaPhiMET>0.3 && MET>=100 && MET<150 && nbSSVM>=1";
  //  TString cut = "minDeltaPhiMET>0.3 && MET>=150 && nbSSVM==1";
  
  const int numCuts = 3;
  TString cutA[numCuts] = {"minDeltaPhiMET<0.3 && MET >=150 &&  nbSSVM==1", "minDeltaPhiMET<0.3 && MET >=150 &&  nbSSVM>=1", "minDeltaPhiMET<0.3 && MET >=150 &&  nbSSVM>=2"};
  TString cutPassA[numCuts] = {"minDeltaPhiMET>=0.3 && MET >=150 &&  nbSSVM==1", "minDeltaPhiMET>=0.3 && MET >=150 &&  nbSSVM>=1", "minDeltaPhiMET>=0.3 && MET >=150 &&  nbSSVM>=2"};
  

  double qcd[2], data[2], other[2], qcdPass[2];

  for(int i = 0; i<numCuts; i++){

    cout << cutA[i] << endl;

    nfromtree(data,datafile,cutA[i]);
    nfromtree(qcd,qcdfile,cutA[i]);
    nfromtree(other,smfile,cutA[i]);
    nfromtree(qcdPass,qcdfile,cutPassA[i]);

    cout<<"Data D: "   << data[0]    << " +/- "<<data[1]    <<endl;
    cout<<"QCD D: "    << qcd[0]     << " +/- "<<qcd[1]     <<endl;
    cout<<"nonQCD D: " << other[0]   << " +/- "<<other[1]   <<endl;
    cout<<"QCD C: "    << qcdPass[0] << " +/- "<<qcdPass[1] <<endl;
    
    backgdSubError(data[0],other[0],other[1],qcd[0],qcd[1], qcdPass[0], qcdPass[1]);
  }//for i

}//end runit()
