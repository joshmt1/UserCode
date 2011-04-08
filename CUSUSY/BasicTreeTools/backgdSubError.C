void backgdSubError( double Ndata, double Nbackground, double Nbackground_err, double NQCD, double NQCD_err) {

  double Ndata_err = sqrt(Ndata);

  if (NQCD<=0) return;
  double S = (Ndata - Nbackground)/NQCD;

  double S_err = sqrt(Ndata_err*Ndata_err + Nbackground_err*Nbackground_err + NQCD_err*NQCD_err*pow(Nbackground-Ndata,2)/pow(NQCD,2))/NQCD;

  cout<<"S = "<<S<<" +/- "<<S_err<<endl;

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

  TString qcdfile="/cu1/joshmt/ABCDtrees/11Mar30/ABCDtree.Baseline0_PF_JERbias_pfMEThigh_PFLepRA20e0mu_minDP_MuonCleaning_NoMET_NoDeltaPhi.ge0b.PythiaPUQCDFlat.root";
  //hack -- I have made a directory with symlinks to only the background I want
  TString smfile="SMnonQCD/*.root";
  TString datafile="/cu1/joshmt/ABCDtrees/11Mar30/ABCDtree.Baseline0_PF_pfMEThigh_PFLepRA20e0mu_minDP_MuonCleaning_NoMET_NoDeltaPhi.ge0b.data-0.root";

  //  TString cut = "minDeltaPhiMET>0.3 &&  nbSSVM>=2";
  TString cut = "minDeltaPhiMET<0.3 && MET >=150 &&  nbSSVM>=1";
  //  TString cut = "minDeltaPhiMET>0.3 && MET>=100 && MET<150 && nbSSVM>=1";
  //  TString cut = "minDeltaPhiMET>0.3 && MET>=150 && nbSSVM==1";

  double qcd[2], data[2], other[2];

  nfromtree(data,datafile,cut);
  nfromtree(qcd,qcdfile,cut);
  nfromtree(other,smfile,cut);
  cout<<data[0]<< " +/- "<<data[1]<<endl;
  cout<<qcd[0]<< " +/- "<<qcd[1]<<endl;
  cout<<other[0]<< " +/- "<<other[1]<<endl;

  backgdSubError(data[0],other[0],other[1],qcd[0],qcd[1]);

}
