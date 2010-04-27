/*
findAllWeights(lumi); //update weight files

makeSimpleTrees(); //makes the simple trees

*/

#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCut.h"

#include "TSystem.h"

#include "/afs/cern.ch/user/j/joshmt/root/util/HistHolder.h"
#include "/afs/cern.ch/user/j/joshmt/root/util/PlotUtil.h"

using namespace std;

#include "eventLoop.C"
// gROOT->LoadMacro("eventLoop.C");

  //======
  const UInt_t kqcd = 1;
  const UInt_t ksig = 2;
  const UInt_t kttbar = 4;
  const UInt_t kww = 3;
  const UInt_t kwz = 5;
  const UInt_t kzz = 7;
  const UInt_t kvqq = 6;
  const UInt_t kwjets = 8;
  const UInt_t kzjets = 9;

  const double pi=4*atan(1);
  //======
//  const TString path= "data/"; //for running on dellcmscornell
// const TString path= "tempdata/"; //for running on elsewhere

TString getPath() {
  TString path=  ( TString(gSystem->Getenv("HOST"))=="dellcmscornell" )
    ? "data/" : "tempdata/";

  return path;
}

void findWeight(TString sampleid, float lumi=100) { //pb^-1

  //  cout<<"[findWeight] "<<lumi<<"\t"<<sampleid<<endl;

  double sigma=0,ngen=0;
  //LM0
  if (sampleid=="LM0") {
    sigma = 38.93;
    ngen = 208629;
  }
  else if (sampleid=="LM1") {
    sigma = 4.888;
    ngen = 215853;
  }
  else if (sampleid=="LM2") {
    sigma = 0.603;
    ngen = 199120;
  }
  else if (sampleid=="LM3") {
    sigma = 3.438;
    ngen = 228875;
  }
  else if (sampleid=="LM4") {
    sigma = 1.879;
    ngen = 252265;
  }
  else if (sampleid=="LM5") {
    sigma = 0.4734;
    ngen = 238880;
  }
  else if (sampleid=="LM6") {
    sigma = 0.3104;
    ngen = 265150;
  }
  else if (sampleid=="LM7") {
    sigma = 1.209;
    ngen = 265923;
  }
  else if (sampleid=="LM8") {
    sigma = 0.73;
    ngen = 216625;
  }
  else if (sampleid=="LM9") { //LM9
    sigma = 7.134;
    ngen = 216722;
  }
  else if (sampleid=="LM10") {
    sigma = 0.04778;
    ngen = 210514;
  }
  else if (sampleid=="LM11") {
    sigma = 0.8236;
    ngen = 267260;
  }
  else if (sampleid=="LM12") {
    sigma = 4.414;
    ngen = 203713;
  }
  else if (sampleid=="LM13") { //LM13
    sigma = 6.899;
    ngen = 223144;
  }
  else if (sampleid=="ZZ") {
    sigma = 4.3;
    ngen = 145368;
  }
  else if (sampleid=="WZ") {
    sigma = 10.5;
    ngen = 114070;
  }
  else if (sampleid=="WW") {
    sigma = 28;
    ngen = 120280;
  }
  else if (sampleid=="WJets-madgraph") {
    sigma = 28000;
    ngen = 10068895;
  }
  else if (sampleid=="ZJets-madgraph") {
    sigma = 2800;
    ngen = 1084921;
  }
  else if (sampleid=="VqqJets-madgraph") {
    sigma = 36;
    ngen = 936242;
  }
  else if (sampleid=="TTBar") { //PYTHIA
    sigma = 165;//this is NLO; it is 94.3 LO
    ngen = 626610;
  }
  else {
    cout<<"Sample not defined!"<<endl;
    return;
  }

  TH1F  Hlumi("Hlumi","luminosity used in weight calculation (pb^-1)",1,0,1);
  Hlumi.SetBinContent(1,lumi);

  TH1D  Hweight("Hweight","weight as a function of pthat",1,0,1000000);
  Hweight.SetBinContent(1,lumi*sigma/ngen);
  
  TFile fw(getPath()+"weight."+sampleid+".root","RECREATE");
  Hlumi.Write();
  Hweight.Write();
  fw.Close();

}


void findQCDWeight(float lumi) {

  const int nptbins = 5; //this had better agree with the next line
  //  Float_t ptbins[]={0,15,30,80,1000000};
  Float_t ptbins[]={0,15,30,80,170,1000000};

  TH1D Hsigmaqcd("Hsigma","cross section of qcd events",nptbins,ptbins);
  TH1D HQCDweight("Hweight","weight of qcd events",nptbins,ptbins);
  // ======= QCD cross sections =======

  const double s_qcd15=8.762e+08;
  const double s_qcd30=6.041e+07;
  const double s_qcd80=9.238e+05;
  const double s_qcd170=25474.9;

  Hsigmaqcd.SetBinContent(1,0);
  Hsigmaqcd.SetBinContent(2,s_qcd15 - s_qcd30);
  Hsigmaqcd.SetBinContent(3,s_qcd30 - s_qcd80);
  Hsigmaqcd.SetBinContent(4,s_qcd80 - s_qcd170);
  Hsigmaqcd.SetBinContent(5,s_qcd170);
  // ==================================

  TH1F Hnqcd("Hnqcd","n of qcd events",nptbins,ptbins);
  TH1F Hnqcd15("Hnqcd15","n of qcd events",nptbins,ptbins);
  TH1F Hnqcd30("Hnqcd30","n of qcd events",nptbins,ptbins);
  TH1F Hnqcd80("Hnqcd80","n of qcd events",nptbins,ptbins);
  TH1F Hnqcd170("Hnqcd170","n of qcd events",nptbins,ptbins);

  const TString treedirectory="makeTopologyNtuple/tree";

  //big question -- can i trust that the trees have 100% of the events?
  //each root file also has an event counter histogram
  //in principle i should check against that

  TChain In_qcd170(treedirectory);
  In_qcd170.Add(getPath()+"QCD_Pt170/Ntuple_QCD_Pt170_*.root");

  TChain In_qcd80(treedirectory);
  In_qcd80.Add(getPath()+"Ntuple_QCD_Pt80_*.root");

  TChain In_qcd30(treedirectory);
  In_qcd30.Add(getPath()+"Ntuple_QCD_Pt30_*.root");

  TChain In_qcd15(treedirectory);
  In_qcd15.Add(getPath()+"Ntuple_QCD_Pt15_*.root");

  In_qcd15.Project("Hnqcd15","processPtHat");
  In_qcd30.Project("Hnqcd30","processPtHat");
  In_qcd80.Project("Hnqcd80","processPtHat");
  In_qcd170.Project("Hnqcd170","processPtHat");

  Hnqcd.Add(&Hnqcd15);
  Hnqcd.Add(&Hnqcd30);
  Hnqcd.Add(&Hnqcd80);
  Hnqcd.Add(&Hnqcd170);

  //  cout<<  Hnqcd.Integral()<<endl;
  //  cout<<  Hnqcd15.Integral() +Hnqcd30.Integral()+Hnqcd80.Integral()<<endl;

  for (int i=1;i<=nptbins; i++) {
    if (Hnqcd.GetBinContent(i) !=0) {
      HQCDweight.SetBinContent(i, lumi*Hsigmaqcd.GetBinContent(i)/Hnqcd.GetBinContent(i));
    }
    else {
      HQCDweight.SetBinContent(i,0);
    }
  }

  TH1F  Hlumi("Hlumi","luminosity used in weight calculation (pb^-1)",1,0,1);
  Hlumi.SetBinContent(1,lumi);
  
  TFile fw(getPath()+"weight.QCD.root","RECREATE");
  Hlumi.Write();
  Hnqcd.Write();
  Hsigmaqcd.Write();
  HQCDweight.Write();
  fw.Close();

}

void findAllWeights(float lumi=100) {
  findWeight("TTBar",lumi);

  findWeight("WJets-madgraph",lumi);
  findWeight("ZJets-madgraph",lumi);

  findWeight("VqqJets-madgraph",lumi);

  findWeight("LM9",lumi);
  findWeight("LM13",lumi);

  findWeight("WW",lumi);
  findWeight("WZ",lumi);
  findWeight("ZZ",lumi);

  findQCDWeight(lumi);

}

void makeSimpleTrees( int bjetcut, TString which="all" ) {

  //bjetcut = 0 means no cut
  //bjetcut = 1 or 2 means == 1 or == 2
  //bjetcut = 3      means >=3

  //write ascii or not
  const bool writeASCII=false;

  gSystem->Load("eventLoop_C.so");
  const TString treedirectory="makeTopologyNtuple/tree";

  //fixme should just do all lm samples in this function
  int LMsample=9;
  TString signalname="LM";
  signalname+=LMsample;
  TString signalfile; signalfile.Form("%sNtuple_LM%d-7TeV*.root",getPath().Data(),LMsample);
  TChain* In_sig = new TChain(treedirectory);
  In_sig->Add(signalfile);

  TChain* In_qcd80 = new TChain(treedirectory);
  In_qcd80->Add(getPath()+"Ntuple_QCD_Pt80_*.root");
  TChain* In_qcd15 =  new TChain(treedirectory);
  In_qcd15->Add(getPath()+"Ntuple_QCD_Pt15_*.root");
  TChain* In_qcd30 =  new TChain(treedirectory);
  In_qcd30->Add(getPath()+"Ntuple_QCD_Pt30_*.root");
  TChain* In_qcd170 =  new TChain(treedirectory);
  In_qcd170->Add(getPath()+"QCD_Pt170/Ntuple_QCD_Pt170_*.root");

  TChain* In_ttbar = new TChain(treedirectory);
  In_ttbar->Add(getPath()+"Ntuple_TTbar*.root");

  TChain* In_ww = new TChain(treedirectory);
  In_ww->Add(getPath()+"Ntuple_WW*.root");

  TChain* In_wz = new TChain(treedirectory);
  In_wz->Add(getPath()+"Ntuple_WZ*.root");

  TChain* In_zz = new TChain(treedirectory);
  In_zz->Add(getPath()+"Ntuple_ZZ*.root");

  TChain* In_zjets = new TChain(treedirectory);
  In_zjets->Add(getPath()+"Ntuple_ZJets*.root");

  TChain* In_wjets = new TChain(treedirectory);
  In_wjets->Add(getPath()+"Ntuple_WJets*.root");

  TChain* In_vqq = new TChain(treedirectory);
  In_vqq->Add(getPath()+"Ntuple_Vqq*.root");

  TString samplename;
  TString outpath=getPath()+"simpletrees/";

  TString nbs;
  nbs.Form(".%db",bjetcut);

  if (which.Contains("all") || which.Contains("sig")) {
    eventLoop Lsig(In_sig); 
    Lsig.writeASCII(writeASCII); Lsig.goodBJets_ = bjetcut;
    //don't forget the special true for signal
    Lsig.Loop(outpath+signalname+nbs+".root",getPath()+"weight."+signalname+".root",true);
  }

  if (which.Contains("all") || which.Contains("qcd80")) {
    eventLoop Lqcd80(In_qcd80); 
    Lqcd80.writeASCII(writeASCII);  Lqcd80.goodBJets_ = bjetcut;
    Lqcd80.Loop(outpath+"QCD80"+nbs+".root",getPath()+"weight.QCD.root");
  }

  if (which.Contains("all") || which.Contains("qcd30")) {
    eventLoop Lqcd30(In_qcd30); 
    Lqcd30.writeASCII(writeASCII);  Lqcd30.goodBJets_ = bjetcut;
    Lqcd30.Loop(outpath+"QCD30"+nbs+".root",getPath()+"weight.QCD.root");
  }

  if (which.Contains("all") || which.Contains("qcd15")) {
    eventLoop Lqcd15(In_qcd15); 
    Lqcd15.writeASCII(writeASCII); Lqcd15.goodBJets_ = bjetcut;
    Lqcd15.Loop(outpath+"QCD15"+nbs+".root",getPath()+"weight.QCD.root");
  }

  if (which.Contains("all") || which.Contains("qcd170")) {
    eventLoop Lqcd170(In_qcd170);
    Lqcd170.writeASCII(writeASCII); Lqcd170.goodBJets_ = bjetcut;
    Lqcd170.Loop(outpath+"QCD170"+nbs+".root",getPath()+"weight.QCD.root");
  }

  if (which.Contains("all") || which.Contains("ttbar")) {
    samplename="TTBar";
    eventLoop Lttbar(In_ttbar); 
    Lttbar.writeASCII(writeASCII); Lttbar.goodBJets_ = bjetcut;
    Lttbar.Loop(outpath+samplename+nbs+".root",getPath()+"weight."+samplename+".root");
  }

  if (which.Contains("all") || which.Contains("ww")) {
    samplename="WW";
    eventLoop Lww(In_ww); 
    Lww.writeASCII(writeASCII); Lww.goodBJets_ = bjetcut;
    Lww.Loop(outpath+samplename+nbs+".root",getPath()+"weight."+samplename+".root");
  }

  if (which.Contains("all") || which.Contains("wz")) {
    samplename="WZ";
    eventLoop Lwz(In_wz); 
    Lwz.writeASCII(writeASCII); Lwz.goodBJets_ = bjetcut;
    Lwz.Loop(outpath+samplename+nbs+".root",getPath()+"weight."+samplename+".root");
  }

  if (which.Contains("all") || which.Contains("zz")) {
    samplename="ZZ";
    eventLoop Lzz(In_zz); 
    Lzz.writeASCII(writeASCII); Lzz.goodBJets_ = bjetcut;
    Lzz.Loop(outpath+samplename+nbs+".root",getPath()+"weight."+samplename+".root");
  }

  if (which.Contains("all") || which.Contains("zjets")) {
    samplename="ZJets-madgraph";
    eventLoop Lzjets(In_zjets); 
    Lzjets.writeASCII(writeASCII);  Lzjets.goodBJets_ = bjetcut;
    Lzjets.Loop(outpath+samplename+nbs+".root",getPath()+"weight."+samplename+".root");
  }

  if (which.Contains("all") || which.Contains("wjets")) {
    samplename="WJets-madgraph";
    eventLoop Lwjets(In_wjets); 
    Lwjets.writeASCII(writeASCII); Lwjets.goodBJets_ = bjetcut;
    Lwjets.Loop(outpath+samplename+nbs+".root",getPath()+"weight."+samplename+".root");
  }

  if (which.Contains("all") || which.Contains("vqq")) {
    samplename="VqqJets-madgraph";
    eventLoop Lvqq(In_vqq); 
    Lvqq.writeASCII(writeASCII); Lvqq.goodBJets_ = bjetcut;
    Lvqq.Loop(outpath+samplename+nbs+".root",getPath()+"weight."+samplename+".root");
  }

}

//this seems kinda roundabout. but it works
//user provides a printf-style string for forming the cut
TCut formCut(TString formstring, Float_t value) {
  TString cut;
  cut.Form(formstring.Data(),value);
  TCut mycut = cut.Data();
  return mycut;
}

//these are defined globally
//so that plots will not go blank at the end of drawme()
TChain* sig =0;
TChain* qcd15 =0;
TChain* qcd30 =0;
TChain* qcd80 =0;
TChain* qcd170 =0;
TChain* ttbar =0;
TChain* vqq =0;
TChain* wjets =0;
TChain* zjets =0;
TChain* ww =0;
TChain* wz =0;
TChain* zz =0;

//needs some modification for new file names!

void searchme( UInt_t LMsample=9) {

  const  TString outpath=getPath()+"simpletrees/";

  TString samplename;
  TString signalname="LM";
  signalname+=LMsample;

  TFile fsig(outpath+signalname+".root");
   sig = (TChain*) fsig.Get("tree");

  TFile fqcd15(outpath+"QCD15.root");
   qcd15 = (TChain*) fqcd15.Get("tree");

  TFile fqcd30(outpath+"QCD30.root");
   qcd30 = (TChain*) fqcd30.Get("tree");

  TFile fqcd80(outpath+"QCD80.root");
   qcd80 = (TChain*) fqcd80.Get("tree");

  TFile fqcd170(outpath+"QCD170.root");
   qcd170 = (TChain*) fqcd170.Get("tree");

  samplename="TTBar";
  TFile fttbar(outpath+samplename+".root");
   ttbar = (TChain*) fttbar.Get("tree");

  samplename="VqqJets-madgraph";
  TFile fvqq(outpath+samplename+".root");
   vqq = (TChain*) fvqq.Get("tree");

  samplename="WJets-madgraph";
  TFile fwjets(outpath+samplename+".root");
   wjets = (TChain*) fwjets.Get("tree");

  samplename="ZJets-madgraph";
  TFile fzjets(outpath+samplename+".root");
   zjets = (TChain*) fzjets.Get("tree");

  samplename="WW";
  TFile fww(outpath+samplename+".root");
   ww = (TChain*) fww.Get("tree");

  samplename="WZ";
  TFile fwz(outpath+samplename+".root");
   wz = (TChain*) fwz.Get("tree");

  samplename="ZZ";
  TFile fzz(outpath+samplename+".root");
   zz = (TChain*) fzz.Get("tree");

  //-----------------------------------------------------
  gROOT->cd(); //important!
  HistHolder histo;

  //0 leptons is already enforced
  //primary vertex is missing from top ntuple (i think)
  //  TString basecut = "1"; //not using this anymore

  int nbins;
  float hmin,hmax;
  TString varname;
  // ========= make histograms ===========

  PlotUtil plots(&histo);
  //plots.setDebug(true);
  //arguments are: identifier, tree, name for humans, color code //last 2 are optional
  plots.addSample("sig",sig,signalname,ksig);
  plots.addSample("qcd",qcd15,"QCD",kqcd); //qcd is a 'magic' identifier!
  plots.addSample("qcd30",qcd30,"QCD1",ksig);
  plots.addSample("qcd80",qcd80,"QCD2",ksig);
  plots.addSample("qcd170",qcd170,"QCD3",ksig);
  plots.addSample("ttbar",ttbar,"TTBar",kttbar);
  plots.addSample("vqq",vqq,"Vqq+Jets",kvqq);
  plots.addSample("wjets",wjets,"W+Jets",kwjets);
  plots.addSample("zjets",zjets,"Z+Jets",kzjets);
  plots.addSample("ww",ww,"WW",kww);
  plots.addSample("wz",wz,"WZ",kwz);
  plots.addSample("zz",zz,"ZZ",kzz);

  varname="MET"; nbins=100; hmin=0; hmax=700;
  plots.createHistos(varname,nbins,hmin,hmax);
  histo.Sumw2();

  TString env_ID(gSystem->Getenv("JOBID"));
  Int_t jobid = env_ID.Atoi();

  TString outfilename;
  outfilename.Form("cutexploration.%d.root",jobid);

  TFile fout(outfilename,"RECREATE");
  Float_t signif,nsig,nqcd,nttbar,nvqq,nwjets,nzjets,nww,nwz,nzz;
  Float_t nsig_err,nqcd_err,nttbar_err,nvqq_err,nwjets_err,nzjets_err,nww_err,nwz_err,nzz_err;
  Float_t minpt1,minpt2,minpt3,minminDeltaPhi,minMET;
  Float_t minDeltaPhi1,minDeltaPhi2;
  Float_t minHT;//,maxmaxDeltaPhi;
  //  Float_t maxEMf, maxb1Eta;
  Float_t minST,maxDeltaPhiMPTMET;
  Float_t min_ngoodjets,nbtags;
  Float_t min_aplanarity;

  TTree cuttree("cuttree","cuttree");
  cuttree.Branch("min_ngoodjets",&min_ngoodjets,"min_ngoodjets");
  cuttree.Branch("nbtags",&nbtags,"nbtags");

  cuttree.Branch("minpt1",&minpt1,"minpt1");
  cuttree.Branch("minpt2",&minpt2,"minpt2");
  cuttree.Branch("minpt3",&minpt3,"minpt3");
  cuttree.Branch("minminDeltaPhi",&minminDeltaPhi,"minminDeltaPhi");
  cuttree.Branch("minDeltaPhi1",&minDeltaPhi1,"minDeltaPhi1");
  cuttree.Branch("minDeltaPhi2",&minDeltaPhi2,"minDeltaPhi2");

  cuttree.Branch("minMET",&minMET,"minMET");
  cuttree.Branch("minHT",&minHT,"minHT");
  //  cuttree.Branch("maxEMf",&maxEMf,"maxEMf");
  //  cuttree.Branch("maxmaxDeltaPhi",&maxmaxDeltaPhi,"maxmaxDeltaPhi");
  //  cuttree.Branch("maxb1Eta",&maxb1Eta,"maxb1Eta");
  cuttree.Branch("minST",&minST,"minST");
  cuttree.Branch("maxDeltaPhiMPTMET",&maxDeltaPhiMPTMET,"maxDeltaPhiMPTMET");
  cuttree.Branch("min_aplanarity",&min_aplanarity,"min_aplanarity");

  cuttree.Branch("signif",&signif,"signif");
  cuttree.Branch("nsig",&nsig,"nsig");
  cuttree.Branch("nqcd",&nqcd,"nqcd");
  cuttree.Branch("nttbar",&nttbar,"nttbar");
  cuttree.Branch("nvqq",&nvqq,"nvqq");
  cuttree.Branch("nwjets",&nwjets,"nwjets");
  cuttree.Branch("nzjets",&nzjets,"nzjets");
  cuttree.Branch("nww",&nww,"nww");
  cuttree.Branch("nwz",&nwz,"nwz");
  cuttree.Branch("nzz",&nzz,"nzz");

  cuttree.Branch("nsig_err",&nsig_err,"nsig_err");
  cuttree.Branch("nqcd_err",&nqcd_err,"nqcd_err");
  cuttree.Branch("nttbar_err",&nttbar_err,"nttbar_err");
  cuttree.Branch("nvqq_err",&nvqq_err,"nvqq_err");
  cuttree.Branch("nwjets_err",&nwjets_err,"nwjets_err");
  cuttree.Branch("nzjets_err",&nzjets_err,"nzjets_err");
  cuttree.Branch("nww_err",&nww_err,"nww_err");
  cuttree.Branch("nwz_err",&nwz_err,"nwz_err");
  cuttree.Branch("nzz_err",&nzz_err,"nzz_err");

  gROOT->cd(); //important!



  TString env_minMET(gSystem->Getenv("MET"));
  minMET = env_minMET.Atoi();

  TString env_minST(gSystem->Getenv("ST"));
  minST = env_minST.Atoi() * 0.1;

  TString env_minHT(gSystem->Getenv("HT"));
  minHT = env_minHT.Atoi();

  TString env_minPT1(gSystem->Getenv("PT1"));
  minpt1 = env_minPT1.Atoi();

  TString env_minPT2(gSystem->Getenv("PT2"));
  minpt2 = env_minPT2.Atoi();

  TString env_minPT3(gSystem->Getenv("PT3"));
  minpt3 = env_minPT3.Atoi();

  cout<<"Starting..."<<endl;

  for (nbtags = 2; nbtags <=3 ; nbtags++) {
    //this is a special cut, in that the functional form of the cut depends on the cut
    TCut cut_nbtags;
    if (nbtags==2) {
      cut_nbtags = formCut("nbSSVM == %f",nbtags);
    }
    else {
      cut_nbtags = formCut("nbSSVM >= %f",nbtags);
    }
    
    for (min_ngoodjets =3; min_ngoodjets<=4; min_ngoodjets++) {
      TCut cut_njets= formCut("ngoodjets >= %f",min_ngoodjets);
      cout<<int(min_ngoodjets)<<endl;      
      //  for ( minHT = 200 ; minHT <=400 ; minHT += 100) {
      TCut cut_ht = formCut("HT > %f",minHT);
      
      //pT
      //    for ( minpt1 = 100 ; minpt1 <= 200; minpt1+=25) {
      TCut cut_pt1 = formCut("jet_pT1 > %f",minpt1);
      
      //      for ( minpt2 = 50 ; minpt2 <= 150; minpt2+=25) {
      TCut cut_pt2 = formCut("jet_pT2 > %f",minpt2);
      
      //for ( minpt3 = 50 ; minpt3 <= 91; minpt3+=20) {
	TCut cut_pt3 = formCut("jet_pT3 > %f",minpt3);
	
// 	for ( maxEMf = 0.99 ; maxEMf <= 1; maxEMf += 0.01) {
// 	  TCut cut_emf1 = formCut("jet_EMf1 <= %f",maxEMf);
// 	  TCut cut_emf2 = formCut("jet_EMf2 <= %f",maxEMf);
// 	  TCut cut_emf3 = formCut("jet_EMf3 <= %f",maxEMf);

	for (min_aplanarity = 0; min_aplanarity <=0.021; min_aplanarity +=0.01) {
	  TCut cut_aplanarity=formCut("aplanarity >= %f",min_aplanarity);

	  for ( maxDeltaPhiMPTMET =1.5; maxDeltaPhiMPTMET<=3.5; maxDeltaPhiMPTMET+=1) {
	    TCut cut_dpmptmet=formCut("(3.14159-DeltaPhi_MPTLoose_MET) < %f",maxDeltaPhiMPTMET);
	    
	    //minDeltaPhi
	    for ( minminDeltaPhi = 0.0; minminDeltaPhi <=0.4 ; minminDeltaPhi+=0.1) {
	      TCut cut_minDeltaPhi = formCut("minDeltaPhi > %f",minminDeltaPhi);

	    for ( minDeltaPhi1 = 0.0; minDeltaPhi1 <=1.01 ; minDeltaPhi1+=0.5) {
	      TCut cut_minDeltaPhi1 = formCut("DeltaPhi1 >= %f",minDeltaPhi1);

	    for ( minDeltaPhi2 = 0.0; minDeltaPhi2 <=1.01 ; minDeltaPhi2+=0.5) {
	      TCut cut_minDeltaPhi2 = formCut("DeltaPhi2 >= %f",minDeltaPhi2);
	      
// 	      for (maxmaxDeltaPhi = 175 ; maxmaxDeltaPhi <=180 ; maxmaxDeltaPhi+=5) {
// 		TCut cut_maxDeltaPhi = formCut("maxDeltaPhi <= (%f*3.14159/180.0)",maxmaxDeltaPhi);
		
		//		for (maxb1Eta = 2; maxb1Eta <=2.41; maxb1Eta+=0.2) {
		  //i think cuts use abs instead of fabs
		//		  TCut cut_maxb1eta = formCut("abs(bjet_eta1) <= %f",maxb1Eta);
		  
		    //		    for (minST = 0; minST<=0.3; minST+=0.1) {
		  TCut cut_st=formCut("ST >= %f",minST);
		  
		  //MET
		  //		      for (  minMET = 100; minMET <= 225 ; minMET+=25) {
		  
		  TCut cut_MET = formCut("MET > %f",minMET); 
		  
		  //put together cuts
		  //as usual, the lack of STL here is annoying....
		  TCut othercut = cut_pt1 && cut_pt2 && cut_pt3 
		    && cut_minDeltaPhi //&& cut_maxDeltaPhi 
		    && cut_minDeltaPhi1 && cut_minDeltaPhi2
		    && cut_MET 
		    && cut_ht 
		    //		    && cut_emf1 && cut_emf2 && cut_emf3
		    //		    && cut_maxb1eta 
		    && cut_dpmptmet &&cut_st
		    && cut_aplanarity
		    && cut_njets &&cut_nbtags;
		  TString cut;
		  //		  cut.Form("((%s) && (%s))*eventweight",basecut.Data(),othercut.GetTitle());
		  cut.Form("(%s)*eventweight",othercut.GetTitle());
  
		  //comment out for batch
		  //		  		  cout<<cut<<endl;
		  cout<<".";

  histo.Reset(); //reset histograms

  plots.fillHistos(varname,cut);
  plots.addQCD(varname);

  // ========= draw histograms ===========

  /* don't need to draw for this exercise
     gStyle->SetOptLogy();
     TCanvas* c1 = new TCanvas("c1","c1");

  TLegend* leg = new TLegend(0.7,0.7,0.95,0.95);
  leg->SetFillColor(0);
  plots.fillLegend(leg); //automatically sets the colors

  plots.drawPlots();
  leg->Draw();
  */

  //in principle i can think of some clever way to get rid of this repetitive code
  nsig = histo["H"+varname+"_sig"]->Integral();
   nsig_err = plots.ErrorOnIntegral( histo["H"+varname+"_sig"] );

  Double_t nbkg=0,nbkg_err=0;

   nqcd = histo["H"+varname+"_qcd"]->Integral();
   nqcd_err = plots.ErrorOnIntegral(histo["H"+varname+"_qcd"]);
  nbkg+=nqcd;
  nbkg_err+=nqcd_err*nqcd_err;

   nttbar = histo["H"+varname+"_ttbar"]->Integral();
   nttbar_err = plots.ErrorOnIntegral(histo["H"+varname+"_ttbar"]);
  nbkg+=nttbar;
  nbkg_err+=nttbar_err*nttbar_err;

   nvqq = histo["H"+varname+"_vqq"]->Integral();
   nvqq_err = plots.ErrorOnIntegral(histo["H"+varname+"_vqq"]);
  nbkg+=nvqq;
  nbkg_err+=nvqq_err*nvqq_err;

   nwjets = histo["H"+varname+"_wjets"]->Integral();
   nwjets_err = plots.ErrorOnIntegral(histo["H"+varname+"_wjets"]);
  nbkg+=nwjets;
  nbkg_err+=nwjets_err*nwjets_err;

   nzjets = histo["H"+varname+"_zjets"]->Integral();
   nzjets_err = plots.ErrorOnIntegral(histo["H"+varname+"_zjets"]);
  nbkg+=nzjets;
  nbkg_err+=nzjets_err*nzjets_err;

   nww = histo["H"+varname+"_ww"]->Integral();
   nww_err = plots.ErrorOnIntegral(histo["H"+varname+"_ww"]);
  nbkg+=nww;
  nbkg_err+=nww_err*nww_err;

   nwz = histo["H"+varname+"_wz"]->Integral();
   nwz_err = plots.ErrorOnIntegral(histo["H"+varname+"_wz"]);
  nbkg+=nwz;
  nbkg_err+=nwz_err*nwz_err;

   nzz = histo["H"+varname+"_zz"]->Integral();
   nzz_err = plots.ErrorOnIntegral(histo["H"+varname+"_zz"]);
  nbkg+=nzz;
  nbkg_err+=nzz_err*nzz_err;

  nbkg_err = sqrt(nbkg_err);

  const  TString pm=" +/- ";
   signif = (nsig+nbkg>0) ? nsig/sqrt(nsig+nbkg) : 0;

   //  cout<<cut_pt1.GetTitle()<<" ; ";
   //  cout<<cut_pt2.GetTitle()<<" ; ";
   //  cout<<cut_pt3.GetTitle()<<" ; ";
   //  cout<<cut_minDeltaPhi.GetTitle()<<" ; ";
   //  cout<<cut_MET.GetTitle()<<" ; ";

   //comment out for batch
   //  cout<<signif<<" sigma ; Signal = "<<nsig<<pm<<nsig_err<<"\t"<<"Background = "<<nbkg<<pm<<nbkg_err<<endl;
  cuttree.Fill();
  
	    } //nbtags
	    } //ngoodjets
	    //		      } //MET
	    //		    }  //ST
	    //		  } //HT
	    //		} //pt1
	    //	      }  //pt2
	    //  } //pt3
	    //	  }  //EMf
	    //	} //b1eta
	    // }
	  }
	}
      }
    }
  }
  cout<<"about to be done!"<<endl;
  fout.cd();
  cuttree.Write();
  fout.Close();
}



void drawme(TString plotvar="MET", int nb=2, UInt_t LMsample=9) {

  TString nbs;
  nbs.Form(".%db",nb);

  if (nb==99) nbs="*";

  const  TString outpath=getPath()+"simpletrees/";

  TString samplename;
  TString signalname="LM";
  signalname+=LMsample;
  /*
  TFile fsig(outpath+signalname+nbs+".root");
   sig = (TChain*) fsig.Get("tree");

  TFile fqcd15(outpath+"QCD15"+nbs+".root");
   qcd15 = (TChain*) fqcd15.Get("tree");

  TFile fqcd30(outpath+"QCD30"+nbs+".root");
   qcd30 = (TChain*) fqcd30.Get("tree");

  TFile fqcd80(outpath+"QCD80"+nbs+".root");
   qcd80 = (TChain*) fqcd80.Get("tree");

  TFile fqcd170(outpath+"QCD170"+nbs+".root");
   qcd170 = (TChain*) fqcd170.Get("tree");

  samplename="TTBar";
  TFile fttbar(outpath+samplename+nbs+".root");
   ttbar = (TChain*) fttbar.Get("tree");

  samplename="VqqJets-madgraph";
  TFile fvqq(outpath+samplename+nbs+".root");
   vqq = (TChain*) fvqq.Get("tree");

  samplename="WJets-madgraph";
  TFile fwjets(outpath+samplename+nbs+".root");
   wjets = (TChain*) fwjets.Get("tree");

  samplename="ZJets-madgraph";
  TFile fzjets(outpath+samplename+nbs+".root");
   zjets = (TChain*) fzjets.Get("tree");

  samplename="WW";
  TFile fww(outpath+samplename+nbs+".root");
   ww = (TChain*) fww.Get("tree");

  samplename="WZ";
  TFile fwz(outpath+samplename+nbs+".root");
   wz = (TChain*) fwz.Get("tree");

  samplename="ZZ";
  TFile fzz(outpath+samplename+nbs+".root");
   zz = (TChain*) fzz.Get("tree");
  */
  ///


  TChain* sig = new TChain("tree");
  sig->Add(outpath+signalname+nbs+".root");

  TChain* qcd15 = new TChain("tree");
  qcd15->Add(outpath+"QCD15"+nbs+".root");

  TChain* qcd30 = new TChain("tree");
  qcd30->Add(outpath+"QCD30"+nbs+".root");

  TChain* qcd80 = new TChain("tree");
  qcd80->Add(outpath+"QCD80"+nbs+".root");

  TChain* qcd170 = new TChain("tree");
  qcd170->Add(outpath+"QCD170"+nbs+".root");

  samplename="TTBar";
  TChain* ttbar = new TChain("tree");
  ttbar->Add(outpath+samplename+nbs+".root");

  samplename="VqqJets-madgraph";
  TChain* vqq = new TChain("tree");
  vqq->Add(outpath+samplename+nbs+".root");

  samplename="WJets-madgraph";
  TChain* wjets = new TChain("tree");
  wjets->Add(outpath+samplename+nbs+".root");

  samplename="ZJets-madgraph";
  TChain* zjets = new TChain("tree");
  zjets->Add(outpath+samplename+nbs+".root");

  samplename="WW";
  TChain* ww = new TChain("tree");
  ww->Add(outpath+samplename+nbs+".root");

  samplename="WZ";
  TChain* wz = new TChain("tree");
  wz->Add(outpath+samplename+nbs+".root");

  samplename="ZZ";
  TChain* zz = new TChain("tree");
  zz->Add(outpath+samplename+nbs+".root");


  //-----------------------------------------------------
  gROOT->cd(); //important!
  HistHolder histo;

  //0 leptons is already enforced
  //primary vertex is missing from top ntuple (i think)
  //  TString basecut = "1"; //not using this anymore

  int nbins;
  float hmin,hmax;
  TString varname;
  // ========= make histograms ===========

  PlotUtil plots(&histo);
  //plots.setDebug(true);
  //arguments are: identifier, tree, name for humans, color code //last 2 are optional
  plots.addSample("sig",sig,signalname,ksig);
  plots.addSample("qcd",qcd15,"QCD",kqcd); //qcd is a 'magic' identifier!
  plots.addSample("qcd30",qcd30,"QCD1",ksig);
  plots.addSample("qcd80",qcd80,"QCD2",ksig);
  plots.addSample("qcd170",qcd170,"QCD3",ksig);
  plots.addSample("ttbar",ttbar,"TTBar",kttbar);
  plots.addSample("vqq",vqq,"Vqq+Jets",kvqq);
  plots.addSample("wjets",wjets,"W+Jets",kwjets);
  plots.addSample("zjets",zjets,"Z+Jets",kzjets);
  plots.addSample("ww",ww,"WW",kww);
  plots.addSample("wz",wz,"WZ",kwz);
  plots.addSample("zz",zz,"ZZ",kzz);

  varname="MET"; nbins=100; hmin=0; hmax=700;
  plots.createHistos(varname,nbins,hmin,hmax);

  varname="HT"; nbins=100; hmin=0; hmax=1000;
  plots.createHistos(varname,nbins,hmin,hmax);

  varname="aplanarity"; nbins=100; hmin=0; hmax=0.1;
  plots.createHistos(varname,nbins,hmin,hmax);

  varname="jet_pT1"; nbins=50; hmin=0; hmax=1000;
  plots.createHistos(varname,nbins,hmin,hmax);

  varname="jet_pT2"; nbins=50; hmin=0; hmax=500;
  plots.createHistos(varname,nbins,hmin,hmax);

  varname="jet_pT3"; nbins=50; hmin=0; hmax=250;
  plots.createHistos(varname,nbins,hmin,hmax);

  varname="minDeltaPhi"; nbins=50; hmin=0; hmax=pi;
  plots.createHistos(varname,nbins,hmin,hmax);

  varname="DeltaPhi1"; nbins=50; hmin=0; hmax=pi;
  plots.createHistos(varname,nbins,hmin,hmax);

  varname="DeltaPhi2"; nbins=50; hmin=0; hmax=pi;
  plots.createHistos(varname,nbins,hmin,hmax);

  varname="DeltaPhi_MPTLoose_MET"; nbins=50; hmin=0; hmax=pi;
  plots.createHistos(varname,nbins,hmin,hmax);

  histo.Sumw2();

  //  Float_t signif,nsig,nqcd,nttbar,nvqq,nwjets,nzjets,nww,nwz,nzz;
  //  Float_t nsig_err,nqcd_err,nttbar_err,nvqq_err,nwjets_err,nzjets_err,nww_err,nwz_err,nzz_err;
  Float_t minpt1=100,minpt2=50,minpt3=50,minminDeltaPhi=0.0,minMET=150;
  Float_t minDeltaPhi1=0.5,minDeltaPhi2=0;
  Float_t minHT=400;//,maxmaxDeltaPhi;
  Float_t maxEMf=1, maxb1Eta=3;
  Float_t minST=0,maxDeltaPhiMPTMET=3.5;
  Float_t min_ngoodjets=3,nbtags=nb;
  Float_t min_aplanarity=0;

  if (nb == 2) {
    minminDeltaPhi=0.1;
    minDeltaPhi1=0;
    maxDeltaPhiMPTMET=1.5;
    min_aplanarity=0.02;
  }
  else if (nb==99) {
   minpt1=100; minpt2=50; minpt3=50; minminDeltaPhi=0.0; minMET=150;
   minDeltaPhi1=0; minDeltaPhi2=0;
   minHT=400;//; maxmaxDeltaPhi;
   maxEMf=1;  maxb1Eta=3;
   minST=0; maxDeltaPhiMPTMET=3.5;
   min_ngoodjets=3; nbtags=nb;
   min_aplanarity=0;
  }

  gROOT->cd(); //important!

  cout<<"Starting..."<<endl;

    //this is a special cut, in that the functional form of the cut depends on the cut
    TCut cut_nbtags;
    if (nbtags==2) {
      cut_nbtags = formCut("nbSSVM == %f",nbtags);
    }
    else if (nbtags==99) {
      cut_nbtags = formCut("nbSSVM >= %f",2);
    }
    else {
      cut_nbtags = formCut("nbSSVM >= %f",nbtags);
    }
    

    TCut cut_njets= formCut("ngoodjets >= %f",min_ngoodjets);
    
    TCut cut_ht = formCut("HT > %f",minHT);
    
    TCut cut_pt1 = formCut("jet_pT1 > %f",minpt1);
    
    TCut cut_pt2 = formCut("jet_pT2 > %f",minpt2);
    
    TCut cut_pt3 = formCut("jet_pT3 > %f",minpt3);
    
    TCut cut_aplanarity=formCut("aplanarity >= %f",min_aplanarity);

    TCut cut_dpmptmet=formCut("(3.14159-DeltaPhi_MPTLoose_MET) < %f",maxDeltaPhiMPTMET);
    
    TCut cut_emf1 = formCut("jet_EMf1 <= %f",maxEMf);
    TCut cut_emf2 = formCut("jet_EMf2 <= %f",maxEMf);
    TCut cut_emf3 = formCut("jet_EMf3 <= %f",maxEMf);
    TCut cut_maxb1eta = formCut("abs(bjet_eta1) <= %f",maxb1Eta);

    //minDeltaPhi
    TCut cut_minDeltaPhi = formCut("minDeltaPhi > %f",minminDeltaPhi);

    TCut cut_minDeltaPhi1 = formCut("DeltaPhi1 >= %f",minDeltaPhi1);
    
    TCut cut_minDeltaPhi2 = formCut("DeltaPhi2 >= %f",minDeltaPhi2);
    
    TCut cut_st=formCut("ST >= %f",minST);
    
    //MET
    TCut cut_MET = formCut("MET > %f",minMET); 

    //gonna have to hack this in one variable ata  time.
    if (plotvar=="HT") cut_ht ="1";
    else if (plotvar=="DeltaPhi_MPTLoose_MET") cut_dpmptmet="1";
    else if (plotvar=="MET") cut_MET="1";
    else if (plotvar=="aplanarity") cut_aplanarity="1";
    else if (plotvar=="jet_pT1") cut_pt1="1";
    else if (plotvar=="jet_pT2") cut_pt2="1";
    else if (plotvar=="jet_pT3") cut_pt3="1";
    else if (plotvar=="minDeltaPhi") cut_minDeltaPhi="1";
    else if (plotvar=="DeltaPhi1") cut_minDeltaPhi1="1";
    else if (plotvar=="DeltaPhi2") cut_minDeltaPhi2="1";
		  
    //put together cuts
    //as usual, the lack of STL here is annoying....
    TCut othercut = cut_pt1 && cut_pt2 && cut_pt3 
      && cut_minDeltaPhi //&& cut_maxDeltaPhi 
      && cut_minDeltaPhi1 && cut_minDeltaPhi2
      && cut_MET 
      && cut_ht 
      && cut_emf1 && cut_emf2 && cut_emf3
      && cut_maxb1eta 
      && cut_dpmptmet &&cut_st
      && cut_aplanarity
      && cut_njets &&cut_nbtags;
    TString cut;
    //		  cut.Form("((%s) && (%s))*eventweight",basecut.Data(),othercut.GetTitle());
    cut.Form("(%s)*eventweight",othercut.GetTitle());
  
    histo.Reset(); //reset histograms
    
    TString drawcommand="";
    if (plotvar=="DeltaPhi_MPTLoose_MET") drawcommand="(3.14159-DeltaPhi_MPTLoose_MET)";
    plots.fillHistos(plotvar,cut,drawcommand);
    plots.addQCD(plotvar);
    
  // ========= draw histograms ===========
    
    // don't need to draw for this exercise
    //gStyle->SetOptLogy();
    TCanvas* c1 = new TCanvas("c1","c1");

    TLegend* leg = new TLegend(0.7,0.7,0.95,0.95);
    leg->SetFillColor(0);
    plots.fillLegend(leg); //automatically sets the colors
    
    plots.drawPlots();
    leg->Draw();
    
}
