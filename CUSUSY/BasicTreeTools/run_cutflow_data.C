#include "TString.h"
#include "TSystem.h"
#include "TChain.h"
#include "TObjArray.h"

#include "basicLoop.C"
/*
these are not needed (assuming this macro is compiled) because of the include above
(Not a solution that I love, but it works)

.L basicLoop.C++

gSystem->Load("basicLoop_C.so");

*/

const TString version = "V00-03-01/DATA/387";
const TString ecalDeadCellPath = "/cu1/joshmt/ECALDeadCellNtuples/DATA/387"; //hard coded for dellcmscornell

void run_cutflow_data()
{

  //TString libname="basicLoop_C.so";
  TString computername = gSystem->Getenv("HOST");
  TString dir = "/cu1/joshmt/";
  if (computername =="JoshPC") {
    dir="~/data/";
  }

  dir += "BasicNtuples/";
  dir += version; dir+="/";
  TChain dummy("dummy");
  TString dirs = dir; dirs+="*";
  dummy.Add(dirs);
  TObjArray* dirlist = dummy.GetListOfFiles();
  int nfiles=dirlist->GetEntries();

  TChain ch("BasicTreeMaker/tree");
  TChain info("BasicTreeMaker/infotree");

  //assemble the whole thing into one TChain!

  for (int ifile=0; ifile<nfiles; ifile++) {
    TString samplefiles = dirlist->At(ifile)->GetTitle();
    
    samplefiles+="/*.root";

    //if (!samplefiles.Contains("JetMETTau")) continue;
    cout<<"About to add files in: "<<samplefiles<<endl;

    ch.Add(samplefiles);
    info.Add(samplefiles);
  }

  //for data, let's just hadd the whole ecal dead cell stuff into one tree
  TFile fEcal(ecalDeadCellPath + "/deadCellFilterProfile.root");
  TTree* ecalTree = fEcal.IsZombie() ? 0 : (TTree*) fEcal.Get("filter");

  basicLoop looper(&ch,&info,ecalTree);
  //looper.setSpecialCutDescription("NewMuonCleaning");

  looper.setCutScheme(basicLoop::kBaseline0);
  looper.setMETType(basicLoop::kpfMET);
  looper.setMETRange(basicLoop::kHigh); //signal region
  //looper.setMETRange(basicLoop::kMedhigh);
  looper.setJetType(basicLoop::kPF);
  looper.setLeptonType(basicLoop::kPFLeptonsRA2); //PF leptons
  looper.setDPType(basicLoop::kminDP);
  looper.setCleaningType(basicLoop::kMuonCleaning);

  looper.setBCut(3); //require 3 b tags so that we make the full cut flow table

  //looper.setEleReq(1); //inverted muon veto

  looper.cutflow();
  //  looper.screendump();
  
}
