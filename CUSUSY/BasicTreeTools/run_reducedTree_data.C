#include "TString.h"
#include "TSystem.h"
#include "TChain.h"
#include "TObjArray.h"

#include "basicLoop.C"
/*
Usage:
root -b -l -q run_reducedTree.C++

*/
const TString version = "V00-03-01/DATA/387";
const TString ecalDeadCellPath = "/cu1/joshmt/ECALDeadCellNtuples/DATA/387"; //hard coded for dellcmscornell

void run_reducedTree_data()
{
  TString computername = gSystem->Getenv("HOST");
  TString dir = "/cu1/joshmt/"; //default is dellcmscornell
  if (computername =="JoshPC") {
    dir="~/data/";
    //  libname="basicLoop_C.dll";
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
  
  for (int ifile=0; ifile<nfiles; ifile++) {
    TString samplefiles = dirlist->At(ifile)->GetTitle();
    samplefiles+="/*.root";
    
    cout<<"About to add files in: "<<samplefiles<<endl;

    ch.Add(samplefiles);
    info.Add(samplefiles);
  }

  //fetch ecal dead cell info
  TFile fEcal(ecalDeadCellPath + "/deadCellFilterProfile.root"); //hadd files together into one
  TTree* ecalTree = fEcal.IsZombie() ? 0 : (TTree*) fEcal.Get("filter");

  basicLoop looper(&ch,&info,ecalTree);
  
  //important! this is where cuts are defined
  //looper.setCutScheme(basicLoop::kRA2); //usually this is kRA2
  looper.setCutScheme(basicLoop::kBaseline0); //usually this is kRA2
  looper.setMETType(basicLoop::kpfMET);
  looper.setLeptonType(basicLoop::kPFLeptonsRA2);
  looper.setJetType(basicLoop::kPF);
  looper.setCleaningType(basicLoop::kMuonEcalCleaning);
  
  //  looper.reducedTree("/cu2/joshmt/");  //go!
  looper.reducedTree("/home/joshmt/");  //go!

}
