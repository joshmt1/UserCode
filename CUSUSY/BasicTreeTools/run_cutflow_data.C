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

const TString version = "V00-01-04/DATA/386";

void run_cutflow_data()
{

  //TString libname="basicLoop_C.so";
  TString computername = gSystem->Getenv("HOST");
  TString dir = "/cu1/joshmt/";
  if (computername =="JoshPC") {
    dir="~/data/";
    //  libname="basicLoop_C.dll";
  }
  //could also add CASTOR

  //gSystem->Load(libname);

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

    //if (!samplefiles.Contains("MultiJet")) continue;
    cout<<"About to add files in: "<<samplefiles<<endl;

    ch.Add(samplefiles);
    info.Add(samplefiles);
  }

  basicLoop looper(&ch,&info);
  
  looper.setCutScheme(basicLoop::kBaseline0);
  looper.setMETType(basicLoop::kpfMET);
  //  looper.setMETRange(basicLoop::kHigh); //signal region
  looper.setMETRange(basicLoop::kMedium); //50 - 100 GeV region
  looper.setJetType(basicLoop::kPF);
  looper.setLeptonType(basicLoop::kPFLeptons); //PF leptons
  //looper.setLeptonType(basicLoop::kNormal); //normal (not PF) leptons
  
  looper.setDPType(basicLoop::kminDP);
  //  looper.setIgnoredCut("cutDeltaPhi"); //was ignoring this cut for preliminary sections of the note
  
  looper.setBCut(3); //require 3 b tags so that we make the full cut flow table

  //looper.setMuonReq(1); //inverted muon veto

  looper.cutflow();
  //  looper.screendump();
  
}
