#include "TString.h"
#include "TSystem.h"
#include "TChain.h"
#include "TObjArray.h"

#include "basicLoop.C"

const TString version = "V00-01-01/DATA/38X";

void run_screendump_data()
{

  //TString libname="basicLoop_C.so";
  TString computername = gSystem->Getenv("HOST");
  TString dir = "/cu1/joshmt/";
  if (computername =="JoshPC") {
    dir="~/data/";
    //  libname="basicLoop_C.dll";
  }
  //could also add CASTOR

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

    if (!samplefiles.Contains("JetMET-")) continue;
    cout<<"About to add files in: "<<samplefiles<<endl;

    ch.Add(samplefiles);
    info.Add(samplefiles);
  }

  basicLoop looper(&ch,&info);
  looper.setCutScheme(basicLoop::kSync1);
  looper.setMETType(basicLoop::kpfMET);
  looper.setJetType(basicLoop::kPF);
  looper.setDPType(basicLoop::kDPSync1);
  
  looper.setBCut(3); //require 3 b tags so that we make the full cut flow table

  looper.screendump();
  
}
