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

const TString version = "V00-01-01"; //01 for most samples

//const TString extrapath = "SUSYPATv8_363"; //pass an empty string unless you need something else
const TString extrapath = "";

void run_cutflow()
{

  //TString libname="basicLoop_C.so";
  TString computername = gSystem->Getenv("HOST");
  TString dir = "/cu1/joshmt/";
  if (computername =="JoshPC") {
    dir="~/data/";
  }
  //could also add CASTOR

  dir += "BasicNtuples/";
  if (extrapath!="") { dir += extrapath; dir += "/";}
  dir += version; dir+="/";
  TChain dummy("dummy");
  TString dirs = dir; dirs+="*";
  dummy.Add(dirs);
  TObjArray* dirlist = dummy.GetListOfFiles();
  int nfiles=dirlist->GetEntries();

  for (int ifile=0; ifile<nfiles; ifile++) {
    TString samplefiles = dirlist->At(ifile)->GetTitle();
    //for the first iteration I had this stupid subdirectory in there
    if (version=="V00-00-01")  samplefiles+="/Spring10/*.root";
    else                       samplefiles+="/*.root";

    cout<<"About to start on files: "<<samplefiles<<endl;

    if (samplefiles.Contains("DATA")) continue; //skip data (use run_cutflow_data.C)
    //if (!(samplefiles.Contains("MSSM") )) continue; //hack to skip some samples
    
    TChain ch("BasicTreeMaker/tree");
    TChain info("BasicTreeMaker/infotree");
    ch.Add(samplefiles);
    info.Add(samplefiles);
    basicLoop looper(&ch,&info);

    /* sync exercise settings
    looper.setCutScheme(basicLoop::kSync1);
    looper.setMETType(basicLoop::kpfMET);
    looper.setJetType(basicLoop::kPF);
    looper.setLeptonType(basicLoop::kNormal);
    looper.setDPType(basicLoop::kDPSync1);
    */
    looper.setCutScheme(basicLoop::kBaseline0);
    looper.setMETType(basicLoop::kpfMET);
    looper.setMETRange(basicLoop::kHigh); //signal region
    //looper.setMETRange(basicLoop::kMedium); //50-100 GeV region
    looper.setJetType(basicLoop::kPF);
    looper.setLeptonType(basicLoop::kPFLeptons);
    //looper.setLeptonType(basicLoop::kNormal);

    looper.setDPType(basicLoop::kminDP);
    //    looper.setIgnoredCut("cutDeltaPhi"); //ignore this cut for preliminary sections of the note

    looper.setBCut(3); //require 3 b tags so that we make the full cut flow table

    //    looper.setMuonReq(1); //inverted muon veto

    looper.cutflow();
  }


}
