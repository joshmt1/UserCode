#include "TString.h"
#include "TSystem.h"
#include "TChain.h"
#include "TObjArray.h"

#include "basicLoop.C"
/*
Usage:
root -b -l -q run_basicLoop.C++

these are not needed (assuming this macro is compiled) because of the include above
(Not a solution that I love, but it works)

.L basicLoop.C++

gSystem->Load("basicLoop_C.so");

*/
const TString version = "V00-01-05/DATA/38X";

void run_ABCDLoop_data()
{

  //the problem with this structure is that it doesn't parallelize
  //parallelization would work well on dellcmscornell

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
 
  //assemble the whole thing into one TChain!
  for (int ifile=0; ifile<nfiles; ifile++) {
    TString samplefiles = dirlist->At(ifile)->GetTitle();
    
    samplefiles+="/*.root";
    
    //    if (!samplefiles.Contains("JetMET-")) continue;
    cout<<"About to add files in: "<<samplefiles<<endl;
    
    ch.Add(samplefiles);
    info.Add(samplefiles);
  }

  basicLoop looper(&ch,&info);  

  //important! this is where cuts are defined
  looper.setCutScheme(basicLoop::kBaseline0);
  
  looper.setMETType(basicLoop::kpfMET);  //other options are basicLoop::ktcMET, basicLoop::kMET, basicLoop::kMHT
  looper.setLeptonType(basicLoop::kPFLeptons); 
  looper.setJetType(basicLoop::kPF); //other options are basicLoop::kCalo
  looper.setBCut(0); //adjust number of b tags here
  
  //careful what is set here!
  looper.setIgnoredCut("cutMET"); //for kRA2
  looper.setIgnoredCut("cutDeltaPhi");
  
  looper.ABCDtree(0);
}


