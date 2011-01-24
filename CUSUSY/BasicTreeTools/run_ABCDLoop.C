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
const TString version = "V00-01-05";

void run_ABCDLoop()
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

  for (int ifile=0; ifile<nfiles; ifile++) {
    TString samplefiles = dirlist->At(ifile)->GetTitle();
    //for the first iteration I had this stupid subdirectory in there
    if (version=="V00-00-01")  samplefiles+="/Spring10/*.root";
    else                       samplefiles+="/*.root";

    cout<<"About to start on files: "<<samplefiles<<endl;

    if(!(samplefiles.Contains("LM13")
	 || samplefiles.Contains("TTbar")
	 || samplefiles.Contains("Zinv")
	 || samplefiles.Contains("WJets")
	 || samplefiles.Contains("ZJets")
	 || samplefiles.Contains("SingleTop")
	 || samplefiles.Contains("QCD")) ) continue; //hack to skip some samples   

    TChain ch("BasicTreeMaker/tree");
    TChain info("BasicTreeMaker/infotree");
    ch.Add(samplefiles);
    info.Add(samplefiles);
    basicLoop looper(&ch,&info);
    
    //important! this is where cuts are defined
    //looper.setCutScheme(basicLoop::kRA2); //usually this is kRA2
    looper.setCutScheme(basicLoop::kBaseline0); //usually this is kRA2

    looper.setMETType(basicLoop::kpfMET);
    looper.setLeptonType(basicLoop::kPFLeptons);
    looper.setJetType(basicLoop::kPF);
    looper.setBCut(0);
    
    //careful what is set here!
    looper.setIgnoredCut("cutMET"); 
    looper.setIgnoredCut("cutDeltaPhi");
    
    looper.ABCDtree();  //go!
  }

}
