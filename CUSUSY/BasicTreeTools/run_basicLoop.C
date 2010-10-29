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
const TString version = "V00-01-01";

void run_basicLoop()
{

  //the problem with this structure is that it doesn't parallelize
  //parallelization would work well on dellcmscornell

  TString computername = gSystem->Getenv("HOST");
  TString dir = "/cu1/joshmt/"; //default is dellcmscornell
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

  for (int ifile=0; ifile<nfiles; ifile++) {
    TString samplefiles = dirlist->At(ifile)->GetTitle();
    //for the first iteration I had this stupid subdirectory in there
    if (version=="V00-00-01")  samplefiles+="/Spring10/*.root";
    else                       samplefiles+="/*.root";


    if (samplefiles.Contains("DATA")) continue; //skip data (use run_basicLoop_data.C)

    //if (!samplefiles.Contains("QCD")) continue; //hack to skip some samples
    if (!(samplefiles.Contains("LM") || samplefiles.Contains("TTbar"))) continue; //hack to skip some samples

    cout<<"About to start on files: "<<samplefiles<<endl;

    TChain ch("BasicTreeMaker/tree");
    TChain info("BasicTreeMaker/infotree");
    ch.Add(samplefiles);
    info.Add(samplefiles);
    basicLoop looper(&ch,&info);
    //important! this is where cuts are defined
    looper.setCutScheme(basicLoop::kRA2); //this is now the only scheme!
    looper.setMETType(basicLoop::kMET);
    looper.setMETRange(basicLoop::kWide);
    looper.setDPType(basicLoop::kminDP);
    //no b tagging cut so that plots can apply it selectively

    //    looper.setCutScheme(basicLoop::kRA2minDP);

    //careful what is set here!
    //looper.setIgnoredCut(basicLoop::cutTrigger);
    //    looper.setIgnoredCut(basicLoop::cutHT);
    //looper.setIgnoredCut(basicLoop::cutMET); //MET
    //looper.setIgnoredCut(basicLoop::cutDeltaPhi); //DeltaPhi
    
    looper.Loop();  //go!
  }

}
