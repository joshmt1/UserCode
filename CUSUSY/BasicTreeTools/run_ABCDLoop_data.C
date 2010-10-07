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
const TString version = "V00-00-05/DATA/38X";

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

  unsigned int currentfileindex=1;
  for (int ifile=0; ifile<nfiles; ifile++) {
    TString samplefiles = dirlist->At(ifile)->GetTitle();
    samplefiles+="/*.root";

    cout<<"About to start on files: "<<samplefiles<<endl;

    //if (!samplefiles.Contains("LM13")) continue; //for V00-00-04b
    
    //for V00-00-04

    TChain ch("BasicTreeMaker/tree");
    ch.Add(samplefiles);
    basicLoop looper(&ch);
    //important! this is where cuts are defined
    looper.setCutScheme(basicLoop::kRA2METminDP);
    //looper.setCutScheme(basicLoop::kRA2METMPT);
    looper.setBCut(0); //adjust number of b tags here
    //careful what is set here!
    looper.setIgnoredCut(basicLoop::cutMET); //MET
    //looper.setIgnoredCut(basicLoop::cutMHT); //for kRA2
    looper.setIgnoredCut(basicLoop::cutDeltaPhi); //DeltaPhi

    looper.ABCDtree(currentfileindex++);  //go!
  }

}
