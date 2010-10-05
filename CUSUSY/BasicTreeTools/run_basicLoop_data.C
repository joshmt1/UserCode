#include "TString.h"
#include "TSystem.h"
#include "TChain.h"
#include "TObjArray.h"

#include "basicLoop.C"
/*
Usage:
root -b -l -q run_basicLoop_data.C++
*/
const TString version = "V00-00-05/DATA/38X";

void run_basicLoop_data()
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

  unsigned int currentfileindex=1;
  for (int ifile=0; ifile<nfiles; ifile++) {
    TString samplefiles = dirlist->At(ifile)->GetTitle();
    samplefiles+="/*.root";

    cout<<"About to start on files: "<<samplefiles<<endl;

    //    if (!samplefiles.Contains("LM13")) continue; //hack to skip some samples

    TChain ch("BasicTreeMaker/tree");
    ch.Add(samplefiles);
    basicLoop looper(&ch);
    //important! this is where cuts are defined
    looper.setCutScheme(basicLoop::kRA2METminDP);
    //    looper.setCutScheme(basicLoop::kRA2minDP);
    //no b tagging cut so that plots can apply it selectively
    //careful what is set here!
    //looper.setIgnoredCut(basicLoop::cutTrigger);
    //    looper.setIgnoredCut(basicLoop::cutHT);
    //looper.setIgnoredCut(basicLoop::cutMET); //MET
    //    looper.setIgnoredCut(basicLoop::cutDeltaPhi); //DeltaPhi
    
    looper.Loop(currentfileindex);  //go!
    currentfileindex++;
  }

}
