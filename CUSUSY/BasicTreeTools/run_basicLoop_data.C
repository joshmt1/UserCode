#include "TString.h"
#include "TSystem.h"
#include "TChain.h"
#include "TObjArray.h"

#include "basicLoop.C"
/*
Usage:
root -b -l -q run_basicLoop_data.C++
*/
const TString version = "V00-01-04/DATA/386";

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
    TChain info("BasicTreeMaker/infotree");
    ch.Add(samplefiles);
    info.Add(samplefiles);
    basicLoop looper(&ch,&info);
    //important! this is where cuts are defined
    looper.setCutScheme(basicLoop::kBaseline0);
    looper.setMETType(basicLoop::kpfMET);
    looper.setMETRange(basicLoop::kWide); //get more events for plotting
    looper.setJetType(basicLoop::kPF);
    looper.setLeptonType(basicLoop::kPFLeptons);
    looper.setDPType(basicLoop::kminDP); //apply the CU minDeltaPhi cut

    looper.setIgnoredCut("cut3Jets"); //N-1

    looper.setBCut(0);  //no b tagging cut so that plots can apply it selectively

    looper.Loop(currentfileindex);  //go!
    currentfileindex++;
  }

}
