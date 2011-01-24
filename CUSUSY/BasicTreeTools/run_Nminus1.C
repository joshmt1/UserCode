#include "TString.h"
#include "TSystem.h"
#include "TChain.h"
#include "TObjArray.h"

#include "basicLoop.C"
/*
by adding the include of the .C file, we just compile all of basicLoop into this so file.
This is an easy solution and works ok for these simple classes.

just do:
root -b -l -q run_Nminus1.C++
*/

const TString version = "V00-01-05";

//const TString extrapath = "SUSYPATv8_363"; //pass an empty string unless you need something else
const TString extrapath = "";

void run_Nminus1()
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
    samplefiles+="/*.root";

    cout<<"About to start on files: "<<samplefiles<<endl;

    if (samplefiles.Contains("DATA")) continue; //skip data

    //if (!(samplefiles.Contains("QCD") )) continue; //hack to skip some samples
    
    TChain ch("BasicTreeMaker/tree");
    TChain info("BasicTreeMaker/infotree");
    ch.Add(samplefiles);
    info.Add(samplefiles);
    basicLoop looper(&ch,&info);

    looper.setCutScheme(basicLoop::kBaseline0);
    looper.setMETType(basicLoop::kpfMET);
    looper.setMETRange(basicLoop::kHigh); //signal region
    //    looper.setMETRange(basicLoop::kMedhigh); //sideband region!
    looper.setJetType(basicLoop::kPF);
    looper.setLeptonType(basicLoop::kPFLeptons);
    looper.setDPType(basicLoop::kminDP);

    looper.Nminus1plots();
  }


}
