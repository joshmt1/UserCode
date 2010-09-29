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

const TString version = "V00-00-04";

void run_cutflow()
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

  for (int ifile=0; ifile<nfiles; ifile++) {
    TString samplefiles = dirlist->At(ifile)->GetTitle();
    //for the first iteration I had this stupid subdirectory in there
    if (version=="V00-00-01")  samplefiles+="/Spring10/*.root";
    else                       samplefiles+="/*.root";

    cout<<"About to start on files: "<<samplefiles<<endl;

    //if (!samplefiles.Contains("v2")) continue; //hack to skip some samples

    TChain ch("BasicTreeMaker/tree");
    ch.Add(samplefiles);
    basicLoop looper(&ch);
    looper.setCutScheme(basicLoop::kRA2tcMET); //this is a critical line! defines the cuts to use!
    looper.setBCut(3); //require 3 b tags so that we make the full cut flow table
    looper.cutflow();
  }

  //the old fashioned way

//   {
//     //QCD-Pt100to250-madgraph
//     TChain* qcd100Chain = new TChain("BasicTreeMaker/tree");
//     qcd100Chain->Add("/cu1/joshmt/BasicNtuples/V00-00-01/QCD-Pt100to250-madgraph/Spring10/BasicNtuple*");
    
//     basicLoop l_qcd100(qcd100Chain);
//     l_qcd100.cutflow();
//   }

//   {
//   //QCD-Pt250to500-madgraph
//   TChain* qcd250Chain = new TChain("BasicTreeMaker/tree");
//   qcd250Chain->Add("/cu1/joshmt/BasicNtuples/V00-00-01/QCD-Pt250to500-madgraph/Spring10/BasicNtuple*");

//   basicLoop l_qcd250(qcd250Chain);
//   l_qcd250.cutflow();
//   }
  

//   //QCD-Pt500to1000-madgraph
//   TChain* qcd500Chain = new TChain("BasicTreeMaker/tree");
//   qcd500Chain->Add("/cu1/joshmt/BasicNtuples/V00-00-01/QCD-Pt500to1000-madgraph/Spring10/BasicNtuple*");

//   basicLoop l_qcd500(qcd500Chain);
//   l_qcd500.cutflow();
  
//   //QCD-Pt1000toInf-madgraph
//   TChain* qcd1000Chain = new TChain("BasicTreeMaker/tree");
//   qcd1000Chain->Add("/cu1/joshmt/BasicNtuples/V00-00-01/QCD-Pt1000toInf-madgraph/Spring10/BasicNtuple*");

//   basicLoop l_qcd1000(qcd1000Chain);
//   l_qcd1000.cutflow();

}
