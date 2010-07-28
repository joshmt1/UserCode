#include "TString.h"
#include "TSystem.h"
#include "TChain.h"
#include "TObjArray.h"

#include "basicLoop.C"
/*

.L basicLoop.C++

gSystem->Load("basicLoop_C.so");

*/
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

  dir += "BasicNtuples/V00-00-01/";
  TChain dummy("dummy");
  TString dirs = dir; dirs+="*";
  dummy.Add(dirs);
  TObjArray* dirlist = dummy.GetListOfFiles();
  int nfiles=dirlist->GetEntries();

  for (int ifile=0; ifile<nfiles; ifile++) {
    TString samplefiles = dirlist->At(ifile)->GetTitle();
    samplefiles+="/Spring10/*.root";

    cout<<"About to start on files: "<<samplefiles<<endl;

    if (samplefiles.Contains("LM")) continue; //hack to skip some samples

    TChain ch("BasicTreeMaker/tree");
    ch.Add(samplefiles);
    cout<<"files added..."<<endl;
    basicLoop looper(&ch);
    cout<<"class constructed"<<endl;
    looper.setCutScheme(basicLoop::kRA2withB);
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
