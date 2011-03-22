#include "TString.h"
#include "TSystem.h"
#include "TChain.h"
#include "TObjArray.h"

#include "basicLoop.C"
/*
Usage:
root -b -l -q run_reducedTree.C++

*/
const TString version = "V00-03-01";
const TString ecalDeadCellPath = "/cu1/joshmt/ECALDeadCellNtuples/"; //hard coded for dellcmscornell

void run_reducedTree()
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
    //before tacking on the *.root, let's grab the full sample name
    TString samplename = samplefiles( samplefiles.Last('/')+1, samplefiles.Length() );
    samplefiles+="/*.root";

    if(samplefiles.Contains("DATA")) continue;
    cout<<"About to start on sample: "<<samplename<<endl;
    if (!samplefiles.Contains("WJets") ) continue;
    //if (samplefiles.Contains("QCD") && !samplefiles.Contains("PU")) continue;

  //want to skip samples that aren't ready, and also ZJets
//     if (!
//  	(samplefiles.Contains("LM")
// 	 // 	 || samplefiles.Contains("PU")
// 	 // 	 || samplefiles.Contains("Single")
//  	 || samplefiles.Contains("TTbarJets")
//  	 )
//  	) continue;
    
    //fetch ecal dead cell info
    TFile fEcal(ecalDeadCellPath + samplename + "/deadCellFilterProfile.root"); //hadd files together into one
    TTree* ecalTree = fEcal.IsZombie() ? 0 : (TTree*) fEcal.Get("filter");

    TChain ch("BasicTreeMaker/tree");
    TChain info("BasicTreeMaker/infotree");
    ch.Add(samplefiles);
    info.Add(samplefiles);
    basicLoop looper(&ch,&info,ecalTree);
    
    //important! this is where cuts are defined
    looper.setCutScheme(basicLoop::kBaseline0); //usually this is kRA2
    looper.setMETType(basicLoop::kpfMET);
    looper.setMETRange(basicLoop::kHigh); //signal region
    //looper.setLeptonType(basicLoop::kPFLeptons);
    looper.setLeptonType(basicLoop::kPFLeptonsRA2);
    looper.setJetType(basicLoop::kPF);
    //looper.setCleaningType(basicLoop::kMuonCleaning);
    looper.setCleaningType(basicLoop::kMuonEcalCleaning);
    looper.setJERType(basicLoop::kJERbias);
    looper.setDPType(basicLoop::kminDP);

    looper.reducedTree("/home/joshmt/");  //go!
  }

}
