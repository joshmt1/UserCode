#include "TString.h"
#include "TSystem.h"
#include "TChain.h"
#include "TObjArray.h"

#include "basicLoop.C"
/*
Usage:
root -b -l -q run_ABCDLoop.C++

Updated for ECAL dead cell cleaning!
*/
const TString version = "V00-03-01";
const TString ecalDeadCellPath = "/cu1/joshmt/ECALDeadCellNtuples/"; //hard coded for dellcmscornell

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
    //before tacking on the *.root, let's grab the full sample name
    TString samplename = samplefiles( samplefiles.Last('/')+1, samplefiles.Length() );
    samplefiles+="/*.root";

    if(samplefiles.Contains("DATA")) continue;
    if(!samplefiles.Contains("PythiaZ2-PU2010")) continue;
    cout<<"About to start on sample: "<<samplename<<endl;

    //fetch ecal dead cell info
    TFile fEcal(ecalDeadCellPath + samplename + "/deadCellFilterProfile.root"); //hadd files together into one
    TTree* ecalTree = fEcal.IsZombie() ? 0 : (TTree*) fEcal.Get("filter");

    TChain ch("BasicTreeMaker/tree");
    TChain info("BasicTreeMaker/infotree");
    ch.Add(samplefiles);
    info.Add(samplefiles);
    basicLoop looper(&ch,&info,ecalTree);
    
    //important! this is where cuts are defined
    //looper.setCutScheme(basicLoop::kRA2); //usually this is kRA2
    looper.setCutScheme(basicLoop::kBaseline0); //usually this is kRA2
    looper.setMETType(basicLoop::kpfMET);
    looper.setDPType(basicLoop::kminDP);
    looper.setLeptonType(basicLoop::kPFLeptonsRA2);
    looper.setJetType(basicLoop::kPF);
    looper.setCleaningType(basicLoop::kMuonCleaning);
    looper.setJERType(basicLoop::kJERbias);

    looper.setBCut(0);
    
    //careful what is set here!
    looper.setIgnoredCut("cutMET"); 
    looper.setIgnoredCut("cutDeltaPhi");
    
    looper.ABCDtree();  //go!
  }

}
