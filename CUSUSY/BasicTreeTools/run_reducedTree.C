#include "TString.h"
#include "TSystem.h"
#include "TChain.h"
#include "TObjArray.h"

#include "basicLoop.C"
/*
Usage:
root -b -l -q run_reducedTree.C++

output directory is defined
*/
const TString version = "V00-03-01";
const TString ecalDeadCellPath = "/cu1/joshmt/ECALDeadCellNtuples/"; //hard coded for dellcmscornell

const TString outputdir = "/cu2/joshmt/V00-03-01_5/";

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

    //    if (samplefiles.Contains("QCD") && !samplefiles.Contains("PU")) continue;
//     if ( samplefiles.Contains("WJetsZ2")) continue;

//     if ( samplefiles.Contains("LM")) continue;

//     if ( samplefiles.Contains("QCD-Pt0")) continue;
//     if ( samplefiles.Contains("QCD-Pt1")) continue;
//     if ( samplefiles.Contains("QCD-Pt2")) continue;
//     if ( samplefiles.Contains("QCD-Pt3")) continue;
//     if ( samplefiles.Contains("QCD-Pt4")) continue;


    //fetch ecal dead cell info
    TFile fEcal(ecalDeadCellPath + samplename + "/deadCellFilterProfile.root"); //hadd files together into one
    TTree* ecalTree = fEcal.IsZombie() ? 0 : (TTree*) fEcal.Get("filter");

    TChain ch("BasicTreeMaker/tree");
    TChain info("BasicTreeMaker/infotree");
    ch.Add(samplefiles);
    info.Add(samplefiles);
    basicLoop looper(&ch,&info,ecalTree);
    
    //important! this is where cuts are defined
    looper.setCutScheme(basicLoop::kBaseline0);
    looper.setMETType(basicLoop::kpfMET);
    looper.setMETRange(basicLoop::kHigh); //signal region
    looper.setLeptonType(basicLoop::kPFLeptonsRA2);
    looper.setJetType(basicLoop::kPF);
    looper.setCleaningType(basicLoop::kMuonCleaning);
    looper.setJERType(basicLoop::kJERbias6);
    looper.setDPType(basicLoop::kminDP);

    looper.reducedTree(outputdir);  //go!
  }

}
