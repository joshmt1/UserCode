#include "TString.h"
#include "TSystem.h"
#include "TChain.h"
#include "TObjArray.h"

#include "basicLoop.C"
/*
use:
root -b -l -q run_cutflow.C++
*/

const TString version = "V00-03-01";

//const TString extrapath = "SUSYPATv8_363"; //pass an empty string unless you need something else
const TString extrapath = "";
const TString ecalDeadCellPath = "/cu1/joshmt/ECALDeadCellNtuples/"; //hard coded for dellcmscornell

void run_cutflow()
{

  //TString libname="basicLoop_C.so";
  TString computername = gSystem->Getenv("HOST");
  TString dir = "/cu1/joshmt/";
  if (computername =="JoshPC") {
    dir="~/data/";
  }
  else if (computername.Contains("lxplus")) {
    dir="/tmp/joshmt/";
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
    //before tacking on the *.root, let's grab the full sample name
    TString samplename = samplefiles( samplefiles.Last('/')+1, samplefiles.Length() );
    samplefiles+="/*.root";

    cout<<"About to start on sample: "<<samplename<<endl;

    if (samplefiles.Contains("DATA")) continue; //skip data (use run_cutflow_data.C)
    
    if (!samplefiles.Contains("LM9") ) continue;
    //if (!samplefiles.Contains("LM13") ) continue;
    //if (!samplefiles.Contains("TTbarJets") ) continue;
    //if (!samplefiles.Contains("WJets") ) continue;
    //if (!samplefiles.Contains("ZJets") ) continue;
    //if (!samplefiles.Contains("SingleTop-sChannel") ) continue;
    //if (!samplefiles.Contains("SingleTop-tWChannel") ) continue;
    //if (!samplefiles.Contains("Zinvisible") ) continue;

    //if (samplefiles.Contains("QCD") && !samplefiles.Contains("PU")) continue;
    //if (!(samplefiles.Contains("QCD") && samplefiles.Contains("PU"))) continue;

    TChain ch("BasicTreeMaker/tree");
    TChain info("BasicTreeMaker/infotree");
    ch.Add(samplefiles);
    info.Add(samplefiles);
    //    cout<<    ch.GetEntries()<<endl; continue; //for quickly checking the number of entries in the samples

    //fetch ecal dead cell info
    TFile fEcal(ecalDeadCellPath + samplename + "/deadCellFilterProfile.root"); //hadd files together into one
    TTree* ecalTree = fEcal.IsZombie() ? 0 : (TTree*) fEcal.Get("filter");

    basicLoop looper(&ch,&info,ecalTree);

    //    looper.setSpecialCutDescription("noJetID");

    looper.setCutScheme(basicLoop::kBaseline0);
    looper.setMETType(basicLoop::kpfMET);
    looper.setMETRange(basicLoop::kHigh); //signal region
    //looper.setMETRange(basicLoop::kMedhigh);
    looper.setJetType(basicLoop::kPF);
    looper.setLeptonType(basicLoop::kPFLeptonsRA2);
    looper.setDPType(basicLoop::kminDP);

    looper.setCleaningType(basicLoop::kMuonCleaning);
    looper.setJERType(basicLoop::kJERbias);

    //looper.setJESType(basicLoop::kJESup);  
    //    looper.setMETuncType(basicLoop::kMETuncUp);

    looper.setBCut(3); //require 3 b tags so that we make the full cut flow table

    //    looper.setMuonReq(1); //inverted muon veto
    //looper.setRequiredCut("cut2SUSYb");

    looper.cutflow(false); //true means write verbose event count files

    ////does this work?
    looper.calculateTagProb();//default
    //looper.setBTagEffType(basicLoop::kBTagEffup);
    //looper.calculateTagProb();
    //looper.setBTagEffType(basicLoop::kBTagEffdown);
    //looper.calculateTagProb();

  }

}
