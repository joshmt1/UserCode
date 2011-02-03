#include "TString.h"
#include "TSystem.h"
#include "TChain.h"
#include "TObjArray.h"

#include <fstream>

#include "basicLoop.C"
/*
use:
root -b -l -q run_cutflow.C++
*/

const TString version = "V00-02-00";

//const TString extrapath = "SUSYPATv8_363"; //pass an empty string unless you need something else
const TString extrapath = "";

void run_cutflow_castor()
{

  int pid = gSystem->GetPid();

  TString    dir="/castor/cern.ch/user/k/kreis/analysis/SUSY/DonPAT/";
  if (extrapath!="") { dir += extrapath; dir += "/";}
  dir += version; dir+="/";

  //for castor we've got to do all sorts of crap
  TString outname_dir = ".castor_tmp_"; outname_dir+=pid; outname_dir+="_dir";
  TString command; command.Form("nsls %s > %s",dir.Data(),outname_dir.Data());
  gSystem->Exec(command); //run nsls
  ifstream dirlistfile(outname_dir.Data());
  char dirname[200];
  while (dirlistfile>>dirname) {
    //so this is the standard loop over samples

    TString samplename(dirname);
    if (samplename.Contains("DATA")) continue; //skip data (use run_cutflow_data.C)

    //put here the standard hacks to skip samples, if desired
    //here we are skipping all samples except ttbar
    //    if (!samplename.Contains("TTbarJets")) continue;
 

    TString samplepath = dir;
    samplepath+= samplename;

    TString outname=".castor_tmp_";  outname+=pid; outname+=samplename;
    command.Form("nsls %s > %s",samplepath.Data(),outname.Data());
    gSystem->Exec(command);

    samplepath+="/";

    //this is then a loop over files in the sample
    //so let's prepare our chains
    TChain ch("BasicTreeMaker/tree");
    TChain info("BasicTreeMaker/infotree");

    ifstream filelistfile(outname.Data());
    char filename[200];
    while (filelistfile>>filename) { //here's the loop

      if (!TString(filename).BeginsWith("BasicN")) continue; //do only BasicN*

      TString fullcastorpath=samplepath;
      fullcastorpath+=filename;
      //      command.Form("stager_qry -M %s > ._castor_tmp_%d_qry",fullcastorpath.Data()); gSystem->Exec(command.Data());

      fullcastorpath.Prepend("rfio_lseek64://");

      ch.Add(fullcastorpath);
      info.Add(fullcastorpath);
    }
    //once we get to here, we've added the full list of files to the Chain
    //at this point we go to the "normal" code

    cout<<"Running on "<<ch.GetListOfFiles()->GetEntries()<<" files from: "<<samplepath<<endl;

    basicLoop looper(&ch,&info);
    //    looper.setSpecialCutDescription("noJetID");

    looper.setCutScheme(basicLoop::kBaseline0);
    looper.setMETType(basicLoop::kpfMET);
    looper.setMETRange(basicLoop::kHigh); //signal region
    //looper.setMETRange(basicLoop::kMedium); //50-100 GeV region
    looper.setJetType(basicLoop::kPF);
    looper.setLeptonType(basicLoop::kPFLeptons);
    looper.setDPType(basicLoop::kminDP);

    looper.setCleaningType(basicLoop::kMuonCleaning);

    looper.setBCut(3); //require 3 b tags so that we make the full cut flow table

    looper.setMuonReq(1); //inverted muon veto

    looper.cutflow(false); //true means write verbose event count files
  }

  command.Form("rm .castor_tmp_%d*",pid);
  gSystem->Exec(command);

}
