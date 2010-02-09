#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"

/*
gSystem->Load("eventLoop_C.so");
init()
selectData(sample);
setCutsStandard()
setCutsDon()
setCutsRA2()
eventLoop()

*/

#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TAxis.h"
#include "TMatrixT.h"
#include "TVectorT.h"
#include <iostream>
#include <vector>

#include <exception>
#include <cmath> 
#include <iomanip>
#include "NtupleMaker/CUJEMtuple/interface/CUselect.h"

#include "/afs/cern.ch/user/j/joshmt/root/util/HistHolder.cxx"

//need to ask darren about this
//#if !defined(__CINT__) && !defined(__MAKECINT__)

//Headers for the data items
#include "ProdTutorial/CUcollections/interface/CUelectron.h"
#include "ProdTutorial/CUcollections/interface/CUevent.h"
#include "ProdTutorial/CUcollections/interface/CUjet.h"
#include "ProdTutorial/CUcollections/interface/CUmcparticle.h"
#include "ProdTutorial/CUcollections/interface/CUmet.h"
#include "ProdTutorial/CUcollections/interface/CUmuon.h"
#include "ProdTutorial/CUcollections/interface/CUphoton.h"
#include "ProdTutorial/CUcollections/interface/CUsupercluster.h"
#include "ProdTutorial/CUcollections/interface/CUtrack.h"
#include "ProdTutorial/CUcollections/interface/CUprimaryvertex.h"
#include "ProdTutorial/CUcollections/interface/CUtrigger.h"
#include "ProdTutorial/CUcollections/interface/CUskimbits.h"

//#endif

//******************************************************************************
typedef std::vector<std::vector<double> > vvector;
  typedef CUelectronCollection::const_iterator      EleIter;
  typedef CUjetCollection::const_iterator           JetIter;
  typedef CUmcparticleCollection::const_iterator    MCparIter;
  typedef CUmuonCollection::const_iterator          MuonIter;
  typedef CUphotonCollection::const_iterator        PhotonIter;
  typedef CUsuperclusterCollection::const_iterator  SCIter;
  typedef CUprimaryvertexCollection::const_iterator  VtxIter;
  typedef CUtrackCollection::const_iterator         TrackIter;
  typedef CUtriggerCollection::const_iterator       TrigIter;

//******************************************************************************

  // begin custom code
using namespace std;
//    void setCutsDon();

//    bool passJetCuts(int index);
//    bool passPreselection( Int_t i1, Int_t i2);

TFile* file;
vector<string> fileNames;
//variables to "turn off" parts of the event loop
bool useTriggerInfo_;
bool useScInfo_;
bool useMcInfo_;
bool useTrackInfo_;

   //variables for cuts
   TString cuts_;
   float minJetPT_;
   float  minLooseJetPT_ ;
   float minEMfrac_;
   float maxEMfrac_;
   float maxJetEta_;
   unsigned int minGoodJets_;
   bool badJetVeto_;
   bool leptonVeto_;
   bool photonVeto_;
   bool requirePrimaryVertex_;
   float minJetPT1_preSel_;
   float minJetPT2_preSel_;
   float minJetPT3_preSel_;

   float minJetPT1_Sel_;
   float minJetPT2_Sel_;
   float minJetPT3_Sel_;

   float maxJetEta_preSel_;
   float minMET_;
   float minDeltaPhi_; //between lead jets and MET
   float minSSV_;
   float minbjets_;
   float minBJetPT_;
   //other stuff
   TString sampleid_;
   // end custom code

//This auto-file naming should be revived, but not yet
void init() {

   //begin custom code
   //this depends on the format of the file path being .../treeSAMPLENAME_INDEX.root
//    if (tree->GetCurrentFile() == 0) { ///for TChain
//      TString firstfile = ((TChain*) tree)->GetListOfFiles()->At(0)->GetTitle();
//      int slashpos = firstfile.Index("/tree");
//      sampleid_=  firstfile(slashpos+5,firstfile.Last('_')-slashpos-5);
//    }
//    else {  //for TTree
//      TString filename= tree->GetCurrentFile()->GetName();
//      sampleid_=  filename(4,filename.Last('_')-4);
//    }
   //end custom code
  sampleid_ = "LM9";
  file = 0;

  cuts_="";
  useTriggerInfo_=false;
  useScInfo_=false;
  useMcInfo_=false;
  useTrackInfo_=false;
  requirePrimaryVertex_=false;

}

void setCutsStandard() {
   // begin custom code
   //object definition cuts
  cuts_="JMT";

  minJetPT_ = 50;
  minLooseJetPT_ = 30;
  minEMfrac_ = 0.05;
  maxEMfrac_ = 0.95;
  maxJetEta_ = 2.4;
  minGoodJets_ = 3;
  badJetVeto_ = true;
  leptonVeto_ = true;
  photonVeto_ = false; //should probably turn this on, after i implement some better selection for photons

  useTrackInfo_=true;
  requirePrimaryVertex_=true;

  //preselection cuts
  minJetPT1_preSel_ = 100; //the value of this cut depends on what analysis you look at
  minJetPT2_preSel_ = minJetPT_;
  maxJetEta_preSel_ = maxJetEta_;
  //MET cut
  minMET_ = 100; //following Don
  minDeltaPhi_ = 0.2; //following Don

  minSSV_ = 1.74; //7TeV
  minbjets_ = 2;
  minBJetPT_ = 50;
  // end custom code
  sampleid_+= "_CutsJMT";

  
}

void setCutsDon() {
  std::cout<<"Setting cuts for 'Don' mode!"<<std::endl;

  minJetPT_=50;
  minEMfrac_=-1; //no cut
  maxEMfrac_=2; //no cut
  maxJetEta_=2.4;
  minGoodJets_=4;
  badJetVeto_=false;
  leptonVeto_=false;
  photonVeto_=false;
    //preselection cuts
  minJetPT1_preSel_=100;
  minJetPT2_preSel_=minJetPT_; //no cut on second jet
  maxJetEta_preSel_ = maxJetEta_;

  minMET_=100;

  minDeltaPhi_ = 0.2;
  minSSV_=1.74;
  minbjets_=3;
  minBJetPT_ = 30;
  sampleid_+= "_CutsDon";
}

void setCutsRA2() {
  std::cout<<"Setting cuts for 'RA2' mode!"<<std::endl;
  cuts_="RA2";

  minJetPT_=50;
  minEMfrac_=0.05;
  maxEMfrac_=0.95;
  maxJetEta_=2.5;

  //preselection cuts
  minJetPT1_preSel_=50;
  minJetPT2_preSel_=50;
  minJetPT3_preSel_=30;
  maxJetEta_preSel_ = 5;

  minGoodJets_=3;

  leptonVeto_=true;
  photonVeto_=false;

  //preselection cuts

  minJetPT1_Sel_=180;
  minJetPT2_Sel_=150;
  minJetPT3_Sel_=50;

  badJetVeto_=false;

  minMET_=200;

  minDeltaPhi_ = 0.3;
  minSSV_=1.74;
  minbjets_=3;
  minBJetPT_ = 50;
  sampleid_+= "_CutsRA2";

}


//bool passJetCuts() {return true;}
//bool passPreselection() {return true;}


bool passJetCuts(JetIter jet) {

  //could add a sanity check of i against njets

  //pT
  if ( jet->pt < minJetPT_ ) return false;

  //eta
  if ( fabs( jet->eta ) > maxJetEta_ ) return false;

  //EM fraction
  if ( jet->EMfrac > maxEMfrac_) return false;
  if ( jet->EMfrac < minEMfrac_) return false;

  return true;
}


bool passLooseJetCuts(JetIter jet) {

  //is this evil?
  float minJetPT = minJetPT_;
  minJetPT_ = minLooseJetPT_;
  bool pass = passJetCuts(jet);
  minJetPT_=minJetPT;
  return pass;

}

bool passBJetCuts(JetIter jet) {

  //largely copied from passJetCuts()

  //pT
  if ( jet->pt < minBJetPT_ ) return false;

  //b tagging with SSV
  if ( jet->btagSecVertex < minSSV_) return false;

  //for now eta and EM are the same as for generic jets
  //eta
  if ( fabs( jet->eta ) > maxJetEta_ ) return false;

  //EM fraction
  if ( jet->EMfrac > maxEMfrac_) return false;
  if ( jet->EMfrac < minEMfrac_) return false;

  return true;
}

//this is the SUS-09-001 preselection
bool passPreselection(JetIter j1, JetIter j2) {
  if ( j1->pt < minJetPT1_preSel_ ) return false;

  if ( j2->pt < minJetPT2_preSel_ ) return false;

  if ( fabs(j1->eta) > maxJetEta_preSel_ ) return false;

  return true;
}

bool passRA2JetPreselection(JetIter j1, JetIter j2, JetIter j3) {
  if ( j1->pt < minJetPT1_preSel_ ) return false;

  if ( j2->pt < minJetPT2_preSel_ ) return false;

  if ( j3->pt < minJetPT3_preSel_ ) return false;

  if ( fabs(j1->eta) > maxJetEta_preSel_ ) return false;
  if ( fabs(j2->eta) > maxJetEta_preSel_ ) return false;
  if ( fabs(j3->eta) > maxJetEta_preSel_ ) return false;

  return true;
}

bool passRA2JetCuts(JetIter j1, JetIter j2, JetIter j3) {
  if ( j1->pt < minJetPT1_Sel_ ) return false;
  //cut on eta and EMfrac
  if ( !passJetCuts(j1) ) return false;

  if ( j2->pt < minJetPT2_Sel_ ) return false;
  if ( !passJetCuts(j2) ) return false;

  if ( j3->pt < minJetPT3_Sel_ ) return false;
  if ( !passJetCuts(j3) ) return false;

  return true;
}

double getDeltaPhi(TrackIter tr, CUmet met) {

  //naive me..i wanted eta for met! that doesn't work!

  return acos(cos( met.phi - tr->phi));
  //  double deltaeta = met.eta - tr->eta;

  //  return sqrt( deltaphi*deltaphi + deltaeta*deltaeta);
}

double getDeltaPhi(JetIter j1, CUmet cumet) {
  double metphi = cumet.phi;
  double j1_phi = j1->phi;
  return acos(cos(metphi-j1_phi));
}

double getDeltaPhi(JetIter j1, JetIter j2) {
  return acos(cos(j2->phi - j1->phi));
}

double getDeltaPhiMC(JetIter j1, JetIter j2) {

//   cout<<"--getDeltaPhiMC--"<<endl;
//   cout<<j1->pt<<endl;
//   cout<<j1->genPartonPhi<<endl;
//   cout<<j2->pt<<endl;
//   cout<<j2->genPartonPhi<<endl;
//   cout<<"----"<<endl;
  return acos(cos(j2->genPartonPhi - j1->genPartonPhi));
}

bool passDeltaPhi(JetIter j1, CUmet cumet) {

  if ( getDeltaPhi(j1,cumet) < minDeltaPhi_ ) return false;

  return true;
}

//end custom code

void selectData(TString sample, TString maxindex="") {

  TString path="rfio:/castor/cern.ch/user/p/puigh/CUSusy/CUJEM/Summer09/7TeV/Output/";
  sampleid_ = sample;
  
  if (sample.Contains("QCD")) path+="QCD-madgraph/";

  if (sample.Contains("LM")) {
    path+="LM/";
    path+=sampleid_;
    path += "_Summer09_7TeV_CUJEM_V09.root";
    cout<<"Adding to list of input files: "<<path<<endl;
    fileNames.push_back(string(path.Data()));
  }
  else  {
    if (maxindex!="") {
      for (int ind=1; ind<=maxindex.Atoi(); ind++) {
	TString mypath=path;
	mypath+=sampleid_;
	mypath += "_Summer09_7TeV_CUJEM_V09_";
	mypath +=ind;
	mypath +=".root";
	cout<<"Adding to list of input files: "<<mypath<<endl;
	fileNames.push_back(string(mypath.Data()));
      }
    }
    else {
      path+=sampleid_;
      path += "_Summer09_7TeV_CUJEM_V09.root";
      cout<<"Adding to list of input files: "<<path<<endl;
      fileNames.push_back(string(path.Data()));
    }
  }

}

void eventLoop() {

  // Load the root file to be analyzed
  //  fwlite::Event ev(&*file);

//   TFile  * file = new TFile("myOutputFile.root");
//   fwlite::Event ev(&*file);


  //vector<string> fileNames; //holds all the file names to be linked
  //fileNames.push_back("/afs/cern.ch/user/p/puigh/public/CUSusy/SampleEDMntuple_V0.root");
  //fileNames.push_back("/castor/cern.ch/user/p/puigh/CUSusy/CUJEM/Summer09/7TeV/Output/");
  //fileNames.push_back("/afs/cern.ch/user/p/puigh/public/CUSusy/SUSYPAT_Data.root");
  
  //creates a ChainEvent allowing files to be linked   
  fwlite::ChainEvent ev(fileNames);   


  HistHolder histo; //structure to hold histograms

  TString outfile = "plots_";
  outfile+=sampleid_;
  outfile+=".root";
  TFile fout(outfile,"RECREATE");
  int nbins=200;
  float min=0;
  //  float max=250;
  float htmin=0;
  float htmax=2250;
  histo.make("H_HT","HT",nbins,htmin,htmax);

  float maxmet=800;
  histo.make2("H_MHTMET","MHT versus MET",nbins,min,maxmet,nbins,min,maxmet);
  histo.find2("H_MHTMET")->SetXTitle("MET");
  histo.find2("H_MHTMET")->SetYTitle("MHT");

  histo.make2("H_MPTMET","MPT versus MET",nbins,min,maxmet,nbins,min,maxmet);
  histo.find2("H_MPTMET")->SetXTitle("MET");
  histo.find2("H_MPTMET")->SetYTitle("MPT");

  float jetptmin=0;
  float jetptmax=1000;
  histo.make("H_jetPT1","Jet PT 1",nbins,jetptmin,jetptmax);
  histo.make("H_jetPT2","Jet PT 2",nbins,jetptmin,jetptmax);
  histo.make("H_jetPT3","Jet PT 3",nbins,jetptmin,jetptmax);
  
  histo.make("H_eta1","eta of lead jet",nbins,-maxJetEta_,maxJetEta_);
  histo.make("H_eta1_b","eta of lead b jet",nbins,-maxJetEta_,maxJetEta_);

  int nbins_jets=10;
  string xtitle="N Jets";
  string ytitle="fraction";
  histo.make("H_NJets_ObjDef","H_NJets_ObjDef",nbins_jets,-0.5,nbins_jets-0.5,xtitle);
  histo.make("H_NJets","N Jets",nbins_jets,-0.5,nbins_jets-0.5,xtitle);
  histo.make("H_NJetsLoose","N Jets (loose pT)",nbins_jets,-0.5,nbins_jets-0.5,xtitle);
  //for plotting the number of events passing the trigger as a function of njets
//   histo.make("H_NpassHLT_Jet110","H_NpassHLT_Jet110",nbins_jets,-0.5,nbins_jets-0.5,xtitle);
//   histo.make("H_NpassHLT_DiJetAve70U","H_NpassHLT_DiJetAve70U",nbins_jets,-0.5,nbins_jets-0.5,xtitle);
//   histo.make("H_NpassHLT_MET60","H_NpassHLT_MET60",nbins_jets,-0.5,nbins_jets-0.5,xtitle);
//   histo.make("H_NpassHLT_HT200","H_NpassHLT_HT200",nbins_jets,-0.5,nbins_jets-0.5,xtitle);
//   histo.make("H_NpassAny","H_NpassAny",nbins_jets,-0.5,nbins_jets-0.5,xtitle);
//   histo.make("H_FpassHLT_Jet110","H_FpassHLT_Jet110",nbins_jets,-0.5,nbins_jets-0.5,xtitle,ytitle);
//   histo.make("H_FpassHLT_DiJetAve70U","H_FpassHLT_DiJetAve70U",nbins_jets,-0.5,nbins_jets-0.5,xtitle,ytitle);
//   histo.make("H_FpassHLT_MET60","H_FpassHLT_MET60",nbins_jets,-0.5,nbins_jets-0.5,xtitle,ytitle);
//   histo.make("H_FpassHLT_HT200","H_FpassHLT_HT200",nbins_jets,-0.5,nbins_jets-0.5,xtitle,ytitle);
//   histo.make("H_FpassAny","H_FpassAny",nbins_jets,-0.5,nbins_jets-0.5,xtitle,ytitle);
  
//  histo.select("*");
//  histo.Sumw2();
  
//   histo.select("H_Fpass");
//   histo.SetMaximum(1);   //doesn't work?
  double pi=4*atan(1.0);
  histo.make("H_DeltaPhi1", "angle between jet1 and MET",nbins,0,pi);
  histo.make("H_DeltaPhi2", "angle between jet2 and MET",nbins,0,pi);
  histo.make("H_DeltaPhi3", "angle between jet3 and MET",nbins,0,pi);

  histo.make("H_Nbjets","Number of b jets",nbins_jets,-0.5,nbins_jets-0.5);

  float maxssv=7;
  histo.make("H_SSV_realb","SSV value for MC-matched b jets",50,0,maxssv);
  histo.make("H_SSV_realc","SSV value for MC-matched c jets",50,0,maxssv);
  histo.make("H_SSV_light","SSV value for MC-matched uds+g jets",50,0,maxssv);
  histo.make("H_SSV_other","SSV value for other jets",50,0,maxssv);
  histo.make("H_SSV_unmatched","SSV value for unmatched jets",50,0,maxssv);

  histo.make("H_otherPartonID","Parton ID for other jets", 200,-100,100);

  histo.make("H_NbjetsMC","Number of MC b jets",nbins_jets,-0.5,nbins_jets-0.5);

  histo.make("H_NbjetsMC_noCuts","Number of MC b jets (no cuts)",nbins_jets,-0.5,nbins_jets-0.5);
  histo.make("H_bjetPT1_MC_noCuts","MC pT of lead b jet (no cuts)",nbins,jetptmin,jetptmax);
  histo.make("H_bjetPT2_MC_noCuts","MC pT of 2nd b jet (no cuts)",nbins,jetptmin,jetptmax);
  histo.make("H_bjetPT3_MC_noCuts","MC pT of 3rd b jet (no cuts)",nbins,jetptmin,jetptmax);

  histo.make("H_DeltaPhi_b1b2_MC_noCuts","angle between b1 and b2",nbins,0,pi);
  histo.make("H_DeltaPhi_b1b3_MC_noCuts","angle between b1 and b3",nbins,0,pi);
  histo.make("H_DeltaPhi_b2b3_MC_noCuts","angle between b2 and b3",nbins,0,pi);

  histo.make("H_DeltaPhi_b1j1_MC_noCuts","angle between b1 and j1",nbins,0,pi);
  histo.make("H_DeltaPhi_b1j2_MC_noCuts","angle between b1 and j2",nbins,0,pi);
  histo.make("H_DeltaPhi_b1j3_MC_noCuts","angle between b1 and j3",nbins,0,pi);

  histo.make("H_MinDeltaPhi_trackMET","min DeltaPhi betwen track and MET",nbins,0,1);
  histo.make("H_MHTratio","MHT50 / MHT 30",nbins,0,2);
  histo.make("H_DeltaPhi_b1b2","angle between b1 and b2",nbins,0,pi);
  histo.make("H_DeltaPhi_j1j2","angle between j1 and j2",nbins,0,pi);
  histo.make("H_DeltaPhi_b1j1","angle between b1 and j1",nbins,0,pi);
  //should implement this
  //histo.make("H_DeltaPhi_b1j1_nonb","angle between b1 and j1 (j1 is not a b)",nbins,0,pi);
  histo.make("H_bjetPT1","pT of lead b jet (no cuts)",nbins,jetptmin,jetptmax);
  histo.make("H_bjetPT2","pT of 2nd b jet (no cuts)",nbins,jetptmin,jetptmax);

  histo.make("H_Meff","scalar sum of jet ET and MET",nbins,htmin,htmax);
  histo.make("H_MPT","missing track pT",nbins,min,maxmet);
  histo.make("H_MPT_wide","missing track pT (wide eta)",nbins,min,maxmet);

  histo.make2("H_MET_Meff","MET versus M_eff",nbins,htmin,htmax,nbins,min,maxmet);

  histo.make("H_ST","transverse sphericity",nbins,0,1);
  histo.make("H_STb","transverse sphericity (b jets only)",nbins,0,1);

  histo.make2("H_DeltaPhi_MET", "jet1-MET angle versus MET",nbins,min,maxmet,nbins,0,pi);


  Long64_t npassAllSelection=0;
  Long64_t npassObjectDefinition=0;
  Long64_t npassNJets=0;
  Long64_t npassPreselection=0;
  Long64_t npassMET=0;
  Long64_t npassDeltaPhi=0;
  Long64_t npassbtag=0;

  Long64_t npassBadJetVeto=0;
  Long64_t npassVertex=0;

  Long64_t npassRA2Preselection=0;
  Long64_t npassLeptonVeto=0;

  /***** Darren's code *****/
  //
  // Histogram Declarations
  //
  TH1F* h_ele_pt  = new TH1F("h_ele_pt","Electron p_{T};p_{T}",           60, 0, 300);
  TH1F* h_jet_pt  = new TH1F("h_jet_pt","Jet p_{T};p_{T}",                60, 0, 300);
  TH1F* h_mcpar_pt = new TH1F("h_mcpar_pt","MC particle p_{T};p_{T}",     60, 0, 300);
  TH1F* h_met_pt = new TH1F("h_met_pt","MET p_{T};p_{T}",                 60, 0, 300);
  TH1F* h_muon_pt = new TH1F("h_muon_pt","Muon p_{T};p_{T}",              60, 0, 300);
  TH1F* h_photon_pt = new TH1F("h_photon_pt","Photon p_{T};p_{T}",        60, 0, 300);
  TH1F* h_sc_pt = new TH1F("h_sc_pt","SC p_{T};p_{T}",                    60, 0, 300);
  TH1F* h_track_pt = new TH1F("h_track_pt","Track p_{T};p_{T}",           60, 0, 300);



  int nevents=0;
  int passAll=0;
  double nevents_wgt=0;
  double passAll_wgt=0;


  std::vector<string> hlt_name;
  std::vector<int>    hlt_pass;
  std::vector<string> l1t_algo_name;
  std::vector<int>    l1t_algo_pass;
  std::vector<string> l1t_tech_name;
  std::vector<int>    l1t_tech_pass;

  hlt_name.clear();
  hlt_pass.clear();
  l1t_algo_name.clear();
  l1t_algo_pass.clear();
  l1t_tech_name.clear();
  l1t_tech_pass.clear();


  std::vector<string> selector_name;
  selector_name.push_back("HLT_Ele15_LW_L1R");
  selector_name.push_back("EleCut");
  selector_name.push_back("JetCut");
  selector_name.push_back("METCut");

  selected.clear();
  for( unsigned int i=0; i<selector_name.size(); i++ ) selected.push_back( std::vector<double>(3) );
  for( unsigned int i=0; i<selected.size(); i++ ){
    for( unsigned int j=0; j<selected[0].size(); j++ ) selected[i][j]=0.;
  }



  // Loop over events
  for( ev.toBegin(); !ev.atEnd(); ++ev) {

    if (nevents%5000 == 0) cout << "Event # "<<nevents <<endl;
    try {

      fwlite::Handle<CUeventCollection> h_event;
      h_event.getByLabel(ev,"CUproducer");

      fwlite::Handle<CUelectronCollection> h_electrons;
      h_electrons.getByLabel(ev,"CUproducer","cleanLayer1Electrons");
      CUelectronCollection const &electrons = *h_electrons;

      fwlite::Handle<CUjetCollection> h_jets;
      h_jets.getByLabel(ev,"CUproducer","cleanLayer1JetsAK5");
      CUjetCollection const &jets = *h_jets;

      fwlite::Handle<CUmetCollection> h_met;
      h_met.getByLabel(ev,"CUproducer","layer1METsAK5");

      fwlite::Handle<CUmcparticleCollection> h_mcparticles;
      h_mcparticles.getByLabel(ev,"CUproducer","MCstatus3");
      CUmcparticleCollection const &mcparticles = *h_mcparticles;

      fwlite::Handle<CUmuonCollection> h_muons;
      h_muons.getByLabel(ev,"CUproducer","cleanLayer1Muons");
      CUmuonCollection const &muons = *h_muons;

      fwlite::Handle<CUphotonCollection> h_photons;
      h_photons.getByLabel(ev,"CUproducer","cleanLayer1Photons");
      CUphotonCollection const &photons = *h_photons;

      fwlite::Handle<CUsuperclusterCollection> h_superclusters;
      h_superclusters.getByLabel(ev,"CUproducer","corHybridSCandMulti5x5WithPreshower");
      CUsuperclusterCollection const &scs = *h_superclusters;

      fwlite::Handle<CUtrackCollection> h_tracks;
      h_tracks.getByLabel(ev,"CUproducer","generalTracks");
      CUtrackCollection const &tracks = *h_tracks;

      fwlite::Handle<CUtriggerCollection> h_hlt;
      h_hlt.getByLabel(ev,"CUproducer","HLT");
      CUtriggerCollection const &hlt = *h_hlt;

      fwlite::Handle<CUtriggerCollection> h_l1t_algo;
      h_l1t_algo.getByLabel(ev,"CUproducer","L1Talgo");
      CUtriggerCollection const &l1t_algo = *h_l1t_algo;

      fwlite::Handle<CUtriggerCollection> h_l1t_tech;
      h_l1t_tech.getByLabel(ev,"CUproducer","L1Ttech");
      CUtriggerCollection const &l1t_tech = *h_l1t_tech;

      fwlite::Handle<CUprimaryvertexCollection> h_privertex;
      h_privertex.getByLabel(ev,"CUproducer","offlinePrimaryVertices");
      CUprimaryvertexCollection const &vertices = *h_privertex;

      fwlite::Handle<CUskimbitsCollection> h_skimbits;
      h_skimbits.getByLabel(ev,"CUproducer","SkimBits");


      double weight = h_event->front().weight;
      double sample = h_event->front().sample;
      double run  = h_event->front().run;
      double evt  = h_event->front().evt;
      double lumi = h_event->front().lumi;
      double bsx  = h_event->front().BSx;
      double bsy  = h_event->front().BSy;
      double bsz  = h_event->front().BSz;
      double numPV = h_event->front().numPV;

      double metPT  = h_met->front().pt;
      double genMET = h_met->front().genPT;
      double sumET  = h_met->front().sumET;

      int Ele15 = h_skimbits->front().HLT_Ele15_LW_L1R;
      int Jet15 = h_skimbits->front().HLT_Jet15U;

      double wgt = weight;

      nevents++;
      nevents_wgt+=wgt;

      if( nevents==1 ){
	hlt_name.resize(hlt.size());
	hlt_pass.resize(hlt.size());
	l1t_algo_name.resize(l1t_algo.size());
	l1t_algo_pass.resize(l1t_algo.size());
	l1t_tech_name.resize(l1t_tech.size());
	l1t_tech_pass.resize(l1t_tech.size());
      }

      if( false && nevents<=10 ){
	std::cout << " ------------------------------------------- " << std::endl;
	std::cout << "   Run: " << run << ", Lumi: " << lumi << ", Event: " << evt << std::endl;
	std::cout << "   sample: " << sample << ", weight of event: " << weight << std::endl;
	std::cout << "   BS at (" << bsx << "," << bsy << "," << bsz << ")" << std::endl;
	std::cout << "   NumPV: " << numPV << std::endl;
	std::cout << "   met: " << metPT << ", genMET: " << genMET << ", sumET: " << sumET << std::endl;
	std::cout << "   Pass Ele15? " << Ele15 << ", Pass Jet15? " << Jet15 << std::endl;
	if( nevents==10 ) std::cout << " ------------------------------------------- " << std::endl;
      }

      h_met_pt->Fill(metPT);

      //if the jets object is smaller than the min number of good jets, just give up
      if (jets.size() < minGoodJets_ ) continue;

      //before doing anything else, look at some MC truth quantities
      int ngoodMCbjets_nocuts=0;
      Float_t bjet_pT1=-1, bjet_pT2=-2, bjet_pT3=-3;
      Float_t jjet_pT1=-1, jjet_pT2=-2, jjet_pT3=-3;
      JetIter bj1=jets.end(),bj2=jets.end(),bj3=jets.end(),jj1=jets.end(),jj2=jets.end(),jj3=jets.end();
      for( JetIter jet = jets.begin(); jet != jets.end(); ++jet ) {

	if (jet->genPartonId == 5 || jet->genPartonId == -5 ) {
	  ngoodMCbjets_nocuts++;
	  if ( jet->genPartonPT > bjet_pT1) {
	    bjet_pT3 = bjet_pT2;
	    bj3=bj2;

	    bjet_pT2 = bjet_pT1;
	    bj2=bj1;

	    bjet_pT1 = jet->genPartonPT;
	    bj1=jet;
	  }
	  else if ( jet->genPartonPT > bjet_pT2) {
	    bjet_pT3 = bjet_pT2;
	    bj3=bj2;

	    bjet_pT2 = jet->genPartonPT;
	    bj2=jet;
	  }
	  else if ( jet->genPartonPT > bjet_pT3) {
	    bjet_pT3 = jet->genPartonPT;
	    bj3=jet;
	  }
	}

	if ( jet->genPartonPT > jjet_pT1) {
	  jjet_pT3 = jjet_pT2;
	  jj3=jj2;
	  
	  jjet_pT2 = jjet_pT1;
	  jj2=jj1;
	  
	  jjet_pT1 = jet->genPartonPT;
	  jj1=jet;
	}
	else if ( jet->genPartonPT > jjet_pT2) {
	  jjet_pT3 = jjet_pT2;
	  jj3=jj2;

	  jjet_pT2 = jet->genPartonPT;
	  jj2=jet;
	}
	else if ( jet->genPartonPT > jjet_pT3) {
	  jjet_pT3 = jet->genPartonPT;
	  jj3=jet;
	}
      }

      histo["H_NbjetsMC_noCuts"]->Fill(ngoodMCbjets_nocuts);
      if (ngoodMCbjets_nocuts >= 2) {
	histo["H_bjetPT1_MC_noCuts"]->Fill(bjet_pT1);
	histo["H_bjetPT2_MC_noCuts"]->Fill(bjet_pT2);
	histo["H_bjetPT3_MC_noCuts"]->Fill(bjet_pT3);

	histo["H_DeltaPhi_b1b2_MC_noCuts"]->Fill(getDeltaPhiMC(bj1,bj2));
	if (ngoodMCbjets_nocuts >= 3) {
	  histo["H_DeltaPhi_b1b3_MC_noCuts"]->Fill(getDeltaPhiMC(bj1,bj3));
	  histo["H_DeltaPhi_b2b3_MC_noCuts"]->Fill(getDeltaPhiMC(bj2,bj3));
	}
      }
      if (ngoodMCbjets_nocuts >= 1) {
	if (jj1 !=jets.end()) histo["H_DeltaPhi_b1j1_MC_noCuts"]->Fill(getDeltaPhiMC(bj1,jj1));
	if (jj2 !=jets.end()) histo["H_DeltaPhi_b1j2_MC_noCuts"]->Fill(getDeltaPhiMC(bj1,jj2));
	if (jj3 !=jets.end()) histo["H_DeltaPhi_b1j3_MC_noCuts"]->Fill(getDeltaPhiMC(bj1,jj3));
      }

      bool badjetveto=false;
      //recalculate MHT and HT using the current jet cuts
      Float_t myHT=0,myHTLoose=0;
      Float_t myMHTx=0,myMHTxLoose=0;
      Float_t myMHTy=0,myMHTyLoose=0;
      UInt_t ngoodjets=0,nloosejets=0;
      Int_t nbjets=0;
      JetIter ijet1= jets.end(), ijet2= jets.end(), ijet3= jets.end();
      Float_t jet_pT1=-1, jet_pT2=-2, jet_pT3=-3;
      JetIter ijet1P= jets.end(), ijet2P= jets.end(), ijet3P= jets.end();
      Float_t jet_pT1P=-1, jet_pT2P=-2, jet_pT3P=-3;

      Float_t meff=0;

      Float_t S11=0,S12=0, S22=0;
      Float_t S11b=0,S12b=0, S22b=0;

      bjet_pT1=-1; bjet_pT2=-2;
      bj1=jets.end(); bj2=jets.end();

      // ============ this is the main loop over jets ===============

      for( JetIter jet = jets.begin(); jet != jets.end(); ++jet ) {
	h_jet_pt->Fill(jet->pt);

	//find the leading 3 jets; don't care about passing jet cuts
	if ( jet->pt > jet_pT1P) {
	  jet_pT3P = jet_pT2P;
	  ijet3P=ijet2P;

	  jet_pT2P = jet_pT1P;
	  ijet2P=ijet1P;

	  jet_pT1P = jet->pt;
	  ijet1P=jet;
	}
	else if ( jet->pt > jet_pT2P ) {
	  jet_pT3P = jet_pT2P;
	  ijet3P=ijet2P;

	  jet_pT2P = jet->pt;
	  ijet2P=jet;
	}
	else if (jet->pt > jet_pT3P) {
	  jet_pT3P = jet->pt;
	  ijet3P=jet;
	}

	bool passCuts = passJetCuts(jet);
	//check if there is a hard jet that fails some other cuts
	if ( (jet->pt > minJetPT_) && !passCuts ) badjetveto=true;
	if ( passCuts ) { //these are the only jets (more or less) to be used in subsequent analysis
	  ngoodjets++;	  
	  
	  meff += jet->et;

	  myHT += jet->pt;
	  myMHTx += jet->px;
	  myMHTy += jet->py;
	  
	  S11 += jet->px * jet->px;
	  S12 += jet->px * jet->py;
	  S22 += jet->py * jet->py;

	  //find the first 3 leading jets
	  if ( jet->pt > jet_pT1) {
	    jet_pT3 = jet_pT2;
	    ijet3=ijet2;
	    
	    jet_pT2= jet_pT1;
	    ijet2=ijet1;

	    jet_pT1 = jet->pt;
	    ijet1=jet;
	  }
	  else if ( jet->pt > jet_pT2 ) {
	    jet_pT3 = jet_pT2;
	    ijet3=ijet2;
	    jet_pT2 = jet->pt;
	    ijet2=jet;
	  }
	  else if (jet->pt > jet_pT3) {
	    jet_pT3 = jet->pt;
	    ijet3=jet;
	  }
	}

	if ( passBJetCuts(jet) ) {
	  nbjets++;
	  S11b += jet->px * jet->px;
	  S12b += jet->px * jet->py;
	  S22b += jet->py * jet->py;

	  //find the lead 2 b jets
	  if ( jet->pt > bjet_pT1 ) {
	    bjet_pT2 = bjet_pT1;
	    bj2=bj1;
	    bjet_pT1 = jet->pt;
	    bj1=jet;
	  }
	  else if (jet->pt > bjet_pT2) {
	    bjet_pT2 = jet->pt;
	    bj2=jet;
	  }
	  
	}

	if (cuts_=="JMT" && passLooseJetCuts(jet) ) {
	  
	  myHTLoose += jet->pt;
	  myMHTxLoose += jet->px;
	  myMHTyLoose += jet->py;
	  
	  nloosejets++;
	}
      }

      if (cuts_=="RA2" && !passRA2JetPreselection(ijet1P,ijet2P,ijet3P)) continue;
      npassRA2Preselection++;

      if (cuts_ !="RA2") {
	if ( ngoodjets<minGoodJets_) continue; //require a minimum number of good jets
	npassNJets++;

	if (badJetVeto_ && badjetveto)  continue;
	npassBadJetVeto++;
      }


      //reject events with no primary vertex
      bool vtxReject=false;
      if (requirePrimaryVertex_) {
	vtxReject=true;
	for ( VtxIter vtx = vertices.begin(); vtx != vertices.end(); ++vtx ) {
	  if (vtx->isValid==1 && vtx->isFake==0) { //correct to use 1 and 0?
	    vtxReject=false;
	    break; //avoid wasting time once we've passed the cut
	  }
	}
      }
      if (vtxReject)  continue;
      npassVertex++;

      Int_t nelectrons=0;
      for( EleIter ele = electrons.begin(); ele != electrons.end(); ++ele ){
	h_ele_pt->Fill(ele->pt);

	if (ele->pt < 15) continue;
	if (fabs(ele->eta) > 2.5) continue;

	if ( ele->tkD0bs > 0.2 ) continue;
		     //&& (fabs(ele->eta)>1.567 || fabs(ele->eta)<1.47) );

	//not sure which trackIso variable to use here
	double reliso = (ele->trackIso + ele->ecalIsoDR04 + ele->hcalIsoDR04 )/ele->pt;

	if ( reliso >0.5 ) continue;

	if ( ele->IDLoose != 1) continue;

	//if we get to here it must be good
	nelectrons++;
      }

      //electron veto
      if ( leptonVeto_ && nelectrons >0 ) continue;

      //note that RA2 also defines a control sample of muons with different cuts
      Int_t nmuons=0;
      for( MuonIter muon = muons.begin(); muon != muons.end(); ++muon ){
	h_muon_pt->Fill(muon->pt);

	if ( muon->IDGMPTight != 1) continue;

	if (muon->pt < 10) continue;
	if (fabs(muon->eta) >2.4) continue;

	double reliso = (muon->trackIsoDR03 + muon->ecalIsoDR03 + muon->hcalIsoDR03 )/muon->pt;
	if ( reliso >0.1) continue;

	if ( muon->comNormChi2 > 10 ) continue; //is this the right thing to use?
	if (muon->tkD0bs >0.2) continue;

	if (muon->tkNumValidHits <11) continue;

	nmuons++;
      }
      if ( leptonVeto_ && nmuons >0) continue; //muon veto

      npassLeptonVeto++;

      Int_t nphotons=0;
      for( PhotonIter photon = photons.begin(); photon != photons.end(); ++photon ){
	h_photon_pt->Fill(photon->et);
	nphotons++;
      }
      if ( photonVeto_ && nphotons >0 ) continue; //photon veto

//       if (passHLT_Jet110 ) histo["H_NpassHLT_Jet110"]->Fill(ngoodjets);
//       if (passHLT_DiJetAve70U ) histo["H_NpassHLT_DiJetAve70U"]->Fill(ngoodjets);
//       if (passHLT_MET60 ) histo["H_NpassHLT_MET60"]->Fill(ngoodjets);
//       if (passHLT_HT200 ) histo["H_NpassHLT_HT200"]->Fill(ngoodjets);
//      if (passHLT_Jet110 || passHLT_DiJetAve70U || passHLT_MET60 || passHLT_HT200) 
// 	histo["H_NpassAny"]->Fill(ngoodjets);

      if (useMcInfo_) {
	for( MCparIter mcpar = mcparticles.begin(); mcpar != mcparticles.end(); ++mcpar ){
	  h_mcpar_pt->Fill(mcpar->pt);
	}
      }

      if (useScInfo_) {
	for( SCIter sc = scs.begin(); sc != scs.end(); ++sc ){
	  h_sc_pt->Fill(sc->et);
	}
      }

      //ok, now we have applied some jets cuts and vetoed leptons

      int hlt_i=0;
      int l1t_algo_i=0;
      int l1t_tech_i=0;
      if (useTriggerInfo_) {
	for( TrigIter hltbit = hlt.begin(); hltbit != hlt.end(); ++hltbit ){
	  hlt_pass[hlt_i] += hltbit->pass;
	  if( nevents==1 ) hlt_name[hlt_i] = hltbit->name;
	  hlt_i++;
	}
	
	for( TrigIter l1t_algobit = l1t_algo.begin(); l1t_algobit != l1t_algo.end(); ++l1t_algobit ){
	  l1t_algo_pass[l1t_algo_i] += l1t_algobit->pass;
	  if( nevents==1 ) l1t_algo_name[l1t_algo_i] = l1t_algobit->name;
	  l1t_algo_i++;
	}
	
	for( TrigIter l1t_techbit = l1t_tech.begin(); l1t_techbit != l1t_tech.end(); ++l1t_techbit ){
	  l1t_tech_pass[l1t_tech_i] += l1t_techbit->pass;
	  if( nevents==1 ) l1t_tech_name[l1t_tech_i] = l1t_techbit->name;
	  l1t_tech_i++;
	}
      }

      //the "object definition" cuts are done by now, but not the "preselection"
      //note that this is the SUSY-09-001 language, not the RA2 language
      histo["H_jetPT1"]->Fill(ijet1P->pt);
      histo["H_jetPT2"]->Fill(ijet2P->pt);
      histo["H_jetPT3"]->Fill(ijet3P->pt);
      histo["H_NJets_ObjDef"]->Fill(ngoodjets);
      npassObjectDefinition++;
      //preselection
      //RA2 uses the three lead jets, period
      if (cuts_=="RA2") if (!(passRA2JetCuts(ijet1P,ijet2P,ijet3P))) continue;
      else      if (!passPreselection(ijet1,ijet2)) continue;
      npassPreselection++;

      //========================<<<<>>>><<<<>>>><<<<>>>><<<<>>>>==========================
      //Now we have applied basic jet cuts and gotten rid of leptons
      //plus we have done some event cleanup
      //now make plots!
      histo["H_NJets"]->Fill(ngoodjets);
      histo["H_NJetsLoose"]->Fill(nloosejets);

      //fill HT for events that pass preselection
      histo["H_HT"]->Fill(myHT);

      histo["H_Meff"]->Fill(meff + metPT);

      histo.find2("H_MHTMET")->Fill(metPT,sqrt(myMHTx*myMHTx + myMHTy*myMHTy));
      histo.find("H_MHTratio")->Fill(sqrt(myMHTx*myMHTx + myMHTy*myMHTy) / sqrt(myMHTxLoose*myMHTxLoose + myMHTyLoose*myMHTyLoose));

      histo.find2("H_MET_Meff")->Fill(meff+metPT,metPT);

      //b jet variables
      if (nbjets>=2) {
	histo.find("H_DeltaPhi_b1b2")->Fill(getDeltaPhi(bj1,bj2));
	histo.find("H_bjetPT2")->Fill(bjet_pT2);

	TMatrixT<float> SmatB(2,2);
	SmatB[0][0]=S11b;
	SmatB[0][1]=S12b;
	SmatB[1][0]=S12b;
	SmatB[1][1]=S22b;
	TVectorT<float> eigenvaluesB;
	TMatrixT<float> eigenvectorsB = SmatB.EigenVectors(eigenvaluesB);
	histo["H_STb"]->Fill(2*eigenvaluesB[1] / (eigenvaluesB[0]+eigenvaluesB[1]));

      }
      if (nbjets>=1) {
	histo.find("H_bjetPT1")->Fill(bjet_pT1);
	histo["H_DeltaPhi_b1j1"]->Fill( getDeltaPhi(ijet1,bj1));
	histo["H_eta1_b"]->Fill( bj1->eta);
      }
      histo["H_Nbjets"]->Fill(nbjets);

      histo["H_eta1"]->Fill( ijet1->eta);
      TMatrixT<float> Smat(2,2);
      Smat[0][0]=S11;
      Smat[0][1]=S12;
      Smat[1][0]=S12;
      Smat[1][1]=S22;
      TVectorT<float> eigenvalues;
      TMatrixT<float> eigenvectors = Smat.EigenVectors(eigenvalues);

      histo["H_ST"]->Fill(2*eigenvalues[1] / (eigenvalues[0]+eigenvalues[1]));

      if (useTrackInfo_) {
	double MPTx=0;
	double MPTy=0;
	double MPTxLoose=0;
	double MPTyLoose=0;
	double minDeltaPhi=99;
	for ( TrackIter track = tracks.begin(); track != tracks.end(); ++track ) {
	  h_track_pt->Fill(track->pt);

	  //there are no isolation variables stored for tracks ... need to compute it myself?
	  //FIXME..should move these cut values into variables
	  if (track->pt > 10) {
	    //should also study DeltaR between tracks and jets
	    double dp=   getDeltaPhi(track,h_met->front());
	    if (dp < minDeltaPhi) minDeltaPhi=dp;
	    //MPT
	    if ( fabs(track->eta) < maxJetEta_ ) {
	      MPTx += track->px;
	      MPTy += track->py;
	    }
	    if ( fabs(track->eta) < 5 ) {
	      MPTxLoose += track->px;
	      MPTyLoose += track->py;
	    }
	  }
	}
	histo["H_MinDeltaPhi_trackMET"]->Fill(minDeltaPhi);
	histo["H_MPT"]->Fill( sqrt( MPTx*MPTx + MPTy*MPTy));
	histo.find2("H_MPTMET")->Fill(metPT, sqrt( MPTx*MPTx + MPTy*MPTy));
	histo["H_MPT_wide"]->Fill( sqrt( MPTxLoose*MPTxLoose + MPTyLoose*MPTyLoose));
      }

      histo["H_DeltaPhi1"]->Fill( getDeltaPhi(ijet1,h_met->front()));
      histo["H_DeltaPhi2"]->Fill( getDeltaPhi(ijet2,h_met->front()));

      histo["H_DeltaPhi_j1j2"]->Fill( getDeltaPhi(ijet1,ijet2));
      if (cuts_=="RA2") histo["H_DeltaPhi3"]->Fill( getDeltaPhi(ijet3,h_met->front()));

      histo.find2("H_DeltaPhi_MET")->Fill( metPT,  getDeltaPhi(ijet1,h_met->front()));

      int ngoodMCbjets=0;
      for ( JetIter jet = jets.begin(); jet != jets.end(); ++jet ) {
	if (passJetCuts(jet)) {
	  if (jet->genPartonId == 5 || jet->genPartonId == -5 ) {
	    histo["H_SSV_realb"]->Fill(jet->btagSecVertex );
	  }
	  else if (jet->genPartonId == 4 || jet->genPartonId == -4) {
	    histo["H_SSV_realc"]->Fill(jet->btagSecVertex );
	  }
	  else if ( fabs(jet->genPartonId) == 3 
		    || fabs(jet->genPartonId) == 2
		    || fabs(jet->genPartonId) == 1
		    || fabs(jet->genPartonId) == 21  ) {
	    histo["H_SSV_light"]->Fill(jet->btagSecVertex );
	  }
	  else if ( fabs(jet->genPartonId) == 999 ) {
	    histo["H_SSV_unmatched"]->Fill(jet->btagSecVertex );
	  }
	  else {
	    histo["H_SSV_other"]->Fill(jet->btagSecVertex );
	    histo["H_otherPartonID"]->Fill( jet->genPartonId);
	  }
	}

	if ( jet->genPartonPT >= minJetPT_ && fabs(jet->genPartonEta)<=maxJetEta_ 
	     && fabs(jet->genPartonId)==5 ) {
	  ngoodMCbjets++;
	}
      }
      histo["H_NbjetsMC"]->Fill(ngoodMCbjets);

      //here is where we start making more cuts. I don't think this matters for the moment.
      //====================^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^=========================
      if (metPT < minMET_) continue;
      npassMET++;

 
      if (!passDeltaPhi(ijet1,h_met->front())) continue;
      if (!passDeltaPhi(ijet2,h_met->front())) continue;
      if (cuts_=="RA2" && !passDeltaPhi(ijet3,h_met->front())) continue;
      npassDeltaPhi++;

      if (nbjets < minbjets_) continue;
      npassbtag++;
      
      npassAllSelection++;

      /*
      //Darren's code
      bool HLT_Ele15_LW_L1R = ( Ele15 );
      bool EleCut = ( nelectrons==1 );
      bool JetCut = ( ngoodjets>2 );
      bool METCut = ( metPT>100 );

      mydecision.push_back(HLT_Ele15_LW_L1R);
      mydecision.push_back(EleCut);
      mydecision.push_back(JetCut);
      mydecision.push_back(METCut);


      if( mydecision.size()!=selector_name.size() ){
	std::cout << " ERROR!! mydecision.size()!=selector_name.size()" << std::endl;
	break;
      }

      selected = CUselect::caculateSelector(mydecision,selected,1);

      mydecision.clear();
      */
    }// end try
    catch(std::exception& e) {
      std::cerr << " ==> caught exception " << e.what() << std::endl;
      continue;
    }

  } // end loop over events

  fout.Write();
  fout.Close();

  std::cout << " *********************************************************** " << std::endl;
  std::cout << "   Total Number of Events = " << nevents_wgt << ", (" << nevents << " unwgt)" << std::endl;
  std::cout << "                  passAll = " << passAll_wgt << ", (" << passAll << " unwgt)" << std::endl;
  std::cout << "               efficiency = " << passAll_wgt/nevents_wgt << ", (" << double(passAll)/double(nevents) << " unwgt)" << std::endl;
  std::cout << " *********************************************************** " << std::endl;

  CUselect::printSelectorResults(selector_name,selected,nevents_wgt);

  std::cout << " *********************************************************** " << std::endl;
  std::cout << " HLT Trigger Efficiency " << std::endl;
  std::cout << " ---------------------- " << std::endl;
  for( unsigned int i=0; i<hlt_pass.size(); i++ ){
    std::cout << "\t" << i << "\t eff = "<< (double)hlt_pass[i]/(double)nevents << "\t" << hlt_name[i] << std::endl;
  }
  std::cout << " *********************************************************** " << std::endl;
  std::cout << " L1T Algo Trigger Efficiency " << std::endl;
  std::cout << " --------------------------- " << std::endl;
  for( unsigned int i=0; i<l1t_algo_pass.size(); i++ ){
    std::cout << "\t" << i << "\t eff = "<< (double)l1t_algo_pass[i]/(double)nevents << "\t" << l1t_algo_name[i] << std::endl;
  }
  std::cout << " *********************************************************** " << std::endl;
  std::cout << " L1T Tech Trigger Efficiency " << std::endl;
  std::cout << " --------------------------- " << std::endl;
  for( unsigned int i=0; i<l1t_tech_pass.size(); i++ ){
    std::cout << "\t" << i << "\t eff = "<< (double)l1t_tech_pass[i]/(double)nevents << "\t" << l1t_tech_name[i] << std::endl;
  }
  std::cout << " *********************************************************** " << std::endl;


  cout<<"------"<<endl;

  cout<<"Number total                            = "<<nevents<<endl;
  cout<<"Number pass RA2 Jet Preselection        = "<<npassRA2Preselection<<"\t"<<100*double(npassRA2Preselection)/double(nevents)<<endl;
  cout<<"Number pass N jets                      = "<<npassNJets<<"\t"<<100*double(npassNJets)/double(nevents)<<endl;
  cout<<"Number pass bad jet veto                = "<<npassBadJetVeto<<"\t"<<100*double(npassBadJetVeto)/double(nevents)<<endl;
  cout<<"Number pass pri vtx requirement         = "<<npassVertex<<"\t"<<100*double(npassVertex)/double(nevents)<<endl;
  cout<<"Number pass Lepton Veto                 = "<<npassLeptonVeto<<"\t"<<100*double(npassLeptonVeto)/double(nevents)<<endl;
  cout<<"Number pass obj def                     = "<<npassObjectDefinition<<"\t"<<100*double(npassObjectDefinition)/double(nevents)<<endl;
  cout<<"Number pass preselection/RA2 Jet Cuts   = "<<npassPreselection<<"\t"<<100*double(npassPreselection)/double(nevents)<<endl;
  cout<<" -- make plots --"<<endl;
  cout<<"Number pass MET Cut                     = "<<npassMET<<"\t"<<100*double(npassMET)/double(nevents)<<endl;
  cout<<"Number pass DeltaPhi                    = "<<npassDeltaPhi<<"\t"<<100*double(npassDeltaPhi)/double(nevents)<<endl;
  cout<<"Number pass # B jets                    = "<<npassbtag<<"\t"<<100*double(npassbtag)/double(nevents)<<endl;
  cout<<"Number pass all cuts                    = "<<npassAllSelection<<"\t"<<100*double(npassAllSelection)/double(nevents)<<endl;


}

