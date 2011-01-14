//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jan 14 12:05:44 2011 by ROOT version 5.22/00d
// from TTree tree/tree
// found on file: /cu1/joshmt/BasicNtuples/V00-01-05/QCD-Pt100to250-madgraph/BasicNtuple_10_1_Vkj.root
//////////////////////////////////////////////////////////

#ifndef basicLoop_h
#define basicLoop_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
// ========================================== begin
//this file is in CVS here: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/joshmt/MiscUtil.cxx?view=log
#include "/afs/cern.ch/user/j/joshmt/root/util/MiscUtil.cxx"
//this code will have to be regenerated when changing the ntuple structure
//custom code is marked with these 'begin' and 'end' markers
// ---- this version is compatible with ntuple tag: V00-01-04 ----
#include <iostream>
#include <vector>
#include <set>

/*
some proto-documentation:

To implement a new Cut Scheme:

-- add a text string to CutSchemeNames_ and a corresponding value to the CutScheme enum
-- Add a block to setCutScheme(). This defines the order of the cut flow
-- Add a block to cutRequired(). This defines whether the cut is actually used in this scheme
-- for any cut where the default value stored in the ntuple is not valid, add if statements in passCut() that
can handle this cut properly. (see note below)


notes on the schemes:
RA2 -- the legacy scheme. still uses some info precomputed in the ntuple.
Sync1 -- cuts for first sync exercise. Recomputes (almost?) every cut.
Baseline0 -- *must* recompute every cut. "Baseline" analysis cuts as laid out by me on the joint Twiki.

Other notes:
a lot of the functions are still doing computations using the tight jet list stored in the ntuple.
This is rather dangerous. In particular, the eta range of the ntuple tight jets does not agree 100% with the kBaseline0 tight jets.

*/

//avoid spaces and funny characters!
const char *CutSchemeNames_[]={"RA2", "Sync1", "Baseline0"};
const char *METTypeNames_[]={"MHT", "MET",  "tcMET", "pfMET"};
const char *METRangeNames_[]={"med",  "high", "wide", "medhigh"}; //'no cut' is not given, because the cut can always be skipped!

const char *jetTypeNames_[]={"calo","PF"}; //no JPT in ntuple for now
const char *leptonTypeNames_[]={"RegLep","PFLep"};

const char *dpTypeNames_[]={"DeltaPhi", "minDP",  "MPT", "DPSync1", "minDPinv"};

//in 1/pb
const double lumi=36.143; //Don's number for 386 Nov4ReReco
// ========================================== end

class basicLoop {
public :
  // ========================================== begin
  enum CutScheme {kRA2=0, kSync1, kBaseline0, nCutSchemes}; //cut schemes must be implemented in setCutScheme,etc and in the list above
  CutScheme theCutScheme_;
  enum METType {kMHT=0, kMET, ktcMET, kpfMET};
  METType theMETType_;
  enum METRange {kMedium=0, kHigh, kWide, kMedhigh};
  METRange theMETRange_;
  enum jetType {kCalo=0, kPF};
  jetType theJetType_;
  enum leptonType {kNormal=0, kPFLeptons};
  leptonType theLeptonType_;
  enum dpType {kDeltaPhi=0, kminDP, kMPT, kDPSync1, kminDPinv};
  dpType theDPType_;
  unsigned int nBcut_;
  //these are an extension of the ele/mu veto
  //If the cutEleVeto or cutMuVeto is set, then there must be exactly this number of e/mu in the event
  //by default, init to 0 (veto)
  //setting to 1 inverts the veto
  int ne_;
  int nmu_;

  bool isData_; //set in ctor

  //any trigger on this list will be considered OK for passing the trigger cut
  //key is the trigger name ; value is the position in the passTrigger variable in the ntuple
  std::map<TString, int>  triggerList_; //contains only triggers used for determining the signal
  std::map<TString, int>  triggerListAll_; //this contains all triggers
  //can use hltPrescale to figure out if the trigger exists in a given event/sample

  std::vector<TString> cutTags_;
  std::map<TString, TString> cutNames_; //key is a cutTag. these should be "human readable" but should not have any spaces
  std::map<TString, int> cutMap_; //key is a cutTag, value is the position in the ntuple CutResults
  std::vector<TString> ignoredCut_; //allow more than 1 ignored cut!
  //if theCutFlow changes, be sure to change cutnames_ as well

  std::set<jmt::eventID> specifiedEvents_;

  enum TopDecayCategory {kTTbarUnknown=0,kAllLeptons=1,kAllHadronic=2,kOneElectron=3,kOneMuon=4,kOneTauE=5,kOneTauMu=6,kOneTauHadronic=7,kAllTau=8,kTauPlusLepton=9, nTopCategories=10};

  //tight jet info
  //no longer in ntuple, so we will create them on the fly for each event
  //I will fill these for the jet type selected by theJetType_
  vector<float>   jetPt;
  vector<float>   jetEta;
  vector<float>   jetPhi;
  vector<int>     jetFlavor;
  vector<float>   jetBTagDisc_trackCountingHighPurBJetTags;
  vector<float>   jetBTagDisc_trackCountingHighEffBJetTags;
  vector<float>   jetBTagDisc_simpleSecondaryVertexHighEffBJetTags;
  vector<float>   jetBTagDisc_simpleSecondaryVertexHighPurBJetTags;
  vector<float>   jetBTagDisc_simpleSecondaryVertexBJetTags;

  //pointers to the jet info of the default jet type
  vector<int>     *tightJetIndex;
  vector<int>     *looseJetIndex;
  vector<float>   *loosejetPt;
  vector<float>   *loosejetEt;
  vector<float>   *loosejetEta;
  vector<float>   *loosejetPhi;
  vector<bool>    *loosejetPassLooseID;
  vector<float>   *loosejetEnergyFracHadronic;
  vector<int>     *loosejetFlavor;
  vector<int>     *loosejetGenParticlePDGId;
  vector<float>   *loosejetInvisibleEnergy;
  vector<float>   *loosejetBTagDisc_trackCountingHighPurBJetTags;
  vector<float>   *loosejetBTagDisc_trackCountingHighEffBJetTags;
  vector<float>   *loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags;
  vector<float>   *loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags;
  vector<float>   *loosejetBTagDisc_simpleSecondaryVertexBJetTags;

  //new in V00-01-02
  vector<int>     *loosejetGenPt;
  vector<int>     *loosejetGenEta;
  vector<int>     *loosejetGenPhi;
  vector<int>     *loosejetNSV;
  vector<int>     *loosejetNTracks; //only 00-01-03
  vector<float>   *loosejetSVUnWeightedMass;
  vector<float>   *loosejetSVWeightedMass;
  vector<float>   *loosejetSVUnWeightedLifetime;
  vector<float>   *loosejetSVWeightedLifetime;
  vector<float>   *loosejetSVUnWeightedCosTheta;
  vector<float>   *loosejetSVWeightedCosTheta;

  vector<float>   *badjetPt;
  vector<float>   *badjetEta;
  vector<float>   *badjetPhi;
  // ========================================== end

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   ULong64_t       runNumber;
   ULong64_t       lumiSection;
   ULong64_t       eventNumber;
   Float_t         bsx;
   Float_t         bsy;
   Float_t         bsz;
   vector<bool>    *cutResults;
   vector<bool>    *passTrigger;
   vector<unsigned int> *hltPrescale;
   Int_t           SUSYtriggerIndex;
   vector<bool>    *pv_isFake;
   vector<float>   *pv_z;
   vector<float>   *pv_rho;
   vector<float>   *pv_chi2;
   vector<float>   *pv_ndof;
   vector<int>     *tightJetIndex_calo;
   vector<int>     *looseJetIndex_calo;
   vector<float>   *loosejetPt_calo;
   vector<float>   *loosejetEt_calo;
   vector<float>   *loosejetEta_calo;
   vector<float>   *loosejetPhi_calo;
   vector<bool>    *loosejetPassLooseID_calo;
   vector<bool>    *loosejetPassTightID_calo;
   vector<float>   *loosejetEnergyFracHadronic_calo;
   vector<int>     *loosejetFlavor_calo;
   vector<int>     *loosejetGenPt_calo;
   vector<int>     *loosejetGenEta_calo;
   vector<int>     *loosejetGenPhi_calo;
   vector<int>     *loosejetGenParticlePDGId_calo;
   vector<float>   *loosejetInvisibleEnergy_calo;
   vector<int>     *loosejetNTracks_calo;
   vector<int>     *loosejetNSV_calo;
   vector<float>   *loosejetSVUnWeightedMass_calo;
   vector<float>   *loosejetSVWeightedMass_calo;
   vector<float>   *loosejetSVUnWeightedLifetime_calo;
   vector<float>   *loosejetSVWeightedLifetime_calo;
   vector<float>   *loosejetSVUnWeightedCosTheta_calo;
   vector<float>   *loosejetSVWeightedCosTheta_calo;
   vector<float>   *loosejetBTagDisc_trackCountingHighPurBJetTags_calo;
   vector<float>   *loosejetBTagDisc_trackCountingHighEffBJetTags_calo;
   vector<float>   *loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags_calo;
   vector<float>   *loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags_calo;
   vector<float>   *loosejetBTagDisc_simpleSecondaryVertexBJetTags_calo;
   vector<float>   *badjetPt_calo;
   vector<float>   *badjetEta_calo;
   vector<float>   *badjetPhi_calo;
   vector<int>     *tightJetIndex_JPT;
   vector<int>     *looseJetIndex_JPT;
   vector<float>   *loosejetPt_JPT;
   vector<float>   *loosejetEt_JPT;
   vector<float>   *loosejetEta_JPT;
   vector<float>   *loosejetPhi_JPT;
   vector<bool>    *loosejetPassLooseID_JPT;
   vector<bool>    *loosejetPassTightID_JPT;
   vector<float>   *loosejetEnergyFracHadronic_JPT;
   vector<int>     *loosejetFlavor_JPT;
   vector<int>     *loosejetGenPt_JPT;
   vector<int>     *loosejetGenEta_JPT;
   vector<int>     *loosejetGenPhi_JPT;
   vector<int>     *loosejetGenParticlePDGId_JPT;
   vector<float>   *loosejetInvisibleEnergy_JPT;
   vector<int>     *loosejetNTracks_JPT;
   vector<int>     *loosejetNSV_JPT;
   vector<float>   *loosejetSVUnWeightedMass_JPT;
   vector<float>   *loosejetSVWeightedMass_JPT;
   vector<float>   *loosejetSVUnWeightedLifetime_JPT;
   vector<float>   *loosejetSVWeightedLifetime_JPT;
   vector<float>   *loosejetSVUnWeightedCosTheta_JPT;
   vector<float>   *loosejetSVWeightedCosTheta_JPT;
   vector<float>   *loosejetBTagDisc_trackCountingHighPurBJetTags_JPT;
   vector<float>   *loosejetBTagDisc_trackCountingHighEffBJetTags_JPT;
   vector<float>   *loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags_JPT;
   vector<float>   *loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags_JPT;
   vector<float>   *loosejetBTagDisc_simpleSecondaryVertexBJetTags_JPT;
   vector<float>   *badjetPt_JPT;
   vector<float>   *badjetEta_JPT;
   vector<float>   *badjetPhi_JPT;
   vector<int>     *tightJetIndex_PF;
   vector<int>     *looseJetIndex_PF;
   vector<float>   *loosejetPt_PF;
   vector<float>   *loosejetEt_PF;
   vector<float>   *loosejetEta_PF;
   vector<float>   *loosejetPhi_PF;
   vector<bool>    *loosejetPassLooseID_PF;
   vector<bool>    *loosejetPassTightID_PF;
   vector<float>   *loosejetEnergyFracHadronic_PF;
   vector<int>     *loosejetFlavor_PF;
   vector<int>     *loosejetGenPt_PF;
   vector<int>     *loosejetGenEta_PF;
   vector<int>     *loosejetGenPhi_PF;
   vector<int>     *loosejetGenParticlePDGId_PF;
   vector<float>   *loosejetInvisibleEnergy_PF;
   vector<int>     *loosejetNTracks_PF;
   vector<int>     *loosejetNSV_PF;
   vector<float>   *loosejetSVUnWeightedMass_PF;
   vector<float>   *loosejetSVWeightedMass_PF;
   vector<float>   *loosejetSVUnWeightedLifetime_PF;
   vector<float>   *loosejetSVWeightedLifetime_PF;
   vector<float>   *loosejetSVUnWeightedCosTheta_PF;
   vector<float>   *loosejetSVWeightedCosTheta_PF;
   vector<float>   *loosejetBTagDisc_trackCountingHighPurBJetTags_PF;
   vector<float>   *loosejetBTagDisc_trackCountingHighEffBJetTags_PF;
   vector<float>   *loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags_PF;
   vector<float>   *loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags_PF;
   vector<float>   *loosejetBTagDisc_simpleSecondaryVertexBJetTags_PF;
   vector<float>   *badjetPt_PF;
   vector<float>   *badjetEta_PF;
   vector<float>   *badjetPhi_PF;
   vector<float>   *veryloosejetPtUncorr;
   vector<float>   *veryloosejetEtaUncorr;
   vector<float>   *veryloosejetPhiUncorr;
   Int_t           nbSSVM;
   Float_t         HT;
   Float_t         MHT;
   Float_t         MHTphi;
   Float_t         DeltaPhi_JetMHT1;
   Float_t         DeltaPhi_JetMHT2;
   Float_t         DeltaPhi_JetMHT3;
   Float_t         caloMET;
   Float_t         caloMETphi;
   Float_t         caloMETsig;
   Float_t         caloGenMET;
   Float_t         caloGenMETphi;
   Float_t         tcMET;
   Float_t         tcMETphi;
   Float_t         tcMETsig;
   Float_t         tcGenMET;
   Float_t         tcGenMETphi;
   Float_t         pfMET;
   Float_t         pfMETphi;
   Float_t         pfMETsig;
   Float_t         pfGenMET;
   Float_t         pfGenMETphi;
   vector<float>   *trackPt;
   vector<float>   *trackEta;
   vector<float>   *trackPhi;
   vector<float>   *tauPt;
   vector<float>   *tauEta;
   vector<float>   *tauPhi;
   vector<float>   *tauTaNC;
   vector<bool>    *tauID_againstElectron;
   vector<bool>    *tauID_againstMuon;
   vector<bool>    *tauID_byIsolation;
   vector<bool>    *tauID_byTaNCfrOnePercent;
   vector<bool>    *tauID_byTaNCfrHalfPercent;
   vector<bool>    *tauID_byTaNCfrQuarterPercent;
   vector<float>   *tauPt_PF;
   vector<float>   *tauEta_PF;
   vector<float>   *tauPhi_PF;
   vector<float>   *tauTaNC_PF;
   vector<bool>    *tauID_againstElectron_PF;
   vector<bool>    *tauID_againstMuon_PF;
   vector<bool>    *tauID_byIsolation_PF;
   vector<bool>    *tauID_byTaNCfrOnePercent_PF;
   vector<bool>    *tauID_byTaNCfrHalfPercent_PF;
   vector<bool>    *tauID_byTaNCfrQuarterPercent_PF;
   vector<bool>    *muonIsGlobalMuonPromptTight;
   vector<float>   *muonPt;
   vector<float>   *muonEta;
   vector<float>   *muonPhi;
   vector<float>   *muonTrackIso;
   vector<float>   *muonEcalIso;
   vector<float>   *muonHcalIso;
   vector<float>   *muonChi2;
   vector<float>   *muonNdof;
   vector<float>   *muonNhits;
   vector<float>   *muonTrackd0;
   vector<float>   *muonTrackPhi;
   vector<bool>    *muonPassID;
   vector<float>   *muonVtx_z;
   vector<float>   *muonEcalVeto;
   vector<float>   *muonHcalVeto;
   Int_t           nMuons;
   vector<float>   *eleEt;
   vector<float>   *eleEta;
   vector<float>   *elePhi;
   vector<float>   *eleTrackIso;
   vector<float>   *eleEcalIso;
   vector<float>   *eleHcalIso;
   vector<float>   *eledB;
   vector<float>   *eleVtx_z;
   vector<float>   *eleIDLoose;
   vector<float>   *eleIDRobustTight;
   vector<bool>    *elePassID;
   Int_t           nElectrons;
   vector<bool>    *muonIsGlobalMuonPromptTight_PF;
   vector<float>   *muonPt_PF;
   vector<float>   *muonEta_PF;
   vector<float>   *muonPhi_PF;
   vector<float>   *muonTrackIso_PF;
   vector<float>   *muonEcalIso_PF;
   vector<float>   *muonHcalIso_PF;
   vector<float>   *muonChi2_PF;
   vector<float>   *muonNdof_PF;
   vector<float>   *muonNhits_PF;
   vector<float>   *muonTrackd0_PF;
   vector<float>   *muonTrackPhi_PF;
   vector<bool>    *muonPassID_PF;
   vector<float>   *muonVtx_z_PF;
   vector<float>   *muonEcalVeto_PF;
   vector<float>   *muonHcalVeto_PF;
   Int_t           nMuons_PF;
   vector<float>   *eleEt_PF;
   vector<float>   *eleEta_PF;
   vector<float>   *elePhi_PF;
   vector<float>   *eleTrackIso_PF;
   vector<float>   *eleEcalIso_PF;
   vector<float>   *eleHcalIso_PF;
   vector<float>   *eledB_PF;
   vector<float>   *eleVtx_z_PF;
   vector<float>   *eleIDLoose_PF;
   vector<float>   *eleIDRobustTight_PF;
   vector<bool>    *elePassID_PF;
   Int_t           nElectrons_PF;
   Int_t           SUSY_nb;
   Float_t         qScale;
   vector<int>     *topDecayCode;
   Int_t           flavorHistory;

   // List of branches
   TBranch        *b_runNumber;   //!
   TBranch        *b_lumiSection;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_bsx;   //!
   TBranch        *b_bsy;   //!
   TBranch        *b_bsz;   //!
   TBranch        *b_cutResults;   //!
   TBranch        *b_passTrigger;   //!
   TBranch        *b_hltPrescale;   //!
   TBranch        *b_SUSYtriggerIndex;   //!
   TBranch        *b_pv_isFake;   //!
   TBranch        *b_pv_z;   //!
   TBranch        *b_pv_rho;   //!
   TBranch        *b_pv_chi2;   //!
   TBranch        *b_pv_ndof;   //!
   TBranch        *b_tightJetIndex_calo;   //!
   TBranch        *b_looseJetIndex_calo;   //!
   TBranch        *b_loosejetPt_calo;   //!
   TBranch        *b_loosejetEt_calo;   //!
   TBranch        *b_loosejetEta_calo;   //!
   TBranch        *b_loosejetPhi_calo;   //!
   TBranch        *b_loosejetPassLooseID_calo;   //!
   TBranch        *b_loosejetPassTightID_calo;   //!
   TBranch        *b_loosejetEnergyFracHadronic_calo;   //!
   TBranch        *b_loosejetFlavor_calo;   //!
   TBranch        *b_loosejetGenPt_calo;   //!
   TBranch        *b_loosejetGenEta_calo;   //!
   TBranch        *b_loosejetGenPhi_calo;   //!
   TBranch        *b_loosejetGenParticlePDGId_calo;   //!
   TBranch        *b_loosejetInvisibleEnergy_calo;   //!
   TBranch        *b_loosejetNTracks_calo;   //!
   TBranch        *b_loosejetNSV_calo;   //!
   TBranch        *b_loosejetSVUnWeightedMass_calo;   //!
   TBranch        *b_loosejetSVWeightedMass_calo;   //!
   TBranch        *b_loosejetSVUnWeightedLifetime_calo;   //!
   TBranch        *b_loosejetSVWeightedLifetime_calo;   //!
   TBranch        *b_loosejetSVUnWeightedCosTheta_calo;   //!
   TBranch        *b_loosejetSVWeightedCosTheta_calo;   //!
   TBranch        *b_loosejetBTagDisc_trackCountingHighPurBJetTags_calo;   //!
   TBranch        *b_loosejetBTagDisc_trackCountingHighEffBJetTags_calo;   //!
   TBranch        *b_loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags_calo;   //!
   TBranch        *b_loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags_calo;   //!
   TBranch        *b_loosejetBTagDisc_simpleSecondaryVertexBJetTags_calo;   //!
   TBranch        *b_badjetPt_calo;   //!
   TBranch        *b_badjetEta_calo;   //!
   TBranch        *b_badjetPhi_calo;   //!
   TBranch        *b_tightJetIndex_JPT;   //!
   TBranch        *b_looseJetIndex_JPT;   //!
   TBranch        *b_loosejetPt_JPT;   //!
   TBranch        *b_loosejetEt_JPT;   //!
   TBranch        *b_loosejetEta_JPT;   //!
   TBranch        *b_loosejetPhi_JPT;   //!
   TBranch        *b_loosejetPassLooseID_JPT;   //!
   TBranch        *b_loosejetPassTightID_JPT;   //!
   TBranch        *b_loosejetEnergyFracHadronic_JPT;   //!
   TBranch        *b_loosejetFlavor_JPT;   //!
   TBranch        *b_loosejetGenPt_JPT;   //!
   TBranch        *b_loosejetGenEta_JPT;   //!
   TBranch        *b_loosejetGenPhi_JPT;   //!
   TBranch        *b_loosejetGenParticlePDGId_JPT;   //!
   TBranch        *b_loosejetInvisibleEnergy_JPT;   //!
   TBranch        *b_loosejetNTracks_JPT;   //!
   TBranch        *b_loosejetNSV_JPT;   //!
   TBranch        *b_loosejetSVUnWeightedMass_JPT;   //!
   TBranch        *b_loosejetSVWeightedMass_JPT;   //!
   TBranch        *b_loosejetSVUnWeightedLifetime_JPT;   //!
   TBranch        *b_loosejetSVWeightedLifetime_JPT;   //!
   TBranch        *b_loosejetSVUnWeightedCosTheta_JPT;   //!
   TBranch        *b_loosejetSVWeightedCosTheta_JPT;   //!
   TBranch        *b_loosejetBTagDisc_trackCountingHighPurBJetTags_JPT;   //!
   TBranch        *b_loosejetBTagDisc_trackCountingHighEffBJetTags_JPT;   //!
   TBranch        *b_loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags_JPT;   //!
   TBranch        *b_loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags_JPT;   //!
   TBranch        *b_loosejetBTagDisc_simpleSecondaryVertexBJetTags_JPT;   //!
   TBranch        *b_badjetPt_JPT;   //!
   TBranch        *b_badjetEta_JPT;   //!
   TBranch        *b_badjetPhi_JPT;   //!
   TBranch        *b_tightJetIndex_PF;   //!
   TBranch        *b_looseJetIndex_PF;   //!
   TBranch        *b_loosejetPt_PF;   //!
   TBranch        *b_loosejetEt_PF;   //!
   TBranch        *b_loosejetEta_PF;   //!
   TBranch        *b_loosejetPhi_PF;   //!
   TBranch        *b_loosejetPassLooseID_PF;   //!
   TBranch        *b_loosejetPassTightID_PF;   //!
   TBranch        *b_loosejetEnergyFracHadronic_PF;   //!
   TBranch        *b_loosejetFlavor_PF;   //!
   TBranch        *b_loosejetGenPt_PF;   //!
   TBranch        *b_loosejetGenEta_PF;   //!
   TBranch        *b_loosejetGenPhi_PF;   //!
   TBranch        *b_loosejetGenParticlePDGId_PF;   //!
   TBranch        *b_loosejetInvisibleEnergy_PF;   //!
   TBranch        *b_loosejetNTracks_PF;   //!
   TBranch        *b_loosejetNSV_PF;   //!
   TBranch        *b_loosejetSVUnWeightedMass_PF;   //!
   TBranch        *b_loosejetSVWeightedMass_PF;   //!
   TBranch        *b_loosejetSVUnWeightedLifetime_PF;   //!
   TBranch        *b_loosejetSVWeightedLifetime_PF;   //!
   TBranch        *b_loosejetSVUnWeightedCosTheta_PF;   //!
   TBranch        *b_loosejetSVWeightedCosTheta_PF;   //!
   TBranch        *b_loosejetBTagDisc_trackCountingHighPurBJetTags_PF;   //!
   TBranch        *b_loosejetBTagDisc_trackCountingHighEffBJetTags_PF;   //!
   TBranch        *b_loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags_PF;   //!
   TBranch        *b_loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags_PF;   //!
   TBranch        *b_loosejetBTagDisc_simpleSecondaryVertexBJetTags_PF;   //!
   TBranch        *b_badjetPt_PF;   //!
   TBranch        *b_badjetEta_PF;   //!
   TBranch        *b_badjetPhi_PF;   //!
   TBranch        *b_veryloosejetPtUncorr;   //!
   TBranch        *b_veryloosejetEtaUncorr;   //!
   TBranch        *b_veryloosejetPhiUncorr;   //!
   TBranch        *b_nbSSVM;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_MHT;   //!
   TBranch        *b_MHTphi;   //!
   TBranch        *b_DeltaPhi_JetMHT1;   //!
   TBranch        *b_DeltaPhi_JetMHT2;   //!
   TBranch        *b_DeltaPhi_JetMHT3;   //!
   TBranch        *b_caloMET;   //!
   TBranch        *b_caloMETphi;   //!
   TBranch        *b_caloMETsig;   //!
   TBranch        *b_caloGenMET;   //!
   TBranch        *b_caloGenMETphi;   //!
   TBranch        *b_tcMET;   //!
   TBranch        *b_tcMETphi;   //!
   TBranch        *b_tcMETsig;   //!
   TBranch        *b_tcGenMET;   //!
   TBranch        *b_tcGenMETphi;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETphi;   //!
   TBranch        *b_pfMETsig;   //!
   TBranch        *b_pfGenMET;   //!
   TBranch        *b_pfGenMETphi;   //!
   TBranch        *b_trackPt;   //!
   TBranch        *b_trackEta;   //!
   TBranch        *b_trackPhi;   //!
   TBranch        *b_tauPt;   //!
   TBranch        *b_tauEta;   //!
   TBranch        *b_tauPhi;   //!
   TBranch        *b_tauTaNC;   //!
   TBranch        *b_tauID_againstElectron;   //!
   TBranch        *b_tauID_againstMuon;   //!
   TBranch        *b_tauID_byIsolation;   //!
   TBranch        *b_tauID_byTaNCfrOnePercent;   //!
   TBranch        *b_tauID_byTaNCfrHalfPercent;   //!
   TBranch        *b_tauID_byTaNCfrQuarterPercent;   //!
   TBranch        *b_tauPt_PF;   //!
   TBranch        *b_tauEta_PF;   //!
   TBranch        *b_tauPhi_PF;   //!
   TBranch        *b_tauTaNC_PF;   //!
   TBranch        *b_tauID_againstElectron_PF;   //!
   TBranch        *b_tauID_againstMuon_PF;   //!
   TBranch        *b_tauID_byIsolation_PF;   //!
   TBranch        *b_tauID_byTaNCfrOnePercent_PF;   //!
   TBranch        *b_tauID_byTaNCfrHalfPercent_PF;   //!
   TBranch        *b_tauID_byTaNCfrQuarterPercent_PF;   //!
   TBranch        *b_muonIsGlobalMuonPromptTight;   //!
   TBranch        *b_muonPt;   //!
   TBranch        *b_muonEta;   //!
   TBranch        *b_muonPhi;   //!
   TBranch        *b_muonTrackIso;   //!
   TBranch        *b_muonEcalIso;   //!
   TBranch        *b_muonHcalIso;   //!
   TBranch        *b_muonChi2;   //!
   TBranch        *b_muonNdof;   //!
   TBranch        *b_muonNhits;   //!
   TBranch        *b_muonTrackd0;   //!
   TBranch        *b_muonTrackPhi;   //!
   TBranch        *b_muonPassID;   //!
   TBranch        *b_muonVtx_z;   //!
   TBranch        *b_muonEcalVeto;   //!
   TBranch        *b_muonHcalVeto;   //!
   TBranch        *b_nMuons;   //!
   TBranch        *b_eleEt;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleTrackIso;   //!
   TBranch        *b_eleEcalIso;   //!
   TBranch        *b_eleHcalIso;   //!
   TBranch        *b_eledB;   //!
   TBranch        *b_eleVtx_z;   //!
   TBranch        *b_eleIDLoose;   //!
   TBranch        *b_eleIDRobustTight;   //!
   TBranch        *b_elePassID;   //!
   TBranch        *b_nElectrons;   //!
   TBranch        *b_muonIsGlobalMuonPromptTight_PF;   //!
   TBranch        *b_muonPt_PF;   //!
   TBranch        *b_muonEta_PF;   //!
   TBranch        *b_muonPhi_PF;   //!
   TBranch        *b_muonTrackIso_PF;   //!
   TBranch        *b_muonEcalIso_PF;   //!
   TBranch        *b_muonHcalIso_PF;   //!
   TBranch        *b_muonChi2_PF;   //!
   TBranch        *b_muonNdof_PF;   //!
   TBranch        *b_muonNhits_PF;   //!
   TBranch        *b_muonTrackd0_PF;   //!
   TBranch        *b_muonTrackPhi_PF;   //!
   TBranch        *b_muonPassID_PF;   //!
   TBranch        *b_muonVtx_z_PF;   //!
   TBranch        *b_muonEcalVeto_PF;   //!
   TBranch        *b_muonHcalVeto_PF;   //!
   TBranch        *b_nMuons_PF;   //!
   TBranch        *b_eleEt_PF;   //!
   TBranch        *b_eleEta_PF;   //!
   TBranch        *b_elePhi_PF;   //!
   TBranch        *b_eleTrackIso_PF;   //!
   TBranch        *b_eleEcalIso_PF;   //!
   TBranch        *b_eleHcalIso_PF;   //!
   TBranch        *b_eledB_PF;   //!
   TBranch        *b_eleVtx_z_PF;   //!
   TBranch        *b_eleIDLoose_PF;   //!
   TBranch        *b_eleIDRobustTight_PF;   //!
   TBranch        *b_elePassID_PF;   //!
   TBranch        *b_nElectrons_PF;   //!
   TBranch        *b_SUSY_nb;   //!
   TBranch        *b_qScale;   //!
   TBranch        *b_topDecayCode;   //!
   TBranch        *b_flavorHistory;   //!

   basicLoop(TTree *tree=0, TTree *infotree=0);    // ========================================== begin, end
   virtual ~basicLoop();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   // ========================================== begin
   virtual void     Loop(unsigned int dataindex=0);
   virtual void     exampleLoop();
   virtual void     screendump();
   virtual void     nbLoop();
   virtual void     ABCDtree(unsigned int dataindex=0);
   virtual void triggerPlot();
   virtual void triggerPlotData();
   virtual void triggerTest();
   virtual void Nminus1plots();

   bool isV00_01_02();
   bool isV00_01_03();
   bool isV00_01_05();

   void fillTightJetInfo();
   void InitJets();

   void printState();

   //   int getCutFlow(TString cut);
   bool setCutScheme(CutScheme cutscheme);
   void setMETType(METType mettype);
   void setJetType(jetType jettype);
   void setLeptonType(leptonType leptontype);
   void setMETRange(METRange metrange);
   void setDPType(dpType dptype);
   void setIgnoredCut(const TString cutTag);
   void resetIgnoredCut() ;
   void setBCut(unsigned int nb);
   void setEleReq(int ne); //change Cut suffix to Req because it is an == requirement, not >=
   void setMuonReq(int nmu);
   TString getCutDescriptionString();
   TString getBCutDescriptionString();

   void specifyEvent(ULong64_t run, ULong64_t lumisection, ULong64_t event);
   bool eventIsSpecified();

   float getMET(); //return MET determined by theMETType_
   float getMETphi(); //return MET determined by theMETType_
   
   void cutflow();
   void cutflowPlotter();
   bool cutRequired(TString cutTag) ;
   bool passCut(TString cutTag) ;
   TString getSampleName(TString inname="") ;
   double getCrossSection(TString inname="") ;
   TString findInputName() ;

   bool passMuVetoRA2() ;
   bool passEleVetoRA2() ;
   bool passMuVetoSync1() ;
   bool passEleVetoSync1() ;

   int countMuSync1() ;
   int countEleSync1() ;

   bool passTauVeto();

   bool passSSVM(int i);

   bool isGoodJet_Sync1(unsigned int ijet);
   bool isGoodJet(unsigned int ijet); //index here is on the loose jet list
   bool isGoodJet30(unsigned int ijet); //index here is on the loose jet list
   bool isLooseJet_Sync1(unsigned int ijet); //looser pt cut
   unsigned int nGoodJets_Sync1();
   unsigned int nGoodJets();
   float jetPtOfN(unsigned int n);
   float jetPhiOfN(unsigned int n);
   float jetEtaOfN(unsigned int n);
   int countBJets_Sync1();
   int countBJets();
   bool passDeltaPhi_Sync1() ;
   float getHT_Sync1();
   float getHT();
   float getMHT();
   float getMHTphi();

   bool passPV();

   bool passHLT(); //bool printedHLT_; 
   TString lastTriggerPass_; //not really filled...but keep it in so that basicLoop.C compiles
   bool passHLT(const TString & triggerName);
   int getHLTPrescale(const TString & triggerName);

   double getDeltaPhiMPTMET();
   double getMinDeltaPhibMET() ;
   double getMinDeltaPhiMET(unsigned int maxjets) ;
   //   double getMinDeltaPhiMHT(unsigned int maxjets) ; //deprecated because MHT is now just another type of MET
   double getDeltaPhib1b2();
   double getUncorrectedHT(const double threshold);
   int getTopDecayCategory();
   double getMinDeltaPhi_bj(unsigned int bindex);
   double getOverallMinDeltaR_bj();
   double getMinDeltaR_bj(unsigned int bindex);
   bool passBCut( unsigned int bindex);
   double getDeltaPhi(double phi1, double phi2);
   // ========================================== end

};

#endif

#ifdef basicLoop_cxx
//====================== begin
basicLoop::basicLoop(TTree *tree, TTree *infotree)
  :  theCutScheme_(kRA2),
     theMETType_(kMET),
     theMETRange_(kHigh),
     theJetType_(kCalo),
     theLeptonType_(kNormal),
     theDPType_(kminDP),
     nBcut_(0),
     ne_(0),
     nmu_(0),
     isData_(false),
     //     printedHLT_(false),
     lastTriggerPass_("")
//====================== end

{
   // ========================================== begin
  //don't do the default root thing of loading a default tree
  if (tree == 0 || infotree==0) {
    cout<<"One of the required input trees has not been provided as an argument to the constructor of basicLoop!"<<endl;
    cout<<"It is required to provide pointers to both the tree and the infotree!"<<endl;
    cout<<"basicLoop::basicLoop(TTree *tree, TTree *infotree)\n---------------------------------"<<endl;
    cout<<"Careful! The class has not been constructed properly!"<<endl;
    return;
   }
  // ========================================== end

   Init(tree);

   // ========================================== begin
   specifiedEvents_.clear();

   triggerList_.clear();
   triggerListAll_.clear();
   if (infotree!=0) {
     std::set<TString> triggersForCut;
     triggersForCut.insert("HLT_HT100U");
     triggersForCut.insert("HLT_HT120U");
     triggersForCut.insert("HLT_HT140U");
     triggersForCut.insert("HLT_HT150U");
     triggersForCut.insert("HLT_HT150U_v3");
     triggersForCut.insert("HLT_HT200");

     Long64_t ninfo = infotree->GetEntries();
     if (ninfo > 0) {
       std::vector<std::string> * triggerList=0;
       infotree->SetBranchAddress("triggerList", &triggerList);
       infotree->GetEntry(0);
       for (unsigned int itrig=0; itrig<triggerList->size(); itrig++) {
	 //	 std::cout<<"in infotree, found trigger: "<<triggerList->at(itrig)<<std::endl;
	 if ( triggersForCut.find( TString(triggerList->at(itrig))) != triggersForCut.end()) {
	   //	   std::cout<<"Will accept events with trigger: "<<triggerList->at(itrig)<<std::endl;
	   triggerList_[  TString(triggerList->at(itrig))] = itrig;
	 }
	 triggerListAll_[ TString(triggerList->at(itrig))] = itrig;
       }
     }
   }
   //this routine initializes the various cutNames_ etc vectors
   setCutScheme(kRA2);   
   setJetType(kCalo);
   //set isData_
   if (  getSampleName(findInputName()) == "data") {cout<<"Sample is real data!"<<endl; isData_=true;}
   // ========================================== end

}

basicLoop::~basicLoop()
{
  // ========================================== begin
  //No! stupid ROOT! this leads to a double deletion!
  //if (!fChain) return;
  // delete fChain->GetCurrentFile();
  // ========================================== end
}

Int_t basicLoop::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
// ========================================== begin
   Int_t n=fChain->GetEntry(entry);

   //this order is critical!
   InitJets();
   fillTightJetInfo();

   return n;
  // ========================================== end
}
Long64_t basicLoop::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void basicLoop::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   cutResults = 0;
   passTrigger = 0;
   hltPrescale = 0;
   pv_isFake = 0;
   pv_z = 0;
   pv_rho = 0;
   pv_chi2 = 0;
   pv_ndof = 0;
   tightJetIndex_calo = 0;
   looseJetIndex_calo = 0;
   loosejetPt_calo = 0;
   loosejetEt_calo = 0;
   loosejetEta_calo = 0;
   loosejetPhi_calo = 0;
   loosejetPassLooseID_calo = 0;
   loosejetPassTightID_calo = 0;
   loosejetEnergyFracHadronic_calo = 0;
   loosejetFlavor_calo = 0;
   loosejetGenPt_calo = 0;
   loosejetGenEta_calo = 0;
   loosejetGenPhi_calo = 0;
   loosejetGenParticlePDGId_calo = 0;
   loosejetInvisibleEnergy_calo = 0;
   loosejetNTracks_calo = 0;
   loosejetNSV_calo = 0;
   loosejetSVUnWeightedMass_calo = 0;
   loosejetSVWeightedMass_calo = 0;
   loosejetSVUnWeightedLifetime_calo = 0;
   loosejetSVWeightedLifetime_calo = 0;
   loosejetSVUnWeightedCosTheta_calo = 0;
   loosejetSVWeightedCosTheta_calo = 0;
   loosejetBTagDisc_trackCountingHighPurBJetTags_calo = 0;
   loosejetBTagDisc_trackCountingHighEffBJetTags_calo = 0;
   loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags_calo = 0;
   loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags_calo = 0;
   loosejetBTagDisc_simpleSecondaryVertexBJetTags_calo = 0;
   badjetPt_calo = 0;
   badjetEta_calo = 0;
   badjetPhi_calo = 0;
   tightJetIndex_JPT = 0;
   looseJetIndex_JPT = 0;
   loosejetPt_JPT = 0;
   loosejetEt_JPT = 0;
   loosejetEta_JPT = 0;
   loosejetPhi_JPT = 0;
   loosejetPassLooseID_JPT = 0;
   loosejetPassTightID_JPT = 0;
   loosejetEnergyFracHadronic_JPT = 0;
   loosejetFlavor_JPT = 0;
   loosejetGenPt_JPT = 0;
   loosejetGenEta_JPT = 0;
   loosejetGenPhi_JPT = 0;
   loosejetGenParticlePDGId_JPT = 0;
   loosejetInvisibleEnergy_JPT = 0;
   loosejetNTracks_JPT = 0;
   loosejetNSV_JPT = 0;
   loosejetSVUnWeightedMass_JPT = 0;
   loosejetSVWeightedMass_JPT = 0;
   loosejetSVUnWeightedLifetime_JPT = 0;
   loosejetSVWeightedLifetime_JPT = 0;
   loosejetSVUnWeightedCosTheta_JPT = 0;
   loosejetSVWeightedCosTheta_JPT = 0;
   loosejetBTagDisc_trackCountingHighPurBJetTags_JPT = 0;
   loosejetBTagDisc_trackCountingHighEffBJetTags_JPT = 0;
   loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags_JPT = 0;
   loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags_JPT = 0;
   loosejetBTagDisc_simpleSecondaryVertexBJetTags_JPT = 0;
   badjetPt_JPT = 0;
   badjetEta_JPT = 0;
   badjetPhi_JPT = 0;
   tightJetIndex_PF = 0;
   looseJetIndex_PF = 0;
   loosejetPt_PF = 0;
   loosejetEt_PF = 0;
   loosejetEta_PF = 0;
   loosejetPhi_PF = 0;
   loosejetPassLooseID_PF = 0;
   loosejetPassTightID_PF = 0;
   loosejetEnergyFracHadronic_PF = 0;
   loosejetFlavor_PF = 0;
   loosejetGenPt_PF = 0;
   loosejetGenEta_PF = 0;
   loosejetGenPhi_PF = 0;
   loosejetGenParticlePDGId_PF = 0;
   loosejetInvisibleEnergy_PF = 0;
   loosejetNTracks_PF = 0;
   loosejetNSV_PF = 0;
   loosejetSVUnWeightedMass_PF = 0;
   loosejetSVWeightedMass_PF = 0;
   loosejetSVUnWeightedLifetime_PF = 0;
   loosejetSVWeightedLifetime_PF = 0;
   loosejetSVUnWeightedCosTheta_PF = 0;
   loosejetSVWeightedCosTheta_PF = 0;
   loosejetBTagDisc_trackCountingHighPurBJetTags_PF = 0;
   loosejetBTagDisc_trackCountingHighEffBJetTags_PF = 0;
   loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags_PF = 0;
   loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags_PF = 0;
   loosejetBTagDisc_simpleSecondaryVertexBJetTags_PF = 0;
   badjetPt_PF = 0;
   badjetEta_PF = 0;
   badjetPhi_PF = 0;
   veryloosejetPtUncorr = 0;
   veryloosejetEtaUncorr = 0;
   veryloosejetPhiUncorr = 0;
   trackPt = 0;
   trackEta = 0;
   trackPhi = 0;
   tauPt = 0;
   tauEta = 0;
   tauPhi = 0;
   tauTaNC = 0;
   tauID_againstElectron = 0;
   tauID_againstMuon = 0;
   tauID_byIsolation = 0;
   tauID_byTaNCfrOnePercent = 0;
   tauID_byTaNCfrHalfPercent = 0;
   tauID_byTaNCfrQuarterPercent = 0;
   tauPt_PF = 0;
   tauEta_PF = 0;
   tauPhi_PF = 0;
   tauTaNC_PF = 0;
   tauID_againstElectron_PF = 0;
   tauID_againstMuon_PF = 0;
   tauID_byIsolation_PF = 0;
   tauID_byTaNCfrOnePercent_PF = 0;
   tauID_byTaNCfrHalfPercent_PF = 0;
   tauID_byTaNCfrQuarterPercent_PF = 0;
   muonIsGlobalMuonPromptTight = 0;
   muonPt = 0;
   muonEta = 0;
   muonPhi = 0;
   muonTrackIso = 0;
   muonEcalIso = 0;
   muonHcalIso = 0;
   muonChi2 = 0;
   muonNdof = 0;
   muonNhits = 0;
   muonTrackd0 = 0;
   muonTrackPhi = 0;
   muonPassID = 0;
   muonVtx_z = 0;
   muonEcalVeto = 0;
   muonHcalVeto = 0;
   eleEt = 0;
   eleEta = 0;
   elePhi = 0;
   eleTrackIso = 0;
   eleEcalIso = 0;
   eleHcalIso = 0;
   eledB = 0;
   eleVtx_z = 0;
   eleIDLoose = 0;
   eleIDRobustTight = 0;
   elePassID = 0;
   muonIsGlobalMuonPromptTight_PF = 0;
   muonPt_PF = 0;
   muonEta_PF = 0;
   muonPhi_PF = 0;
   muonTrackIso_PF = 0;
   muonEcalIso_PF = 0;
   muonHcalIso_PF = 0;
   muonChi2_PF = 0;
   muonNdof_PF = 0;
   muonNhits_PF = 0;
   muonTrackd0_PF = 0;
   muonTrackPhi_PF = 0;
   muonPassID_PF = 0;
   muonVtx_z_PF = 0;
   muonEcalVeto_PF = 0;
   muonHcalVeto_PF = 0;
   eleEt_PF = 0;
   eleEta_PF = 0;
   elePhi_PF = 0;
   eleTrackIso_PF = 0;
   eleEcalIso_PF = 0;
   eleHcalIso_PF = 0;
   eledB_PF = 0;
   eleVtx_z_PF = 0;
   eleIDLoose_PF = 0;
   eleIDRobustTight_PF = 0;
   elePassID_PF = 0;
   topDecayCode = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("lumiSection", &lumiSection, &b_lumiSection);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("bsx", &bsx, &b_bsx);
   fChain->SetBranchAddress("bsy", &bsy, &b_bsy);
   fChain->SetBranchAddress("bsz", &bsz, &b_bsz);
   fChain->SetBranchAddress("cutResults", &cutResults, &b_cutResults);
   fChain->SetBranchAddress("passTrigger", &passTrigger, &b_passTrigger);
   fChain->SetBranchAddress("hltPrescale", &hltPrescale, &b_hltPrescale);
   fChain->SetBranchAddress("SUSYtriggerIndex", &SUSYtriggerIndex, &b_SUSYtriggerIndex);
   fChain->SetBranchAddress("pv_isFake", &pv_isFake, &b_pv_isFake);
   fChain->SetBranchAddress("pv_z", &pv_z, &b_pv_z);
   fChain->SetBranchAddress("pv_rho", &pv_rho, &b_pv_rho);
   fChain->SetBranchAddress("pv_chi2", &pv_chi2, &b_pv_chi2);
   fChain->SetBranchAddress("pv_ndof", &pv_ndof, &b_pv_ndof);
   fChain->SetBranchAddress("tightJetIndex_calo", &tightJetIndex_calo, &b_tightJetIndex_calo);
   fChain->SetBranchAddress("looseJetIndex_calo", &looseJetIndex_calo, &b_looseJetIndex_calo);
   fChain->SetBranchAddress("loosejetPt_calo", &loosejetPt_calo, &b_loosejetPt_calo);
   fChain->SetBranchAddress("loosejetEt_calo", &loosejetEt_calo, &b_loosejetEt_calo);
   fChain->SetBranchAddress("loosejetEta_calo", &loosejetEta_calo, &b_loosejetEta_calo);
   fChain->SetBranchAddress("loosejetPhi_calo", &loosejetPhi_calo, &b_loosejetPhi_calo);
   fChain->SetBranchAddress("loosejetPassLooseID_calo", &loosejetPassLooseID_calo, &b_loosejetPassLooseID_calo);
   if (isV00_01_02()) fChain->SetBranchAddress("loosejetPassTightID_calo", &loosejetPassTightID_calo, &b_loosejetPassTightID_calo);
   fChain->SetBranchAddress("loosejetEnergyFracHadronic_calo", &loosejetEnergyFracHadronic_calo, &b_loosejetEnergyFracHadronic_calo);
   fChain->SetBranchAddress("loosejetFlavor_calo", &loosejetFlavor_calo, &b_loosejetFlavor_calo);
   if (isV00_01_02()) { //make backwards compatible
   fChain->SetBranchAddress("loosejetGenPt_calo", &loosejetGenPt_calo, &b_loosejetGenPt_calo);
   fChain->SetBranchAddress("loosejetGenEta_calo", &loosejetGenEta_calo, &b_loosejetGenEta_calo);
   fChain->SetBranchAddress("loosejetGenPhi_calo", &loosejetGenPhi_calo, &b_loosejetGenPhi_calo);
   fChain->SetBranchAddress("loosejetGenParticlePDGId_calo", &loosejetGenParticlePDGId_calo, &b_loosejetGenParticlePDGId_calo);
   }
   fChain->SetBranchAddress("loosejetInvisibleEnergy_calo", &loosejetInvisibleEnergy_calo, &b_loosejetInvisibleEnergy_calo);
   if (isV00_01_02()) { //make backwards compatible
   fChain->SetBranchAddress("loosejetNSV_calo", &loosejetNSV_calo, &b_loosejetNSV_calo);
   if (isV00_01_03())   fChain->SetBranchAddress("loosejetNTracks_calo", &loosejetNTracks_calo, &b_loosejetNTracks_calo);
   fChain->SetBranchAddress("loosejetSVUnWeightedMass_calo", &loosejetSVUnWeightedMass_calo, &b_loosejetSVUnWeightedMass_calo);
   fChain->SetBranchAddress("loosejetSVWeightedMass_calo", &loosejetSVWeightedMass_calo, &b_loosejetSVWeightedMass_calo);
   fChain->SetBranchAddress("loosejetSVUnWeightedLifetime_calo", &loosejetSVUnWeightedLifetime_calo, &b_loosejetSVUnWeightedLifetime_calo);
   fChain->SetBranchAddress("loosejetSVWeightedLifetime_calo", &loosejetSVWeightedLifetime_calo, &b_loosejetSVWeightedLifetime_calo);
   fChain->SetBranchAddress("loosejetSVUnWeightedCosTheta_calo", &loosejetSVUnWeightedCosTheta_calo, &b_loosejetSVUnWeightedCosTheta_calo);
   fChain->SetBranchAddress("loosejetSVWeightedCosTheta_calo", &loosejetSVWeightedCosTheta_calo, &b_loosejetSVWeightedCosTheta_calo);
   }
   fChain->SetBranchAddress("loosejetBTagDisc_trackCountingHighPurBJetTags_calo", &loosejetBTagDisc_trackCountingHighPurBJetTags_calo, &b_loosejetBTagDisc_trackCountingHighPurBJetTags_calo);
   fChain->SetBranchAddress("loosejetBTagDisc_trackCountingHighEffBJetTags_calo", &loosejetBTagDisc_trackCountingHighEffBJetTags_calo, &b_loosejetBTagDisc_trackCountingHighEffBJetTags_calo);
   fChain->SetBranchAddress("loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags_calo", &loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags_calo, &b_loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags_calo);
   fChain->SetBranchAddress("loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags_calo", &loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags_calo, &b_loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags_calo);
   fChain->SetBranchAddress("loosejetBTagDisc_simpleSecondaryVertexBJetTags_calo", &loosejetBTagDisc_simpleSecondaryVertexBJetTags_calo, &b_loosejetBTagDisc_simpleSecondaryVertexBJetTags_calo);
   fChain->SetBranchAddress("badjetPt_calo", &badjetPt_calo, &b_badjetPt_calo);
   fChain->SetBranchAddress("badjetEta_calo", &badjetEta_calo, &b_badjetEta_calo);
   fChain->SetBranchAddress("badjetPhi_calo", &badjetPhi_calo, &b_badjetPhi_calo);
   fChain->SetBranchAddress("tightJetIndex_PF", &tightJetIndex_PF, &b_tightJetIndex_PF);
   fChain->SetBranchAddress("looseJetIndex_PF", &looseJetIndex_PF, &b_looseJetIndex_PF);
   fChain->SetBranchAddress("loosejetPt_PF", &loosejetPt_PF, &b_loosejetPt_PF);
   fChain->SetBranchAddress("loosejetEt_PF", &loosejetEt_PF, &b_loosejetEt_PF);
   fChain->SetBranchAddress("loosejetEta_PF", &loosejetEta_PF, &b_loosejetEta_PF);
   fChain->SetBranchAddress("loosejetPhi_PF", &loosejetPhi_PF, &b_loosejetPhi_PF);
   fChain->SetBranchAddress("loosejetPassLooseID_PF", &loosejetPassLooseID_PF, &b_loosejetPassLooseID_PF);
   if (isV00_01_02())fChain->SetBranchAddress("loosejetPassTightID_PF", &loosejetPassTightID_PF, &b_loosejetPassTightID_PF);
   fChain->SetBranchAddress("loosejetEnergyFracHadronic_PF", &loosejetEnergyFracHadronic_PF, &b_loosejetEnergyFracHadronic_PF);
   fChain->SetBranchAddress("loosejetFlavor_PF", &loosejetFlavor_PF, &b_loosejetFlavor_PF);
   if (isV00_01_02()) {
   fChain->SetBranchAddress("loosejetGenPt_PF", &loosejetGenPt_PF, &b_loosejetGenPt_PF);
   fChain->SetBranchAddress("loosejetGenEta_PF", &loosejetGenEta_PF, &b_loosejetGenEta_PF);
   fChain->SetBranchAddress("loosejetGenPhi_PF", &loosejetGenPhi_PF, &b_loosejetGenPhi_PF);
   }
   fChain->SetBranchAddress("loosejetGenParticlePDGId_PF", &loosejetGenParticlePDGId_PF, &b_loosejetGenParticlePDGId_PF);
   fChain->SetBranchAddress("loosejetInvisibleEnergy_PF", &loosejetInvisibleEnergy_PF, &b_loosejetInvisibleEnergy_PF);
   if (isV00_01_02()) {
   fChain->SetBranchAddress("loosejetNSV_PF", &loosejetNSV_PF, &b_loosejetNSV_PF);
   if (isV00_01_03())   fChain->SetBranchAddress("loosejetNTracks_PF", &loosejetNTracks_PF, &b_loosejetNTracks_PF);
   fChain->SetBranchAddress("loosejetSVUnWeightedMass_PF", &loosejetSVUnWeightedMass_PF, &b_loosejetSVUnWeightedMass_PF);
   fChain->SetBranchAddress("loosejetSVWeightedMass_PF", &loosejetSVWeightedMass_PF, &b_loosejetSVWeightedMass_PF);
   fChain->SetBranchAddress("loosejetSVUnWeightedLifetime_PF", &loosejetSVUnWeightedLifetime_PF, &b_loosejetSVUnWeightedLifetime_PF);
   fChain->SetBranchAddress("loosejetSVWeightedLifetime_PF", &loosejetSVWeightedLifetime_PF, &b_loosejetSVWeightedLifetime_PF);
   fChain->SetBranchAddress("loosejetSVUnWeightedCosTheta_PF", &loosejetSVUnWeightedCosTheta_PF, &b_loosejetSVUnWeightedCosTheta_PF);
   fChain->SetBranchAddress("loosejetSVWeightedCosTheta_PF", &loosejetSVWeightedCosTheta_PF, &b_loosejetSVWeightedCosTheta_PF);
   }
   fChain->SetBranchAddress("loosejetBTagDisc_trackCountingHighPurBJetTags_PF", &loosejetBTagDisc_trackCountingHighPurBJetTags_PF, &b_loosejetBTagDisc_trackCountingHighPurBJetTags_PF);
   fChain->SetBranchAddress("loosejetBTagDisc_trackCountingHighEffBJetTags_PF", &loosejetBTagDisc_trackCountingHighEffBJetTags_PF, &b_loosejetBTagDisc_trackCountingHighEffBJetTags_PF);
   fChain->SetBranchAddress("loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags_PF", &loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags_PF, &b_loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags_PF);
   fChain->SetBranchAddress("loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags_PF", &loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags_PF, &b_loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags_PF);
   fChain->SetBranchAddress("loosejetBTagDisc_simpleSecondaryVertexBJetTags_PF", &loosejetBTagDisc_simpleSecondaryVertexBJetTags_PF, &b_loosejetBTagDisc_simpleSecondaryVertexBJetTags_PF);
   fChain->SetBranchAddress("badjetPt_PF", &badjetPt_PF, &b_badjetPt_PF);
   fChain->SetBranchAddress("badjetEta_PF", &badjetEta_PF, &b_badjetEta_PF);
   fChain->SetBranchAddress("badjetPhi_PF", &badjetPhi_PF, &b_badjetPhi_PF);
   fChain->SetBranchAddress("veryloosejetPtUncorr", &veryloosejetPtUncorr, &b_veryloosejetPtUncorr);
   fChain->SetBranchAddress("veryloosejetEtaUncorr", &veryloosejetEtaUncorr, &b_veryloosejetEtaUncorr);
   fChain->SetBranchAddress("veryloosejetPhiUncorr", &veryloosejetPhiUncorr, &b_veryloosejetPhiUncorr);
   fChain->SetBranchAddress("nbSSVM", &nbSSVM, &b_nbSSVM);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("MHT", &MHT, &b_MHT);
   fChain->SetBranchAddress("MHTphi", &MHTphi, &b_MHTphi);
   fChain->SetBranchAddress("DeltaPhi_JetMHT1", &DeltaPhi_JetMHT1, &b_DeltaPhi_JetMHT1);
   fChain->SetBranchAddress("DeltaPhi_JetMHT2", &DeltaPhi_JetMHT2, &b_DeltaPhi_JetMHT2);
   fChain->SetBranchAddress("DeltaPhi_JetMHT3", &DeltaPhi_JetMHT3, &b_DeltaPhi_JetMHT3);
   fChain->SetBranchAddress("caloMET", &caloMET, &b_caloMET);
   fChain->SetBranchAddress("caloMETphi", &caloMETphi, &b_caloMETphi);
   fChain->SetBranchAddress("caloMETsig", &caloMETsig, &b_caloMETsig);
   if (isV00_01_02())   fChain->SetBranchAddress("caloGenMET", &caloGenMET, &b_caloGenMET);
   if (isV00_01_02())   fChain->SetBranchAddress("caloGenMETphi", &caloGenMETphi, &b_caloGenMETphi);
   fChain->SetBranchAddress("tcMET", &tcMET, &b_tcMET);
   fChain->SetBranchAddress("tcMETphi", &tcMETphi, &b_tcMETphi);
   fChain->SetBranchAddress("tcMETsig", &tcMETsig, &b_tcMETsig);
   if (isV00_01_02())    fChain->SetBranchAddress("tcGenMET", &tcGenMET, &b_tcGenMET);
   if (isV00_01_02())    fChain->SetBranchAddress("tcGenMETphi", &tcGenMETphi, &b_tcGenMETphi);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMETphi", &pfMETphi, &b_pfMETphi);
   fChain->SetBranchAddress("pfMETsig", &pfMETsig, &b_pfMETsig);
   if (isV00_01_02())    fChain->SetBranchAddress("pfGenMET", &pfGenMET, &b_pfGenMET);
   if (isV00_01_02())    fChain->SetBranchAddress("pfGenMETphi", &pfGenMETphi, &b_pfGenMETphi);
   fChain->SetBranchAddress("trackPt", &trackPt, &b_trackPt);
   fChain->SetBranchAddress("trackEta", &trackEta, &b_trackEta);
   fChain->SetBranchAddress("trackPhi", &trackPhi, &b_trackPhi);
   fChain->SetBranchAddress("muonIsGlobalMuonPromptTight", &muonIsGlobalMuonPromptTight, &b_muonIsGlobalMuonPromptTight);
   fChain->SetBranchAddress("muonPt", &muonPt, &b_muonPt);
   fChain->SetBranchAddress("muonEta", &muonEta, &b_muonEta);
   fChain->SetBranchAddress("muonTrackIso", &muonTrackIso, &b_muonTrackIso);
   fChain->SetBranchAddress("muonEcalIso", &muonEcalIso, &b_muonEcalIso);
   fChain->SetBranchAddress("muonHcalIso", &muonHcalIso, &b_muonHcalIso);
   fChain->SetBranchAddress("muonChi2", &muonChi2, &b_muonChi2);
   fChain->SetBranchAddress("muonNdof", &muonNdof, &b_muonNdof);
   fChain->SetBranchAddress("muonNhits", &muonNhits, &b_muonNhits);
   fChain->SetBranchAddress("muonTrackd0", &muonTrackd0, &b_muonTrackd0);
   fChain->SetBranchAddress("muonTrackPhi", &muonTrackPhi, &b_muonTrackPhi);
   fChain->SetBranchAddress("muonPassID", &muonPassID, &b_muonPassID);
   if (isV00_01_02())    fChain->SetBranchAddress("muonVtx_z", &muonVtx_z, &b_muonVtx_z);
   fChain->SetBranchAddress("muonEcalVeto", &muonEcalVeto, &b_muonEcalVeto);
   fChain->SetBranchAddress("muonHcalVeto", &muonHcalVeto, &b_muonHcalVeto);
   fChain->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
   fChain->SetBranchAddress("eleEt", &eleEt, &b_eleEt);
   fChain->SetBranchAddress("eleEta", &eleEta, &b_eleEta);
   fChain->SetBranchAddress("eleTrackIso", &eleTrackIso, &b_eleTrackIso);
   fChain->SetBranchAddress("eleEcalIso", &eleEcalIso, &b_eleEcalIso);
   fChain->SetBranchAddress("eleHcalIso", &eleHcalIso, &b_eleHcalIso);
   if (isV00_01_02())    fChain->SetBranchAddress("eledB", &eledB, &b_eledB);
   if (isV00_01_02())    fChain->SetBranchAddress("eleVtx_z", &eleVtx_z, &b_eleVtx_z);
   fChain->SetBranchAddress("eleIDLoose", &eleIDLoose, &b_eleIDLoose);
   fChain->SetBranchAddress("eleIDRobustTight", &eleIDRobustTight, &b_eleIDRobustTight);
   fChain->SetBranchAddress("elePassID", &elePassID, &b_elePassID);
   fChain->SetBranchAddress("nElectrons", &nElectrons, &b_nElectrons);
   fChain->SetBranchAddress("muonIsGlobalMuonPromptTight_PF", &muonIsGlobalMuonPromptTight_PF, &b_muonIsGlobalMuonPromptTight_PF);
   fChain->SetBranchAddress("muonPt_PF", &muonPt_PF, &b_muonPt_PF);
   fChain->SetBranchAddress("muonEta_PF", &muonEta_PF, &b_muonEta_PF);
   fChain->SetBranchAddress("muonTrackIso_PF", &muonTrackIso_PF, &b_muonTrackIso_PF);
   fChain->SetBranchAddress("muonEcalIso_PF", &muonEcalIso_PF, &b_muonEcalIso_PF);
   fChain->SetBranchAddress("muonHcalIso_PF", &muonHcalIso_PF, &b_muonHcalIso_PF);
   fChain->SetBranchAddress("muonChi2_PF", &muonChi2_PF, &b_muonChi2_PF);
   fChain->SetBranchAddress("muonNdof_PF", &muonNdof_PF, &b_muonNdof_PF);
   fChain->SetBranchAddress("muonNhits_PF", &muonNhits_PF, &b_muonNhits_PF);
   fChain->SetBranchAddress("muonTrackd0_PF", &muonTrackd0_PF, &b_muonTrackd0_PF);
   fChain->SetBranchAddress("muonTrackPhi_PF", &muonTrackPhi_PF, &b_muonTrackPhi_PF);
   fChain->SetBranchAddress("muonPassID_PF", &muonPassID_PF, &b_muonPassID_PF);
   if (isV00_01_02())    fChain->SetBranchAddress("muonVtx_z_PF", &muonVtx_z_PF, &b_muonVtx_z_PF);
   fChain->SetBranchAddress("muonEcalVeto_PF", &muonEcalVeto_PF, &b_muonEcalVeto_PF);
   fChain->SetBranchAddress("muonHcalVeto_PF", &muonHcalVeto_PF, &b_muonHcalVeto_PF);
   fChain->SetBranchAddress("nMuons_PF", &nMuons_PF, &b_nMuons_PF);
   fChain->SetBranchAddress("eleEt_PF", &eleEt_PF, &b_eleEt_PF);
   fChain->SetBranchAddress("eleEta_PF", &eleEta_PF, &b_eleEta_PF);
   fChain->SetBranchAddress("eleTrackIso_PF", &eleTrackIso_PF, &b_eleTrackIso_PF);
   fChain->SetBranchAddress("eleEcalIso_PF", &eleEcalIso_PF, &b_eleEcalIso_PF);
   fChain->SetBranchAddress("eleHcalIso_PF", &eleHcalIso_PF, &b_eleHcalIso_PF);
   if (isV00_01_02())    fChain->SetBranchAddress("eledB_PF", &eledB_PF, &b_eledB_PF);
   if (isV00_01_02())    fChain->SetBranchAddress("eleVtx_z_PF", &eleVtx_z_PF, &b_eleVtx_z_PF);
   fChain->SetBranchAddress("eleIDLoose_PF", &eleIDLoose_PF, &b_eleIDLoose_PF);
   fChain->SetBranchAddress("eleIDRobustTight_PF", &eleIDRobustTight_PF, &b_eleIDRobustTight_PF);
   fChain->SetBranchAddress("elePassID_PF", &elePassID_PF, &b_elePassID_PF);
   fChain->SetBranchAddress("nElectrons_PF", &nElectrons_PF, &b_nElectrons_PF);
   fChain->SetBranchAddress("SUSY_nb", &SUSY_nb, &b_SUSY_nb);
   fChain->SetBranchAddress("qScale", &qScale, &b_qScale);
   fChain->SetBranchAddress("topDecayCode", &topDecayCode, &b_topDecayCode);
   fChain->SetBranchAddress("flavorHistory", &flavorHistory, &b_flavorHistory);

   if (isV00_01_05()) {
   fChain->SetBranchAddress("tauPt", &tauPt, &b_tauPt);
   fChain->SetBranchAddress("tauEta", &tauEta, &b_tauEta);
   fChain->SetBranchAddress("tauPhi", &tauPhi, &b_tauPhi);
   fChain->SetBranchAddress("tauTaNC", &tauTaNC, &b_tauTaNC);
   fChain->SetBranchAddress("tauID_againstElectron", &tauID_againstElectron, &b_tauID_againstElectron);
   fChain->SetBranchAddress("tauID_againstMuon", &tauID_againstMuon, &b_tauID_againstMuon);
   fChain->SetBranchAddress("tauID_byIsolation", &tauID_byIsolation, &b_tauID_byIsolation);
   fChain->SetBranchAddress("tauID_byTaNCfrOnePercent", &tauID_byTaNCfrOnePercent, &b_tauID_byTaNCfrOnePercent);
   fChain->SetBranchAddress("tauID_byTaNCfrHalfPercent", &tauID_byTaNCfrHalfPercent, &b_tauID_byTaNCfrHalfPercent);
   fChain->SetBranchAddress("tauID_byTaNCfrQuarterPercent", &tauID_byTaNCfrQuarterPercent, &b_tauID_byTaNCfrQuarterPercent);
   fChain->SetBranchAddress("tauPt_PF", &tauPt_PF, &b_tauPt_PF);
   fChain->SetBranchAddress("tauEta_PF", &tauEta_PF, &b_tauEta_PF);
   fChain->SetBranchAddress("tauPhi_PF", &tauPhi_PF, &b_tauPhi_PF);
   fChain->SetBranchAddress("tauTaNC_PF", &tauTaNC_PF, &b_tauTaNC_PF);
   fChain->SetBranchAddress("tauID_againstElectron_PF", &tauID_againstElectron_PF, &b_tauID_againstElectron_PF);
   fChain->SetBranchAddress("tauID_againstMuon_PF", &tauID_againstMuon_PF, &b_tauID_againstMuon_PF);
   fChain->SetBranchAddress("tauID_byIsolation_PF", &tauID_byIsolation_PF, &b_tauID_byIsolation_PF);
   fChain->SetBranchAddress("tauID_byTaNCfrOnePercent_PF", &tauID_byTaNCfrOnePercent_PF, &b_tauID_byTaNCfrOnePercent_PF);
   fChain->SetBranchAddress("tauID_byTaNCfrHalfPercent_PF", &tauID_byTaNCfrHalfPercent_PF, &b_tauID_byTaNCfrHalfPercent_PF);
   fChain->SetBranchAddress("tauID_byTaNCfrQuarterPercent_PF", &tauID_byTaNCfrQuarterPercent_PF, &b_tauID_byTaNCfrQuarterPercent_PF);
   }

   Notify();
}

Bool_t basicLoop::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void basicLoop::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

// ========================================== begin

bool basicLoop::isV00_01_02() {

  if (  (findInputName()).Contains("/V00-01-02/") ) return true;
  if (  (findInputName()).Contains("/V00-01-03/") ) return true;
  if (  (findInputName()).Contains("/V00-01-04/") ) return true;
  if (  (findInputName()).Contains("/V00-01-05/") ) return true;
  return false;

}

bool basicLoop::isV00_01_03() {

  if (  (findInputName()).Contains("/V00-01-03/") ) return true;
  if (  (findInputName()).Contains("/V00-01-04/") ) return true;
  if (  (findInputName()).Contains("/V00-01-05/") ) return true;
  return false;

}

bool basicLoop::isV00_01_05() {
  if (  (findInputName()).Contains("/V00-01-05/") ) return true;
  return false;
}

bool basicLoop::setCutScheme(CutScheme cutscheme) {

  if (cutscheme == nCutSchemes) return false;
  theCutScheme_ = cutscheme;

  cutNames_.clear();
  cutTags_.clear();
  cutMap_.clear();

  //cutNames and cutTags are required to 'line up' with each other.
  //one is for printability, the other is for unambiguous typing

  //this is for the kRA2 cut scheme
  //we want to run this no matter what in order to fill the cutMap correctly
  //the cutMap is invariant because the ntuple is always filled the same way

  cutTags_.push_back("cutInclusive"); cutNames_[ cutTags_.back()] = "Inclusive";
  cutTags_.push_back("cutTrigger");  cutNames_[cutTags_.back()]="Trigger";
  cutTags_.push_back("cutPV");  cutNames_[cutTags_.back()]="PV";
  cutTags_.push_back("cut3Jets");  cutNames_[cutTags_.back()]=">=3Jets";
  cutTags_.push_back("cutJetPt1");  cutNames_[cutTags_.back()]="JetPt1";
  cutTags_.push_back("cutJetPt2");  cutNames_[cutTags_.back()]="JetPt2";
  cutTags_.push_back("cutJetPt3");  cutNames_[cutTags_.back()]="JetPt3";
  cutTags_.push_back("cutHT");  cutNames_[cutTags_.back()]="HT";
  cutTags_.push_back("cutMET");  cutNames_[cutTags_.back()]="MET";
  cutTags_.push_back("cutMHT");  cutNames_[cutTags_.back()]="MHT";
  cutTags_.push_back("cutMuVeto");  cutNames_[cutTags_.back()]="MuVeto";
  cutTags_.push_back("cutEleVeto");  cutNames_[cutTags_.back()]="EleVeto";
  cutTags_.push_back("cutDeltaPhi");  cutNames_[cutTags_.back()]="DeltaPhi";
  cutTags_.push_back("cut1b");  cutNames_[cutTags_.back()]=">=1b";
  cutTags_.push_back("cut2b");  cutNames_[cutTags_.back()]=">=2b";
  cutTags_.push_back("cut3b");  cutNames_[cutTags_.back()]=">=3b";

  //read the above comment...this must be here (not lower)
  for (unsigned int i=0;i<cutTags_.size(); i++) {
    cutMap_[cutTags_[i]] = i;
  }
  
  if (theCutScheme_==kSync1) {
    cutNames_.clear();
    cutTags_.clear();
    //now the tags
    cutTags_.push_back("cutInclusive");cutNames_[ cutTags_.back()] = "Inclusive";
    cutTags_.push_back("cutTrigger"); cutNames_[cutTags_.back()]="Trigger";
    cutTags_.push_back("cutPV"); cutNames_[cutTags_.back()]="PV";
    cutTags_.push_back("cut3Jets"); cutNames_[cutTags_.back()]=">=3Jets";
    cutTags_.push_back("cutJetPt1"); cutNames_[cutTags_.back()]="JetPt1";
    cutTags_.push_back("cutJetPt2"); cutNames_[cutTags_.back()]="JetPt2";
    cutTags_.push_back("cutJetPt3");cutNames_[cutTags_.back()]="JetPt3";

    /*
doing it this way is a dirty hack, but it is so much easier than implementing a new cutScheme definition.
    */
    const bool modb=false;
    if (modb) {
      cutTags_.push_back("cut1b"); cutNames_[cutTags_.back()]=">=1b";
      cutTags_.push_back("cut2b"); cutNames_[cutTags_.back()]=">=2b";
      cutTags_.push_back("cut3b"); cutNames_[cutTags_.back()]=">=3b";
    }

    cutTags_.push_back("cutMuVeto"); cutNames_[cutTags_.back()]="MuVeto";
    cutTags_.push_back("cutEleVeto"); cutNames_[cutTags_.back()]="EleVeto";

    cutTags_.push_back("cutHT"); cutNames_[cutTags_.back()]="HT";

    cutTags_.push_back("cutMET"); cutNames_[cutTags_.back()]="MET";
    cutTags_.push_back("cutMHT"); cutNames_[cutTags_.back()]="MHT";

    cutTags_.push_back("cutDeltaPhi"); cutNames_[cutTags_.back()]="DeltaPhi";
    if (!modb) {
      cutTags_.push_back("cut1b"); cutNames_[cutTags_.back()]=">=1b";
      cutTags_.push_back("cut2b"); cutNames_[cutTags_.back()]=">=2b";
      cutTags_.push_back("cut3b"); cutNames_[cutTags_.back()]=">=3b";
    }

  }
  else if (theCutScheme_ == kBaseline0) {
    cutNames_.clear();
    cutTags_.clear();
    //now the tags
    // important...this is the correct order
    cutTags_.push_back("cutInclusive");cutNames_[ cutTags_.back()] = "Inclusive";
    cutTags_.push_back("cutTrigger"); cutNames_[cutTags_.back()]="Trigger";
    cutTags_.push_back("cutPV"); cutNames_[cutTags_.back()]="PV";
    cutTags_.push_back("cutHT"); cutNames_[cutTags_.back()]="HT";
    cutTags_.push_back("cut3Jets"); cutNames_[cutTags_.back()]=">=3Jets";
    cutTags_.push_back("cutJetPt1");  cutNames_[cutTags_.back()]="JetPt1";
    cutTags_.push_back("cutJetPt2");  cutNames_[cutTags_.back()]="JetPt2";
    cutTags_.push_back("cutJetPt3");  cutNames_[cutTags_.back()]="JetPt3";

    cutTags_.push_back("cutEleVeto");  cutNames_[cutTags_.back()]="EleVeto";
    cutTags_.push_back("cutMuVeto");  cutNames_[cutTags_.back()]="MuVeto";
    cutTags_.push_back("cutTauVeto");  cutNames_[cutTags_.back()]="TauVeto";
    
    cutTags_.push_back("cutMET");  cutNames_[cutTags_.back()]="MET";

    cutTags_.push_back("cutDeltaPhi"); cutNames_[cutTags_.back()]="DeltaPhi";
    cutTags_.push_back("cut1b"); cutNames_[cutTags_.back()]=">=1b";
    cutTags_.push_back("cut2b"); cutNames_[cutTags_.back()]=">=2b";
    cutTags_.push_back("cut3b"); cutNames_[cutTags_.back()]=">=3b";
    

    //modified order (for comparisons with Don)
    /*
    cout<<"Using modified Baseline0 cut order!"<<endl;
    cutTags_.push_back("cutInclusive");cutNames_[ cutTags_.back()] = "Inclusive";
    cutTags_.push_back("cutTrigger"); cutNames_[cutTags_.back()]="Trigger";
    cutTags_.push_back("cutPV"); cutNames_[cutTags_.back()]="PV";
    cutTags_.push_back("cut3Jets"); cutNames_[cutTags_.back()]=">=3Jets";
    cutTags_.push_back("cutHT"); cutNames_[cutTags_.back()]="HT";

    cutTags_.push_back("cutJetPt1");  cutNames_[cutTags_.back()]="JetPt1";
    cutTags_.push_back("cutJetPt2");  cutNames_[cutTags_.back()]="JetPt2";
    cutTags_.push_back("cutJetPt3");  cutNames_[cutTags_.back()]="JetPt3";

    cutTags_.push_back("cutMET");  cutNames_[cutTags_.back()]="MET";
    cutTags_.push_back("cutDeltaPhi"); cutNames_[cutTags_.back()]="DeltaPhi";

    cutTags_.push_back("cut1b"); cutNames_[cutTags_.back()]=">=1b";
    cutTags_.push_back("cut2b"); cutNames_[cutTags_.back()]=">=2b";
    cutTags_.push_back("cut3b"); cutNames_[cutTags_.back()]=">=3b";

    cutTags_.push_back("cutEleVeto");  cutNames_[cutTags_.back()]="EleVeto";
    cutTags_.push_back("cutMuVeto");  cutNames_[cutTags_.back()]="MuVeto";
    cutTags_.push_back("cutTauVeto");  cutNames_[cutTags_.back()]="TauVeto";
    */
  }

  //debug
  //  for (unsigned int i=0;i<cutTags_.size(); i++) {
  //    cout<<cutTags_[i]<<"\t"<<cutNames_[ cutTags_[i]]<<endl;
  //  }

  return true;
}

bool basicLoop::cutRequired(const TString cutTag) { //should put an & in here to improve performance


  //check if we are *ignoring* this cut (for N-1 studies, etc)
  for (unsigned int i = 0; i< ignoredCut_.size() ; i++) {
    if ( cutTag == ignoredCut_.at(i) ) return false;
  }

  bool cutIsRequired=false;

  //RA2
  if (theCutScheme_ == kRA2 ) {
    if      (cutTag == "cutInclusive")  cutIsRequired =  true;
    else if (cutTag == "cutTrigger")  cutIsRequired =  true;
    else if (cutTag == "cutPV")  cutIsRequired =  true;
    else if (cutTag == "cut3Jets")  cutIsRequired =  true;
    else if (cutTag == "cutJetPt1")  cutIsRequired =  false;
    else if (cutTag == "cutJetPt2")  cutIsRequired =  false;
    else if (cutTag == "cutJetPt3")  cutIsRequired =  false;
    else if (cutTag == "cutHT")  cutIsRequired =  true;
    //MET
    else if (cutTag == "cutMET")  cutIsRequired =  (theMETType_== kMET || theMETType_==ktcMET || theMETType_==kpfMET);
    //MHT
    else if (cutTag == "cutMHT")  cutIsRequired =  (theMETType_==kMHT);
    else if (cutTag == "cutMuVeto") cutIsRequired =  true;
    else if (cutTag == "cutEleVeto") cutIsRequired =  true;
    else if (cutTag == "cutDeltaPhi") cutIsRequired =  true;
    else if (cutTag == "cut1b") cutIsRequired =  nBcut_ >=1;
    else if (cutTag == "cut2b") cutIsRequired =  nBcut_ >=2;
    else if (cutTag == "cut3b") cutIsRequired =  nBcut_ >=3;
    else assert(0);
  }
  else if (theCutScheme_==kSync1) {
    if      (cutTag == "cutInclusive")  cutIsRequired =  true;
    else if (cutTag == "cutTrigger")  cutIsRequired =  false;
    else if (cutTag == "cutPV")  cutIsRequired =  true;
    else if (cutTag == "cut3Jets")  cutIsRequired =  true;
    else if (cutTag == "cutJetPt1")  cutIsRequired =  true;
    else if (cutTag == "cutJetPt2")  cutIsRequired =  true;
    else if (cutTag == "cutJetPt3")  cutIsRequired =  true;
    else if (cutTag == "cutHT")  cutIsRequired =  true;
    //MET
    else if (cutTag == "cutMET")  cutIsRequired =  (theMETType_== kMET || theMETType_==ktcMET ||theMETType_==kpfMET);
    //MHT
    else if (cutTag == "cutMHT")  cutIsRequired =  (theMETType_==kMHT);
    else if (cutTag == "cutMuVeto") cutIsRequired =  true;
    else if (cutTag == "cutEleVeto") cutIsRequired =  true;
    else if (cutTag == "cutDeltaPhi") cutIsRequired =  true;
    else if (cutTag == "cut1b") cutIsRequired =  nBcut_ >=1;
    else if (cutTag == "cut2b") cutIsRequired =  nBcut_ >=2;
    else if (cutTag == "cut3b") cutIsRequired =  nBcut_ >=3;
    else assert(0);
  }
  else if (theCutScheme_==kBaseline0) {
    if      (cutTag == "cutInclusive")  cutIsRequired =  true;
    else if (cutTag == "cutTrigger")  cutIsRequired = true;
    else if (cutTag == "cutPV")  cutIsRequired =  true;
    else if (cutTag == "cutHT")  cutIsRequired =  true;
    else if (cutTag == "cut3Jets")  cutIsRequired =  true;
    else if (cutTag == "cutJetPt1")  cutIsRequired =  false;
    else if (cutTag == "cutJetPt2")  cutIsRequired =  false;
    else if (cutTag == "cutJetPt3")  cutIsRequired =  false;
    //roll MHT into MET
    else if (cutTag == "cutMET")  cutIsRequired =  (theMETType_== kMET || 
						    theMETType_==ktcMET ||
						    theMETType_==kpfMET ||
						    theMETType_==kMHT);
    else if (cutTag == "cutMuVeto") cutIsRequired =  true;
    else if (cutTag == "cutEleVeto") cutIsRequired =  true;
    else if (cutTag == "cutTauVeto") cutIsRequired =  false; //not required for now
    else if (cutTag == "cutDeltaPhi") cutIsRequired =  true;
    else if (cutTag == "cut1b") cutIsRequired =  nBcut_ >=1;
    else if (cutTag == "cut2b") cutIsRequired =  nBcut_ >=2;
    else if (cutTag == "cut3b") cutIsRequired =  nBcut_ >=3;
    else assert(0);
  }
  else assert(0);

  return cutIsRequired;
}


bool basicLoop::passMuVetoSync1() {
  //this is less fast than returning the veto as soon as we get to '1', but avoids code duplication

  //in principle I can add a special argument to countMuSync1 to return as soon as it gets to 1,
  //but let's leave that optimization for later

  //using nmu_ allows for inverting the veto
  return (countMuSync1() == nmu_);
}

int basicLoop::countMuSync1() {

  int ngoodmu=0;

  unsigned int nmu = 0;
  if (theLeptonType_ == kNormal) nmu=muonPt->size();  
  else if (theLeptonType_ ==kPFLeptons) nmu=muonPt_PF->size();
  else {assert(0);}

  //need to be AllGlobalMuons
  //this is all I stored in the ntuple
  for ( unsigned int i = 0; i< nmu; i++) {
    float pt=0,eta=0,reliso=0;
    
    //load information depending on lepton type
    if (theLeptonType_ == kNormal) {
      pt = muonPt->at(i);
      eta = muonEta->at(i);
      reliso = (muonTrackIso->at(i) + muonHcalIso->at(i) + muonEcalIso->at(i))/muonPt->at(i);
    }
    else if (theLeptonType_ ==kPFLeptons) {
      pt = muonPt_PF->at(i);
      eta = muonEta_PF->at(i);
      reliso = (muonTrackIso_PF->at(i) + muonHcalIso_PF->at(i) + muonEcalIso_PF->at(i))/muonPt_PF->at(i);
    }
    else {assert(0);}

    //now make cuts
    if ( pt < 10 ) continue;
    if ( fabs(eta) > 2.5 ) continue;
    if ( reliso > 0.2) continue;

    ++ngoodmu;
  }
  
  return ngoodmu;
}

bool basicLoop::passTauVeto() {
  //TO DO
  return true;
}

bool basicLoop::passEleVetoSync1() {
  //see comments in muon section

  //using ne_ allows for inverting the veto
  return (countEleSync1() == ne_);
}

int basicLoop::countEleSync1() {

  int ngoodele=0;

  unsigned int nele = 0;
  if (theLeptonType_ == kNormal) nele=eleEt->size();  
  else if (theLeptonType_ ==kPFLeptons) nele=eleEt_PF->size();
  else {assert(0);}

  for (unsigned int i=0; i < nele; i++) {
    float et=0,eta=0,reliso=0;
    
    //load information depending on lepton type
    if (theLeptonType_ == kNormal) {
      et = eleEt->at(i);
      eta = eleEta->at(i);
      reliso = (eleTrackIso->at(i) + eleHcalIso->at(i) + eleEcalIso->at(i))/eleEt->at(i);
    }
    else if (theLeptonType_ ==kPFLeptons) {
      et = eleEt_PF->at(i);
      eta = eleEta_PF->at(i);
      reliso = (eleTrackIso_PF->at(i) + eleHcalIso_PF->at(i) + eleEcalIso_PF->at(i))/eleEt_PF->at(i);
    }
    else {assert(0);}

    //make cuts    
    if ( et < 15 ) continue;
    if ( fabs(eta) > 2.5 ) continue;
    if ( reliso > 0.2) continue;
    
    //if any electron passes all of these cuts, then it is good
    ++ngoodele;
  }

  return ngoodele;

}

bool basicLoop::passMuVetoRA2() {

  /* debugging
  const TString sp=" ";
  cout<<"[passMuVetoRA2] "<<flush;
  cout<<sp<<muonIsGlobalMuonPromptTight->size();
  cout<<sp<<muonIsAllGlobalMuons->size();
  cout<<sp<<  muonPt->size();
  cout<<sp<<  muonEta->size();
  cout<<sp<<  muonTrackIso->size();
  cout<<sp<<  muonEcalIso->size();
  cout<<sp<<  muonHcalIso->size();
  cout<<sp<<  muonChi2->size();
  cout<<sp<<  muonNdof->size();
  cout<<sp<<  muonNhits->size();
  cout<<sp<<  muonTrackd0->size();
  cout<<sp<<  muonTrackPhi->size();
  cout<<sp<<  muonPassID->size();
  cout<<sp<<   muonEcalVeto->size();
  cout<<sp<<   muonHcalVeto->size()<<endl;
  */  

  //loop over regular (not PF) muons
  for ( unsigned int i = 0; i< muonPt->size(); i++) {
    
    //due to a bug, this isn't going to work
    //   if ( !muonIsGlobalMuonPromptTight->at(i) ) continue;
    //here is a sorry attempt to cover part of my ass
    if ( muonIsGlobalMuonPromptTight->size() ==0) continue;

    if ( muonPt->at(i) <= 10 ) continue;
    if ( fabs(muonEta->at(i)) > 2.4 ) continue;

    float ndof= muonNdof->at(i);
    if (ndof<=0) continue;
    if ( muonChi2->at(i) / ndof >= 10) continue;

    if ( fabs(muonTrackd0->at(i)) > 0.2) continue;

    if (muonNhits->at(i) <11) continue;

    if ( (muonTrackIso->at(i) + muonHcalIso->at(i) + muonEcalIso->at(i))/muonPt->at(i) > 0.1) continue;

    //if any muon passes all of these cuts, then veto
    return false;
  }

  return true;
}

bool basicLoop::passEleVetoRA2() {

  //loop over regular (not PF) electrons
  for (unsigned int i=0; i< eleIDLoose->size(); i++) {

    if ( !(eleIDLoose->at(i) >0) ) continue;

    if ( eleEt->at(i) < 15 ) continue;
    if ( fabs(eleEta->at(i)) > 2.4 ) continue;

    //rel iso of 0.1 and d0<0.2 are in here
    if (!elePassID->at(i)) continue;

    //if ( (eleTrackIso->at(i) + eleHcalIso->at(i) + eleEcalIso->at(i))/eleEt->at(i) > 0.1) continue;
    //i don't have d0 in the ntuple

    //if any electron passes all of these cuts, then veto
    return false;
  }
  return true;
}

bool basicLoop::passDeltaPhi_Sync1() {

  //MET and jet type will adjust with the global settings!

  int njet=0;

  bool pass=true;

  for (unsigned int i=0; i<loosejetPt->size(); i++) {

    if (isGoodJet_Sync1(i)) {
      njet++;

      if (njet == 1) {
	if ( getDeltaPhi(loosejetPhi->at(i),getMETphi()) <=0.3) pass=false;
      }
      else if (njet == 2) {
	if ( getDeltaPhi(loosejetPhi->at(i),getMETphi()) <=0.35) pass=false;
      }

    }
  }
  return pass;
}

bool basicLoop::passHLT() { //do I pass the OR of what is in triggerList_

  for (std::map<TString, int>::const_iterator itrig=triggerList_.begin(); 
       itrig!=triggerList_.end(); ++itrig) {

    if ( passTrigger->at( itrig->second) ) {
      //      lastTriggerPass_ = itrig->first; //may want to remove this code sooner or later
      //      if (!printedHLT_) {
      //	cout<<"Pass: "<<itrig->first<<" "<<hltPrescale->at(itrig->second)<<endl;
      //	printedHLT_=true;
      //      }
      return true;
    }
  }
  return false;
}

bool basicLoop::passHLT(const TString & triggerName) { //do I pass the trigger given as the argument?

  if (triggerListAll_.find(triggerName) != triggerListAll_.end() ) {
    int index = triggerListAll_[triggerName];
    //    cout<<triggerName<<"\t"<<passTrigger->at(index)<<"\t"<<hltPrescale->at(index)<<endl;
    return passTrigger->at(index);
  }
  //only get here if we don't find the trigger
  cout<<triggerName<<" not found!"<<endl;
  return false;
  
}

int basicLoop::getHLTPrescale(const TString & triggerName) {

  if (triggerListAll_.find(triggerName) != triggerListAll_.end() ) {
    int index = triggerListAll_[triggerName];
    //    cout<<triggerName<<"\t"<<passTrigger->at(index)<<"\t"<<hltPrescale->at(index)<<endl;
    return hltPrescale->at(index);
  }
  //only get here if we don't find the trigger
  cout<<triggerName<<" not found!"<<endl;
  return -1;
  
}

bool basicLoop::passCut(const TString cutTag) {
  //implement special exceptions here

  const bool debug=false;  if (debug) cout<<cutTag<<endl;

  //trigger cut -- same for all cut schemes
  if (cutTag=="cutTrigger" ) return passHLT();

  if (theCutScheme_ == kBaseline0) {
    if (cutTag=="cutPV") return passPV();

    if (cutTag=="cutHT") return getHT()>300;

    if (cutTag=="cut3Jets") return (nGoodJets() >= 3);

    if (cutTag=="cutJetPt1") return jetPtOfN(1)>100; //this cut is not actually required at the moment

    //deltaPhi done below (only updated for kminDP !)
    //lepton vetoes done below
    //MET done below
    //b tagging done below
  }

  //lepton vetoes
  if (cutTag=="cutMuVeto" && theCutScheme_==kRA2) return passMuVetoRA2();
  else if (cutTag=="cutMuVeto" && theCutScheme_==kSync1) return passMuVetoSync1();
  else if (cutTag=="cutMuVeto" && theCutScheme_==kBaseline0) return passMuVetoSync1(); //passMuVetoDon()

  if (cutTag=="cutEleVeto" && theCutScheme_==kRA2) return passEleVetoRA2();
  else if (cutTag=="cutEleVeto" && theCutScheme_==kSync1) return passEleVetoSync1();
  else if (cutTag=="cutEleVeto" && theCutScheme_==kBaseline0) return passEleVetoSync1(); //passEleVetoDon()

  if (cutTag=="cutTauVeto" && theCutScheme_==kBaseline0) return passTauVeto();

  if (cutTag=="cut3Jets" && theCutScheme_==kSync1) return (nGoodJets_Sync1() >= 3);

  if (theCutScheme_==kSync1 ) {
    if (cutTag=="cutJetPt1") return jetPtOfN(1)>100;
    if (cutTag=="cutJetPt2") return jetPtOfN(2)>50;
    if (cutTag=="cutJetPt3") return jetPtOfN(3)>50;;
  }

  if (cutTag=="cut3Jets" && theJetType_==kPF && theCutScheme_==kRA2) return (tightJetIndex->size() >= 3);

  if (cutTag=="cutHT" && theJetType_==kPF && theCutScheme_==kRA2) {
    //calculate HT
    double pfHT = 0;
    for (unsigned int i=0; i<tightJetIndex->size(); i++){
      pfHT+=jetPt.at(i);
    }
    return (pfHT>300.);
  }
  else if (cutTag=="cutHT" && theCutScheme_==kSync1) return getHT_Sync1() >300;

  //MET
  //it is now an anachronism that we treat MET and MHT as different cut flow steps
  //this block of code should now evaluate *the same* for both steps
  if (cutTag == "cutMET" || cutTag=="cutMHT") {
    float mymet = 0;
    if      ( theMETType_ == kMET )   mymet=caloMET;
    else if ( theMETType_ == ktcMET ) mymet=tcMET;
    else if ( theMETType_ == kpfMET ) mymet=pfMET;
    //FIXME...loose jet def'n has changed, so MHT stored in ntuple is probably not ok
    else if ( theMETType_ == kMHT)    mymet=MHT;
    else {assert(0);}

    if (theMETRange_ == kMedium) return (mymet>=50 && mymet<100); //redefining kMedium as 50-100!
    //on the next line we are redoing the default that is stored in the ntuple. I think this is ok
    else if (theMETRange_ == kHigh) return (mymet >=150);
    else if (theMETRange_ == kWide) return (mymet >=50);
    else if (theMETRange_ == kMedhigh) return (mymet>=100 && mymet<150);
    else {assert(0);}
  }
  
  // -- in case we are using MET instead of MHT, we need to alter the DeltaPhi cuts --
  if (cutTag == "cutDeltaPhi" && 
      (theDPType_ == kDeltaPhi && theMETType_!=kMHT)) {
    //      ( theCutScheme_ == kRA2MET ||theCutScheme_ == kRA2tcMET ) ) {
    
    float phi_of_MET = getMETphi();
    
    // FIXME loose jet def'n has changed! this needs a careful update
    int nloosejets = loosejetPhi_calo->size();
    //need to calculate DeltaPhi between jets and MET
    //for RA2 and MHT, this is done with *loose* jets!
    double dp0=0,dp1=0,dp2=0;
    if (nloosejets>0) {dp0=getDeltaPhi( loosejetPhi_calo->at(0), phi_of_MET);} else {return false;}
    if (nloosejets>1) {dp1=getDeltaPhi( loosejetPhi_calo->at(1), phi_of_MET);} else {return false;}
    if (nloosejets>2) {dp2=getDeltaPhi( loosejetPhi_calo->at(2), phi_of_MET);} else {return false;}
    //here is the implementation of the DeltaPhi cuts
    if ( dp0 >0.3 && dp1 >0.5 && dp2 >0.3 ) { return true; } else {return false;}
  }
  else if (cutTag == "cutDeltaPhi" && //replace normal DeltaPhi cut with minDeltaPhi cut
	   ( theDPType_ == kminDP ) ) {
    
    if (theMETType_ ==kMET ||theMETType_==ktcMET ||theMETType_==kpfMET ||theMETType_ == kMHT) return ( getMinDeltaPhiMET(3) >= 0.3 );
    else {assert(0);}
  }
  else if (cutTag == "cutDeltaPhi" && //replace normal DeltaPhi cut with minDeltaPhi cut
	   ( theDPType_ == kminDPinv ) ) { //this is the *inverted* minDeltaPhi cut (for control region)
    
    if (theMETType_ ==kMET ||theMETType_==ktcMET ||theMETType_==kpfMET ||theMETType_ == kMHT) return ( getMinDeltaPhiMET(3) < 0.3 );
    else {assert(0);}
  }
  else if (cutTag == "cutDeltaPhi" && //replace normal DeltaPhi cut with DeltaPhi(MPT,MET) cut
	   ( theDPType_ == kMPT ) ) {
    if (theMETType_ == kMET)    return ( getDeltaPhiMPTMET() < 2.0 );
    else {cout<<"DeltaPhiMPTMET not implemented for that MET type!"<<endl; assert(0);}
  }
  else if (cutTag == "cutDeltaPhi" && (theDPType_==kDPSync1)) return passDeltaPhi_Sync1();

  //compensate for bugs in nbSSVM variable (fixed in this code)
  if (theCutScheme_==kRA2 || theCutScheme_==kBaseline0) {
    if (cutTag == "cut1b") return nbSSVM >=1;
    if (cutTag == "cut2b") return nbSSVM >=2;
    if (cutTag == "cut3b") return nbSSVM >=3;
  }
  else if (theCutScheme_==kSync1) {
    if (cutTag == "cut1b" || cutTag == "cut2b" || cutTag == "cut3b") {
      int nb = countBJets_Sync1();
      if (cutTag == "cut1b") return nb >=1;
      if (cutTag == "cut2b") return nb >=2;
      if (cutTag == "cut3b") return nb >=3;
    }
  }

  //we want to avoid using the precomputed cut flow for the basline0 scheme
  if (theCutScheme_==kBaseline0 && cutTag!="cutInclusive") {
    cout<<"[passCut] should not have reached this point for Baseline0 scheme! "<<cutTag<<endl;;
    assert(0);
  }

  //in case this is not an exception, return the cut result stored in the ntuple
  int cutIndex = cutMap_[cutTag];
  return cutResults->at(cutIndex);
}

bool basicLoop::passPV() {


  //i have now understood the discrepancies between the PVSelector and the hand-calculated result
  //the PV selector only looks at the first PV. I think this is wrong.
  //so this function should be used instead

  //update...Don endorses the PVSelector method.
  bool ntupleResult = cutResults->at( cutMap_["cutPV"]);
  return ntupleResult;

  //never run the code down here!

  bool pass=false;
  for (unsigned int ipv = 0; ipv<pv_isFake->size(); ipv++) {
    if ( pv_isFake->at(ipv) ) continue;
    if ( fabs(pv_z->at(ipv)) > 24 ) continue;
    if ( fabs(pv_rho->at(ipv)) > 2 ) continue;
    if ( pv_ndof->at(ipv) <= 4 ) continue;

    pass=true; break;
  }

  /*
  if (ntupleResult != pass) {
    cout<<"PV conflict! ntupleResult = "<<ntupleResult<<" passPV results = "<<pass<<endl;
    for (unsigned int ipv = 0; ipv<pv_isFake->size(); ipv++) {
      cout<<ipv<<"\t"
	  <<pv_isFake->at(ipv) <<" "
	  <<fabs(pv_z->at(ipv))<<" "
	  <<fabs(pv_rho->at(ipv))<<" "
	  <<pv_ndof->at(ipv)<<endl;
    }
  }
  */

  return pass; 
}

int basicLoop::countBJets() {
  int nb=0;

  for (unsigned int i=0; i<loosejetPt->size(); i++) {
    if (isGoodJet30( i) ) {
      if ( passSSVM(i) ) nb++;
    }
  }
  return nb;
}

bool basicLoop::passSSVM(int i) {
  return ( loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags->at(i) >= 1.74 
	   || loosejetBTagDisc_simpleSecondaryVertexBJetTags->at(i) >=1.74 );
}

int basicLoop::countBJets_Sync1() {

  int nb=0;

  for (unsigned int i=0; i<loosejetBTagDisc_trackCountingHighPurBJetTags->size(); i++) {

    if (loosejetPt->at(i) <30) continue;
    if (fabs(loosejetEta->at(i)) >2.4) continue;

    if ( loosejetBTagDisc_trackCountingHighPurBJetTags->at(i) >3) nb++;

  }
  return nb;
}

Int_t basicLoop::Cut(Long64_t entry)
{
  //be careful...not every function that makes cuts uses this! e.g. ::cutflow()

  for (unsigned int i=0; i< cutTags_.size(); i++) {
    if (cutRequired( cutTags_[i] ) && !passCut( cutTags_[i]) ) return -1;
  }
  
  return 1;
  
}

double basicLoop::getUncorrectedHT(const double threshold) {

  /*
I have realized that if I want to do a crude trigger simulation, then I need to store every jet
with Uncorrect Pt>20 GeV, not just jets passing my 'loose' cuts

   vector<float>   *loosejetPtUncorr;
   vector<float>   *loosejetEtaUncorr;
   vector<float>   *loosejetPhiUncorr;*/
  double ht=0;
  for (unsigned int ijet=0; ijet<veryloosejetPtUncorr->size(); ijet++) {

    if ( fabs(veryloosejetEtaUncorr->at(ijet)) > 5 ) continue; //what is the eta range actually used by the HLT?

    double thispt = veryloosejetPtUncorr->at(ijet);
    if (thispt>=threshold)   ht+= thispt;
  }
  
   return ht;
}

float basicLoop::getMET() {
  //for just one number i am happy to do this by value

  //adjust for global MET type
  if (theMETType_== kMET) {
    return caloMET;
  }
  else if (theMETType_ == ktcMET) {
    return tcMET;
  }
  else if (theMETType_ == kpfMET) {
    return pfMET;
  }
  else if (theMETType_ == kMHT) {
    return getMHT();
  }
  
  assert(0);

  return 0;
}

float basicLoop::getMETphi() {
  //for just one number i am happy to do this by value

  //adjust for global MET type
  if (theMETType_== kMET) {
    return caloMETphi;
  }
  else if (theMETType_ == ktcMET) {
    return tcMETphi;
  }
  else if (theMETType_ == kpfMET) {
    return pfMETphi;
  }
  else if (theMETType_ == kMHT) {
    return getMHTphi();
  }
  
  assert(0);

  return 0;
}

double basicLoop::getMinDeltaPhiMET(unsigned int maxjets) {
  /*
code is now bifurcated to:
    use tight jets of the default type;
    unless we are using kBaseline0, in which case we use isGoodJet() to find the good jets
  */

  double mindp=99;

  if (theCutScheme_==kBaseline0) {
    unsigned int ngood=0;
    //get the minimum angle between the first n jets and MET
    for (unsigned int i=0; i< loosejetPhi->size(); i++) {

      if (isGoodJet(i)) {
	++ngood;
	double dp =  getDeltaPhi( loosejetPhi->at(i) , getMETphi());
	if (dp<mindp) mindp=dp;
	if (ngood >= maxjets) break;
      }
    }
  }
  //here is the "legacy" way of doing things
  else {
    unsigned int njets=  jetPhi.size();
    if (njets < maxjets) maxjets = njets;
    
    //get the minimum angle between the first n jets and MET
    for (unsigned int i=0; i< maxjets; i++) {
      
      double dp =  getDeltaPhi( jetPhi.at(i) , getMETphi());
      if (dp<mindp) mindp=dp;
      
    }
  }

  return mindp;
}


double basicLoop::getMinDeltaPhibMET() {
  //get the minimum angle between a b jet and MET
  /* uses default tight jets and default MET */

  double mindp=99;
  for (unsigned int i=0; i<jetPhi.size(); i++) {
    if (passBCut(i) ) {
      double dp =  getDeltaPhi( jetPhi.at(i), getMETphi());
      if (dp<mindp) mindp=dp;
    }
  }
  return mindp;
}

double basicLoop::getDeltaPhib1b2() {
  //get the angle between the lead 2 b jets
  /* legacy code uses tight jets ; kBaseline0 uses loose jets and b cuts*/

  std::vector<float> phis;

  //careful...we're using tight jets for one case and loose jets for the other
  unsigned int maxj = theCutScheme_==kBaseline0 ? loosejetPhi->size() : jetPhi.size();

  for (unsigned int i=0; i< maxj; i++) {
    bool isGoodB=false;
    if (theCutScheme_==kBaseline0)  isGoodB =  isGoodJet30(i) && passSSVM(i);
    else                            isGoodB = passBCut(i);
     
    if (isGoodB ) {
      
      float thisJetPhi = theCutScheme_==kBaseline0 ? loosejetPhi->at(i) : jetPhi.at(i) ;
      phis.push_back( thisJetPhi);

      if (phis.size() == 2) break;
    }
  }

  //this is then invariant between cut schemes and such
  return getDeltaPhi(phis.at(0),phis.at(1));
}

double basicLoop::getDeltaPhiMPTMET() {
  /* uses default MET */

   //find MPT
   double MPTx=0;
   double MPTy=0;
   //we have already made minimal cuts on track pt, eta
   for (unsigned int i=0; i< trackPt->size(); i++) {
     
     MPTx -= trackPt->at(i) * cos(trackPhi->at(i));
     MPTy -= trackPt->at(i) * sin(trackPhi->at(i));
   }

   double MPTphi = atan2(MPTy,MPTx);

   return getDeltaPhi(getMETphi(), MPTphi);
}

bool basicLoop::passBCut( unsigned int bindex) {
  //bindex is taken to be the index of a tight jet

  //use tagger simpleSecondaryVertexHighEffBJetTags
  //which in older samples is called simpleSecondaryVertexBJetTags

  bool pass=false;
  float oldval=jetBTagDisc_simpleSecondaryVertexBJetTags.at(bindex);
  float newval=jetBTagDisc_simpleSecondaryVertexHighEffBJetTags.at(bindex);

  //if the tagger is not valid, there seems to be a value -1000 returned by CMSSW
  //so clearly that will fail the cut
  if ( (oldval >=1.74) || (newval >=1.74)) pass=true;

  return  pass;
}

double basicLoop::getDeltaPhi(double phi1, double phi2) {

  return acos(cos(phi1-phi2));
}

double basicLoop::getMinDeltaPhi_bj(unsigned int bindex) {
  /*
this function assumes that bindex is of a b jet. it doesn't verify it
  */

  //using tight jets

  double bphi = jetPhi.at(bindex);

  double minDeltaPhi=99;
  //loop over the jets
  for (unsigned int jindex=0; jindex<jetPhi.size(); jindex++) {
    if ( !passBCut(jindex) ) { //only look at non-b jets
      double dp = getDeltaPhi(bphi, jetPhi.at(jindex));
      if (dp < minDeltaPhi) minDeltaPhi = dp;
    }
  }

  return minDeltaPhi;
}

double basicLoop::getOverallMinDeltaR_bj() {

  double minDeltaR_bj=999;
  //note that all tight jet vectors should have the same size
  for (unsigned int ib = 0; ib< jetPhi.size(); ib++) {
    if ( passBCut(ib)) { //refind the b jets
      double mdr=getMinDeltaR_bj(ib); 
      if (mdr<minDeltaR_bj) minDeltaR_bj=mdr;
    }
  }

  return minDeltaR_bj;
}

double basicLoop::getMinDeltaR_bj(unsigned int bindex) {

  double beta=jetEta.at(bindex);
  double bphi=jetPhi.at(bindex);

  double minDeltaR = 999;

  //loop over the jets
  for (unsigned int jindex=0; jindex<jetPhi.size(); jindex++) {
    if ( !passBCut(jindex) ) { //only look at non-b jets
      double dr = pow(getDeltaPhi(bphi, jetPhi.at(jindex)),2) + pow(beta - jetEta.at(jindex),2);
      dr = sqrt(dr);
      if (dr < minDeltaR) minDeltaR = dr;
    }
  }

  return minDeltaR;
}

float basicLoop::getMHT() {
  //use isGoodJet() to recalculate it for the specified jet type

  double mhtx=0;
  double mhty=0;

  for (unsigned int i=0; i<loosejetPt->size(); i++) {
    if (isGoodJet( i ) ) {

      mhtx -= loosejetPt->at(i) * cos(loosejetPhi->at(i));
      mhty -= loosejetPt->at(i) * sin(loosejetPhi->at(i));
    }
  }
  
  return sqrt(mhtx*mhtx + mhty*mhty);
}

float basicLoop::getMHTphi() {
  //use isGoodJet() to recalculate it for the specified jet type

  double mhtx=0;
  double mhty=0;

  for (unsigned int i=0; i<loosejetPt->size(); i++) {
    if (isGoodJet( i ) ) {

      mhtx -= loosejetPt->at(i) * cos(loosejetPhi->at(i));
      mhty -= loosejetPt->at(i) * sin(loosejetPhi->at(i));
    }
  }

  return atan2(mhty,mhtx);
}

float basicLoop::getHT() {
  //use isGoodJet() to recalculate it for the specified jet type

  float ht=0;

  for (unsigned int i=0; i<loosejetPt->size(); i++) {
    if (isGoodJet( i ) ) ht+= loosejetPt->at(i);
  }

  return ht;
}

float basicLoop::getHT_Sync1() {
  //it is unclear to me if people want to do this with my proposal
  //(identical to the tight jets in the ntuple)
  //or using the good jet def'n as given by isGoodJet_Sync1()

  //can run it both ways (for LM0) and compare!
  float ht=0;

  //easy def'n using the ntuple cuts
  //  for (unsigned int i=0; i<jetPt.size() ; i++) ht += jetPt.at(i);

  for (unsigned int i=0; i<loosejetPt->size(); i++) {
    if (isGoodJet_Sync1( i ) ) ht+= loosejetPt->at(i);
  }

  return ht;
}

bool basicLoop::isGoodJet(unsigned int ijet) {

  if ( loosejetPt->at(ijet) <50) return false;
  if ( fabs(loosejetEta->at(ijet)) > 2.4) return false;
  if ( !(loosejetPassLooseID->at(ijet)) ) return false;

  return true;
}

bool basicLoop::isGoodJet30(unsigned int ijet) {

  if ( loosejetPt->at(ijet) <30) return false;
  if ( fabs(loosejetEta->at(ijet)) > 2.4) return false;
  if ( !(loosejetPassLooseID->at(ijet)) ) return false;

  return true;
}

bool basicLoop::isGoodJet_Sync1(unsigned int ijet) {

  if ( loosejetPt->at(ijet) <50) return false;
  
  if ( fabs(loosejetEta->at(ijet)) > 2.4) return false;
  
  if (loosejetEnergyFracHadronic->at(ijet) <0.1) return false;

  return true;
}

bool basicLoop::isLooseJet_Sync1(unsigned int ijet) {

  if ( loosejetPt->at(ijet) <30) return false;
  
  if ( fabs(loosejetEta->at(ijet)) > 2.4) return false;
  
  if (loosejetEnergyFracHadronic->at(ijet) <0.1) return false;

  return true;
}

float basicLoop::jetPtOfN(unsigned int n) {
  //updated to use Sync1 or Baseline0 cuts

  unsigned int ngood=0;
  for (unsigned int i=0; i<loosejetPt->size(); i++) {

    bool pass=false;
    if (theCutScheme_==kSync1) pass = isGoodJet_Sync1(i);
    else  if (theCutScheme_==kBaseline0) pass = isGoodJet(i);
    else {assert(0);}

    if (pass ) {
      ngood++;
      if (ngood==n) return  loosejetPt->at(i);
    }
  }
  return 0;
}

float basicLoop::jetPhiOfN(unsigned int n) {
  //updated to use Sync1 or Baseline0 cuts

  unsigned int ngood=0;
  for (unsigned int i=0; i<loosejetPt->size(); i++) {

    bool pass=false;
    if (theCutScheme_==kSync1) pass = isGoodJet_Sync1(i);
    else  if (theCutScheme_==kBaseline0) pass = isGoodJet(i);
    else {assert(0);}

    if (pass ) {
      ngood++;
      if (ngood==n) return  loosejetPhi->at(i);
    }
  }
  return 0;
}

float basicLoop::jetEtaOfN(unsigned int n) {
  //updated to use Sync1 or Baseline0 cuts

  unsigned int ngood=0;
  for (unsigned int i=0; i<loosejetPt->size(); i++) {

    bool pass=false;
    if (theCutScheme_==kSync1) pass = isGoodJet_Sync1(i);
    else  if (theCutScheme_==kBaseline0) pass = isGoodJet(i);
    else {assert(0);}

    if (pass ) {
      ngood++;
      if (ngood==n) return  loosejetEta->at(i);
    }
  }
  return 0;
}

unsigned int basicLoop::nGoodJets_Sync1() {
  
  unsigned int njets=0;

  for (unsigned int i=0; i<loosejetPt->size(); i++) {
    
    if (isGoodJet_Sync1(i) )   njets++;
  }
  return njets;
}

unsigned int basicLoop::nGoodJets() {
  
  unsigned int njets=0;

  for (unsigned int i=0; i<loosejetPt->size(); i++) {
    
    if (isGoodJet(i) )   njets++;
  }
  return njets;
}


int basicLoop::getTopDecayCategory() {
  /*
    this is a kludge on top of a kludge, I suppose.
    reclassifying one set of arbitrary integers with another set of arbitrary integers.
    ugly, but i don't have a better idea.
    i will use enums to make things a little more clear
  */

  //here are the codes defined in the ntuple
  enum TopDecayCodes {kUnknown=-1, kNoB = 0, kHadronic = 1, kElectron = 2, kMuon=3,kTauHadronic=4,kTauElectron=5,kTauMuon=6,kTauMisc=7 };
  //  TopDecayCodes mycode=kUnknown;

  //the set of codes returned here is defined in the class defn at the top of this file

  int code=-1;

  unsigned int ntop=  topDecayCode->size();

  if (ntop==2) { //this is what we expect for ttbar
    int code1 = topDecayCode->at(0);
    int code2 = topDecayCode->at(1);

    //i'm pretty sure that the 'misc' category is only filled by hadronic decays
    if (code1==kTauMisc) code1=kHadronic;
    if (code2==kTauMisc) code2=kHadronic;
    
    if ( (code1==kElectron ||code1==kMuon) && (code2==kElectron ||code2==kMuon) ) code=kAllLeptons;
    else if ( code1==kHadronic && code2==kHadronic ) code=kAllHadronic;
    else if ( (code1==kElectron && code2==kHadronic) || (code2==kElectron && code1==kHadronic) ) code=kOneElectron;
    else if ( (code1==kMuon && code2==kHadronic) || (code2==kMuon && code1==kHadronic) ) code=kOneMuon;
    else if ( (code1==kTauElectron && code2==kHadronic) || (code2==kTauElectron && code1==kHadronic) ) code=kOneTauE;
    else if ( (code1==kTauMuon && code2==kHadronic) || (code2==kTauMuon && code1==kHadronic) ) code=kOneTauMu;
    else if ( (code1==kTauHadronic && code2==kHadronic) || (code2==kTauHadronic && code1==kHadronic) ) code=kOneTauHadronic;
    //this logic depends on all of the highest categories being tau
    else if ( code1>=kTauHadronic && code2>=kTauHadronic ) code=kAllTau;
    else if ( code1>=kTauHadronic && (code2==kElectron || code2==kMuon) ) code=kTauPlusLepton;
    else if ( code2>=kTauHadronic && (code1==kElectron || code1==kMuon) ) code=kTauPlusLepton;
    else {
      code=kTTbarUnknown;
      std::cout<<"TopCode = "<<code1<<"\t"<<code2<<std::endl;
    }
  }
  //  else if (ntop==1) { //maybe we expect this for single top?
  //
  //  }
  //  else {
  //    
  //  }

  return code;
}

void basicLoop::fillTightJetInfo() {

  //first clear old values
  jetPt.clear();
  jetEta.clear();
  jetPhi.clear();
  jetFlavor.clear();
  jetBTagDisc_trackCountingHighPurBJetTags.clear();
  jetBTagDisc_trackCountingHighEffBJetTags.clear();
  jetBTagDisc_simpleSecondaryVertexHighEffBJetTags.clear();
  jetBTagDisc_simpleSecondaryVertexHighPurBJetTags.clear();
  jetBTagDisc_simpleSecondaryVertexBJetTags.clear();

  //now fill with new values
  for (unsigned int i=0; i<tightJetIndex->size(); i++) {
    int j = tightJetIndex->at(i);

    jetPt.push_back( loosejetPt->at( j ) );
    jetEta.push_back( loosejetEta->at( j ) );
    jetPhi.push_back( loosejetPhi->at( j ) );
    jetFlavor.push_back( loosejetFlavor->at( j ) );

    jetBTagDisc_trackCountingHighPurBJetTags.push_back( loosejetBTagDisc_trackCountingHighPurBJetTags->at( j ) );
    jetBTagDisc_trackCountingHighEffBJetTags.push_back( loosejetBTagDisc_trackCountingHighEffBJetTags->at( j ) );

    jetBTagDisc_simpleSecondaryVertexHighPurBJetTags.push_back( loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags->at( j ) );
    jetBTagDisc_simpleSecondaryVertexHighEffBJetTags.push_back( loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags->at( j ) );

    jetBTagDisc_simpleSecondaryVertexBJetTags.push_back( loosejetBTagDisc_simpleSecondaryVertexBJetTags->at( j ) );
  }

  //recover the nbSSVM variable (filled incorrectly in the ntuple)
  //  cout<<nbSSVM<<" ";
  nbSSVM=0;
  if (theCutScheme_==kBaseline0) nbSSVM = countBJets();
  else {
    for ( unsigned int i=0; i<jetPt.size(); i++) {
      if (passBCut(i)) nbSSVM++;
    }
  }

}

//when new jet variables are added to the ntuple then need to update this!
void basicLoop::InitJets() {

  if (theJetType_==kCalo ) {
    tightJetIndex = tightJetIndex_calo;
    looseJetIndex = looseJetIndex_calo;
    loosejetPt = loosejetPt_calo;
    loosejetEt = loosejetEt_calo;
    loosejetEta = loosejetEta_calo;
    loosejetPhi = loosejetPhi_calo;
    loosejetPassLooseID = loosejetPassLooseID_calo;
    loosejetEnergyFracHadronic = loosejetEnergyFracHadronic_calo;
    loosejetFlavor = loosejetFlavor_calo;
    loosejetGenParticlePDGId = loosejetGenParticlePDGId_calo;
    loosejetInvisibleEnergy = loosejetInvisibleEnergy_calo;
    loosejetBTagDisc_trackCountingHighPurBJetTags = loosejetBTagDisc_trackCountingHighPurBJetTags_calo;
    loosejetBTagDisc_trackCountingHighEffBJetTags = loosejetBTagDisc_trackCountingHighEffBJetTags_calo;
    loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags = loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags_calo;
    loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags = loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags_calo;
    loosejetBTagDisc_simpleSecondaryVertexBJetTags = loosejetBTagDisc_simpleSecondaryVertexBJetTags_calo;
    badjetPt = badjetPt_calo;
    badjetEta = badjetEta_calo;
    badjetPhi = badjetPhi_calo;

  //new in V00-01-02
    //this is ok in any version (the pointer is 0 in older ntuple versions)
    loosejetGenPt = loosejetGenPt_calo;
    loosejetGenEta = loosejetGenEta_calo;
    loosejetGenPhi = loosejetGenPhi_calo;
    loosejetNSV = loosejetNSV_calo;
    loosejetNTracks = loosejetNTracks_calo;
    loosejetSVUnWeightedMass = loosejetSVUnWeightedMass_calo;
    loosejetSVWeightedMass = loosejetSVWeightedMass_calo;
    loosejetSVUnWeightedLifetime = loosejetSVUnWeightedLifetime_calo;
    loosejetSVWeightedLifetime = loosejetSVWeightedLifetime_calo;
    loosejetSVUnWeightedCosTheta = loosejetSVUnWeightedCosTheta_calo;
    loosejetSVWeightedCosTheta = loosejetSVWeightedCosTheta_calo;
  }
  else if (theJetType_==kPF) {
    tightJetIndex = tightJetIndex_PF;
    looseJetIndex = looseJetIndex_PF;
    loosejetPt = loosejetPt_PF;
    loosejetEt = loosejetEt_PF;
    loosejetEta = loosejetEta_PF;
    loosejetPhi = loosejetPhi_PF;
    loosejetPassLooseID = loosejetPassLooseID_PF;
    loosejetEnergyFracHadronic = loosejetEnergyFracHadronic_PF;
    loosejetFlavor = loosejetFlavor_PF;
    loosejetGenParticlePDGId = loosejetGenParticlePDGId_PF;
    loosejetInvisibleEnergy = loosejetInvisibleEnergy_PF;
    loosejetBTagDisc_trackCountingHighPurBJetTags = loosejetBTagDisc_trackCountingHighPurBJetTags_PF;
    loosejetBTagDisc_trackCountingHighEffBJetTags = loosejetBTagDisc_trackCountingHighEffBJetTags_PF;
    loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags = loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags_PF;
    loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags = loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags_PF;
    loosejetBTagDisc_simpleSecondaryVertexBJetTags = loosejetBTagDisc_simpleSecondaryVertexBJetTags_PF;
    badjetPt = badjetPt_PF;
    badjetEta = badjetEta_PF;
    badjetPhi = badjetPhi_PF;
  //new in V00-01-02
    loosejetGenPt = loosejetGenPt_PF;
    loosejetGenEta = loosejetGenEta_PF;
    loosejetGenPhi = loosejetGenPhi_PF;
    loosejetNSV = loosejetNSV_PF;
    loosejetNTracks = loosejetNTracks_PF;
    loosejetSVUnWeightedMass = loosejetSVUnWeightedMass_PF;
    loosejetSVWeightedMass = loosejetSVWeightedMass_PF;
    loosejetSVUnWeightedLifetime = loosejetSVUnWeightedLifetime_PF;
    loosejetSVWeightedLifetime = loosejetSVWeightedLifetime_PF;
    loosejetSVUnWeightedCosTheta = loosejetSVUnWeightedCosTheta_PF;
    loosejetSVWeightedCosTheta = loosejetSVWeightedCosTheta_PF;
  }
  else {assert(0);}

}


TString basicLoop::getSampleName(TString inname) {
  if (inname=="") inname=findInputName();

  if (inname.Contains("/TTbarJets/") )                     return "TTbarJets";
  else if (inname.Contains("/LM0/"))                       return "LM0";
  else if (inname.Contains("/LM1/"))                       return "LM1";
  else if (inname.Contains("/LM2/"))                       return "LM2";
  else if (inname.Contains("/LM3/"))                       return "LM3";
  else if (inname.Contains("/LM4/"))                       return "LM4";
  else if (inname.Contains("/LM5/"))                       return "LM5";
  else if (inname.Contains("/LM6/"))                       return "LM6";
  else if (inname.Contains("/LM7/"))                       return "LM7";
  else if (inname.Contains("/LM8/"))                       return "LM8";
  else if (inname.Contains("/LM9/"))                       return "LM9";
  else if (inname.Contains("/LM9t175/"))                   return "LM9t175";
  else if (inname.Contains("/LM9p/"))                      return "LM9p";
  else if (inname.Contains("/LM10/"))                      return "LM10";
  else if (inname.Contains("/LM11/"))                      return "LM11";
  else if (inname.Contains("/LM12/"))                      return "LM12";
  else if (inname.Contains("/LM13/"))                      return "LM13";
  else if (inname.Contains("/MoreMSSM/"))                  return "mMSSM";
  else if (inname.Contains("/MoreMSSMv2/"))                return "mMSSMv2";
  else if (inname.Contains("/MoreMSSMv3/"))                return "mMSSMv3";
  else if (inname.Contains("/SingleTop-sChannel/"))        return "SingleTop-sChannel";
  else if (inname.Contains("/SingleTop-tChannel/"))        return "SingleTop-tChannel";
  else if (inname.Contains("/SingleTop-tWChannel/"))       return "SingleTop-tWChannel";
  else if (inname.Contains("/QCD-Pt1000toInf-madgraph/"))  return "QCD1000";
  else if (inname.Contains("/QCD-Pt100to250-madgraph/"))   return "QCD100";
  else if (inname.Contains("/QCD-Pt250to500-madgraph/"))   return "QCD250";
  else if (inname.Contains("/QCD-Pt500to1000-madgraph/"))  return "QCD500";
  else if (inname.Contains("/WJets/"))                     return "WJets";
  else if (inname.Contains("/ZJets/"))                     return "ZJets";
  else if (inname.Contains("/Zinvisible/"))                return "Zinvisible";

  else if (inname.Contains("/DATA/"))                      return "data";
  
  std::cout<<"Cannot find sample name for this sample!"<<std::endl;
  
  return "";

}

double basicLoop::getCrossSection( TString inname) {
  if (inname=="") inname=findInputName();

  //  https://twiki.cern.ch/twiki/bin/view/CMS/ProductionReProcessingSpring10
  const double bf = 0.32442;

  //those marked with !!! in the comment have been synchronized with our note
  if (inname.Contains("/TTbarJets/") )                     return 165; //!!! NNLO
  
  else if (inname.Contains("/LM0/"))                       return 38.93;
  else if (inname.Contains("/LM1/"))                       return 4.888;
  else if (inname.Contains("/LM2/"))                       return 0.6027;
  else if (inname.Contains("/LM3/"))                       return 3.438;
  else if (inname.Contains("/LM4/"))                       return 1.879;
  else if (inname.Contains("/LM5/"))                       return 0.4734;
  else if (inname.Contains("/LM6/"))                       return 0.3104;
  else if (inname.Contains("/LM7/"))                       return 1.209;
  else if (inname.Contains("/LM8/"))                       return 0.7300;
  else if (inname.Contains("/LM9/"))                       return 7.134; //!!! LO
  else if (inname.Contains("/LM9t175/"))                   return 4.241;
  else if (inname.Contains("/LM9p/"))                      return 1.653;
  else if (inname.Contains("/LM10/"))                      return 0.04778;
  else if (inname.Contains("/LM11/"))                      return 0.8236;
  else if (inname.Contains("/LM12/"))                      return 4.414;
  else if (inname.Contains("/LM13/"))                      return 6.899; //!!! LO
  
  //not updated from tom
  else if (inname.Contains("/LM9t175/"))                   return 4.241;
  else if (inname.Contains("/LM9p/"))                      return 1.653;

  ///these are from Tom
  /*
  else if (inname.Contains("/LM0/"))                       return 54.89;
  else if (inname.Contains("/LM1/"))                       return 6.550;
  else if (inname.Contains("/LM2/"))                       return 0.8015;
  else if (inname.Contains("/LM3/"))                       return 4.813;
  else if (inname.Contains("/LM4/"))                       return 2.537;
  else if (inname.Contains("/LM5/"))                       return 0.634;
  else if (inname.Contains("/LM6/"))                       return 0.4035;
  else if (inname.Contains("/LM7/"))                       return 1.342;
  else if (inname.Contains("/LM8/"))                       return 1.0293;
  else if (inname.Contains("/LM9/"))                       return 10.558;
  else if (inname.Contains("/LM10/"))                      return 0.04778;
  else if (inname.Contains("/LM11/"))                      return 1.1119;
  else if (inname.Contains("/LM12/"))                      return 5.915;
  else if (inname.Contains("/LM13/"))                      return 9.797;
  */

  else if (inname.Contains("/MoreMSSM/"))                  return 1.73;
  else if (inname.Contains("/MoreMSSMv2/"))                return 2.1;
  else if (inname.Contains("/MoreMSSMv3/"))                return 2.6;
  else if (inname.Contains("/SingleTop-sChannel/"))        return bf*4.6; //!!! NNNLO
  else if (inname.Contains("/SingleTop-tChannel/"))        return bf*64.6; //!!! NLO
  else if (inname.Contains("/SingleTop-tWChannel/"))       return 10.6; //!!! NLO
  //QCD are all LO
  else if (inname.Contains("/QCD-Pt100to250-madgraph/"))   return 7e6; //!!!
  else if (inname.Contains("/QCD-Pt250to500-madgraph/"))   return 171000; //!!!
  else if (inname.Contains("/QCD-Pt500to1000-madgraph/"))  return 5200; //!!!
  else if (inname.Contains("/QCD-Pt1000toInf-madgraph/"))  return 83; //!!!

  else if (inname.Contains("/WJets/"))                     return 31314; //!!! NNLO

  //also known as DYJetsToLL with m_ll > 50
  else if (inname.Contains("/ZJets/"))                     return 3048; //!!! NNLO
  //i don't think we've processed DY->ll (10<m_ll<50)

  else if (inname.Contains("/Zinvisible/"))                return 4500; //!!! LO
  //tom had 4900...this might be NLO
  
  else if (inname.Contains("/DATA/"))                return -2;

  std::cout<<"Cannot find cross section for this sample!"<<std::endl; 
  return -1;
  
}

TString basicLoop::findInputName() {
   // ===== do something tricky to figure out which sample this is ======
  if (fChain==0) return "errorNoChain";
   TString inname="";
   if (fChain->InheritsFrom(TChain::Class())) {
     if (((TChain*) fChain)->GetListOfFiles()->GetEntries() <1) {
       std::cout<<"Chain seems to have no files!"<<std::endl; return "error";
     }
     inname = ((TChain*) fChain)->GetListOfFiles()->At(0)->GetTitle();
   }
   else {
     inname=  fChain->GetCurrentFile()->GetName();
   }
   //this is the file name of the first file in the chain

   return inname;
}

void basicLoop::setMETType(METType mettype) {

  theMETType_ = mettype;
}

void basicLoop::setJetType(jetType jettype) {

  theJetType_ = jettype;

}

void basicLoop::setLeptonType(leptonType leptontype) {

  theLeptonType_ = leptontype;

}

void basicLoop::setEleReq(int ne) {
  ne_ = ne;
}

void basicLoop::setMuonReq(int nmu) {
  nmu_ = nmu;
}

void basicLoop::setMETRange(METRange metrange) {

  theMETRange_ = metrange;
}

void basicLoop::setDPType(dpType dptype) {

  theDPType_ = dptype;

}

void basicLoop::setBCut(unsigned int nb) {
  nBcut_=nb;
}

void basicLoop::setIgnoredCut(const TString cutTag) {

  bool ok=false;
  for (unsigned int i=0; i<cutTags_.size() ; i++) {
    if (cutTags_[i] == cutTag) {ok=true; break;}
  }
  if (!ok) {cout<<"Invalid cutIndex"<<endl; return;}

  ignoredCut_.push_back(cutTag);

}

void basicLoop::resetIgnoredCut() {
  ignoredCut_.clear();
}

void basicLoop::specifyEvent(ULong64_t run, ULong64_t lumisection, ULong64_t event) {

  jmt::eventID myevent;
  myevent.run = run;
  myevent.ls=lumisection;
  myevent.ev = event;
  specifiedEvents_.insert( myevent);

}

bool basicLoop::eventIsSpecified() {
  //check if the current event is on the specified list
  //if the list is empty, just return true

  //this implementation does not allow for a wildcard-type search (e.g. all of one lumi section)
  //can think about it....

  if (specifiedEvents_.empty()) return true;

  jmt::eventID thisevent;
  thisevent.run = runNumber;
  thisevent.ls = lumiSection;
  thisevent.ev = eventNumber;

  if ( specifiedEvents_.find( thisevent) != specifiedEvents_.end()) return true;

  return false;
}


TString basicLoop::getCutDescriptionString() {

  TString cuts = CutSchemeNames_[theCutScheme_];
  cuts += "_";
  cuts += jetTypeNames_[theJetType_];
  cuts += "_";
  cuts += METTypeNames_[theMETType_];
  cuts += METRangeNames_[theMETRange_];
  cuts += "_";
  cuts += leptonTypeNames_[theLeptonType_];
  cuts += ne_; 
  cuts += "e";
  cuts += nmu_; 
  cuts += "mu";
  cuts += "_";
  cuts+= dpTypeNames_[theDPType_];
  for (unsigned int icut=0; icut<ignoredCut_.size() ; icut++) {
    cuts+="_No";
    //it would be more robust to use .find() instead of []
    cuts += jmt::fortranize( cutNames_[ignoredCut_.at(icut)]);
  }
  return cuts; 
}

TString basicLoop::getBCutDescriptionString() {

  TString thecut="ge";
  thecut += nBcut_;
  thecut+="b";
  return thecut;
}


void basicLoop::printState() {
  cout<<"Weights are for L = "<<lumi<<" pb^-1"<<endl;
  cout<<"Cut scheme set to:    "<<CutSchemeNames_[theCutScheme_]<<endl;
  cout<<"jet type set to:      "<<jetTypeNames_[theJetType_]<<endl;
  cout<<"lepton type set to:   "<<leptonTypeNames_[theLeptonType_]<<endl;
  cout<<"Requiring exactly     "<<ne_<<" electrons"<<endl;
  cout<<"Requiring exactly     "<<nmu_<<" muons"<<endl;
  cout<<"MET type set to:      "<<METTypeNames_[theMETType_]<<endl;
  cout<<"MET range set to:     "<<METRangeNames_[theMETRange_]<<endl;
  cout<<"DeltaPhi type set to: "<<dpTypeNames_[theDPType_]<<endl;
  cout<<"Requiring at least    "<<nBcut_<<" b tags"<<endl;
  for (unsigned int i = 0; i< ignoredCut_.size() ; i++) {
    cout<<"Will ignore cut:    "<<cutNames_[ignoredCut_.at(i)]<<endl;
  }
  
}

// ========================================== end


#endif // #ifdef basicLoop_cxx
