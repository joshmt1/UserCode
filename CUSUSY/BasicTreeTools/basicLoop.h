//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar 21 14:11:02 2011 by ROOT version 5.22/00d
// from TTree tree/tree
// found on file: /cu1/joshmt/BasicNtuples/V00-03-01/LM13/BasicNtuple.root
//////////////////////////////////////////////////////////

#ifndef basicLoop_h
#define basicLoop_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
// ========================================== begin
#include <TDatime.h>
#include <TH1.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "TMatrixT.h"
#include "TMatrixDEigen.h"
//this file is in CVS here: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/joshmt/MiscUtil.cxx?view=log
#include "MiscUtil.cxx"
//this code will have to be regenerated when changing the ntuple structure
//custom code is marked with these 'begin' and 'end' markers
// ---- this version is compatible with ntuple tag: V00-03-01 ------
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
RA2 -- the legacy scheme. no longer trustworthy
Sync1 -- cuts for first sync exercise. Recomputes (almost?) every cut.
Baseline0 -- *must* recompute every cut (except the PV, which I take from the ntuple). "Baseline" analysis cuts as laid out by me on the joint Twiki.

Other notes:
a lot of the functions are still doing computations using the tight jet list stored in the ntuple.
This is rather dangerous. In particular, the eta range of the ntuple tight jets does not agree 100% with the kBaseline0 tight jets.

== work in progress ==

* how to update for JES uncertainties?

I now get jet pt directly:
  loosejetPt->at(i)
Remember that loosejetPt is a pointer to the actual ntuple vector loosejetPt_PF

The obvious way is to rescale each of these directly:
  jesfactor * loosejetPt->at(i)
But then every time the user wants the Pt, he has to remember this factor.

Another approach is, instead of making loosejetPt a pointer to the loosejetPt_PF vector, make a copy when loading the event:
  loosejetPt->push_back( jesfactor*loosejetPt_PF->at(i) )
This would be transparent to the user. But that is more CPU intensive (quantitatively I'm not sure...)

Another approach is to replace all calls to loosejetPt->at(i) with:
  getLooseJetPt(i)
where the function returns jesfactor * loosejetPt->at(i)
Then the user only has to remember this abstract function call. The user must remember not to 
access the loosjetPt vector directly.
----> I think this is the best approach!

-- plans for next major rewrite: --
I want to unify all functions such as getHT() so that they return the value determined by the main configuraton enums.
Then allow the user to override the default enums with arguments to the function. For instance,
getMET() will return the MET corresponding to theMETType_, but the user can always ask for:
getMET( kMHT) in order to get the MHT value.

This can ever be done with cut schemes, so that:
isGoodJet(ijet), for example, will change behavior based on cut scheme, but can be overridden with:
isGoodJet(ijet, aCutScheme, aJetType);

aCutScheme determines which cuts are applied and with what values. aJetType determines what type of jets to use (PF, calo, etc)

This will take some care to implement, so for now I leave it on the drawing board.

*/

//avoid spaces and funny characters!
const char *CutSchemeNames_[]={"RA2", "Sync1", "Baseline0"};
const char *METTypeNames_[]={"MHT", "MET",  "tcMET", "pfMET"};
const char *METRangeNames_[]={"med",  "high", "wide", "medhigh"}; //'no cut' is not given, because the cut can always be skipped!

const char *jetTypeNames_[]={"calo","PF", "JPT"};
const char *leptonTypeNames_[]={"RegLep","PFLep","PFLepRA2"};

const char *dpTypeNames_[]={"DeltaPhi", "minDP",  "MPT", "DPSync1", "minDPinv", "minDPAll30"};

const char *jesTypeNames_[] = {"JES0","JESup","JESdown"};
const char *jerTypeNames_[] = {"JER0","JERbias","JERup"}; //this one is weird -- 0==0==down, bias=0.1, up=0.2
const char *unclusteredMetUncNames_[] = {"METunc0","METdown","METup"};

const char *btagEffTypeNames_[] = {"BTagEff0","BTagEffup","BTagEffdown"};

const char *tailCleaningNames_[] = {"NoCleaning","MuonCleaning","MuonEcalCleaning"};

//in 1/pb
const double lumi=36.146; //386 Nov4ReReco Datasets - PATIFIED WITH 387

const double mW_ = 80.399;
const double mtop_ = 172.0;

// ========================================== end

class basicLoop {
public :
  // ========================================== begin
  //TO DO (sometime when there is time to validate it): change first item of each enum to be a different integer
  //i think this provides safety against the user passing the wrong type of enum to the wrong set method
  enum CutScheme {kRA2=0, kSync1, kBaseline0, nCutSchemes}; //cut schemes must be implemented in setCutScheme,etc and in the list above
  CutScheme theCutScheme_;
  enum METType {kMHT=0, kMET, ktcMET, kpfMET};
  METType theMETType_;
  enum METRange {kMedium=0, kHigh, kWide, kMedhigh};
  METRange theMETRange_;
  enum jetType {kCalo=0, kPF, kJPT};
  jetType theJetType_;
  enum leptonType {kNormal=0, kPFLeptons, kPFLeptonsRA2};
  leptonType theLeptonType_;
  enum dpType {kDeltaPhi=0, kminDP, kMPT, kDPSync1, kminDPinv, kminDPAll30};
  dpType theDPType_;
  enum JESType {kJES0=0,kJESup,kJESdown};
  JESType theJESType_;
  enum JERType {kJER0=0,kJERbias,kJERup};
  JERType theJERType_;
  enum METuncType {kMETunc0=0,kMETuncDown,kMETuncUp};
  METuncType theMETuncType_;
  enum BTagEffType {kBTagEff0=0,kBTagEffup,kBTagEffdown};
  BTagEffType theBTagEffType_;
  enum tailCleaningType {kNoCleaning=0, kMuonCleaning, kMuonEcalCleaning};
  tailCleaningType theCleaningType_;

  unsigned int nBcut_;
  //these are an extension of the ele/mu veto
  //If the cutEleVeto or cutMuVeto is set, then there must be exactly this number of e/mu in the event
  //by default, init to 0 (veto)
  //setting to 1 inverts the veto
  int ne_;
  int nmu_;

  bool isData_; //set in ctor

  bool  realDatasetNames_;

  //any trigger on this list will be considered OK for passing the trigger cut
  //key is the trigger name ; value is the position in the passTrigger variable in the ntuple
  std::map<TString, int>  triggerList_; //contains only triggers used for determining the signal
  std::map<TString, int>  triggerListAll_; //this contains all triggers
  //can use hltPrescale to figure out if the trigger exists in a given event/sample

  std::vector<TString> cutTags_;
  std::map<TString, TString> cutNames_; //key is a cutTag. these should be "human readable" but should not have any spaces
  std::vector<TString> ignoredCut_; //allow more than 1 ignored cut!
  std::vector<TString> requiredCut_; //new feature to *turn on* a cut that is usually not required by a given cut scheme
  //if theCutFlow changes, be sure to change cutnames_ as well

  std::set<jmt::eventID> specifiedEvents_;
  std::set<jmt::eventID> ecalVetoEvents_;
  bool loadedEcalTree_;

  enum TopDecayCategory {kTTbarUnknown=0,kAllLeptons=1,kAllHadronic=2,kOneElectron=3,kOneMuon=4,kOneTauE=5,kOneTauMu=6,kOneTauHadronic=7,kAllTau=8,kTauPlusLepton=9, nTopCategories=10};

  TDatime* starttime_;
  TString specialCutDescription_;

  //removing code that fills "tight jet" info for each event.

  //pointers to the jet info of the default jet type
  vector<int>     *tightJetIndex;
  vector<int>     *looseJetIndex;
  vector<float>   *loosejetPt;
  vector<float>   *loosejetPz;
  vector<float>   *loosejetPtUncorr;
  vector<float>   *loosejetEt;
  vector<float>   *loosejetE;
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

  //new in V00-02-xx
   vector<float>   *loosejetJECUncPlus;
   vector<float>   *loosejetJECUncMinus;

   Int_t nbSSVM;
   float WCosHel_,topCosHel_,bestWMass_,bestTopMass_;
   float eleet1_,muonpt1_; //FIXME this is a hack that I don't like!
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
   vector<bool>    *passTrigger;
   vector<unsigned int> *hltPrescale;
   Bool_t          pv_pass;
   vector<bool>    *pv_isFake;
   vector<float>   *pv_z;
   vector<float>   *pv_rho;
   vector<float>   *pv_chi2;
   vector<float>   *pv_ndof;
   vector<float>   *muonPtDiff_PF;
   Bool_t          passesBadPFMuonFilter;
   Bool_t          passesInconsistentMuonPFCandidateFilter;
   vector<int>     *tightJetIndex_calo;
   vector<int>     *looseJetIndex_calo;
   vector<float>   *loosejetPt_calo;
   vector<float>   *loosejetPz_calo;
   vector<float>   *loosejetPtUncorr_calo;
   vector<float>   *loosejetEt_calo;
   vector<float>   *loosejetE_calo;
   vector<float>   *loosejetEta_calo;
   vector<float>   *loosejetPhi_calo;
   vector<bool>    *loosejetPassLooseID_calo;
   vector<bool>    *loosejetPassTightID_calo;
   vector<float>   *loosejetEnergyFracHadronic_calo;
   vector<float>   *loosejetJECUncPlus_calo;
   vector<float>   *loosejetJECUncMinus_calo;
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
   vector<int>     *tightJetIndex_PF;
   vector<int>     *looseJetIndex_PF;
   vector<float>   *loosejetPt_PF;
   vector<float>   *loosejetPz_PF;
   vector<float>   *loosejetPtUncorr_PF;
   vector<float>   *loosejetEt_PF;
   vector<float>   *loosejetE_PF;
   vector<float>   *loosejetEta_PF;
   vector<float>   *loosejetPhi_PF;
   vector<bool>    *loosejetPassLooseID_PF;
   vector<bool>    *loosejetPassTightID_PF;
   vector<float>   *loosejetEnergyFracHadronic_PF;
   vector<float>   *loosejetJECUncPlus_PF;
   vector<float>   *loosejetJECUncMinus_PF;
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
   Float_t         pf1MET;
   Float_t         pf1METphi;
   Float_t         pf1METsig;
   Float_t         pf1GenMET;
   Float_t         pf1GenMETphi;
   vector<float>   *trackPt;
   vector<float>   *trackEta;
   vector<float>   *trackPhi;
   Float_t         SumPtOverHT;
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
   vector<bool>    *muonIsRA2;
   vector<bool>    *muonIsGlobalMuon;
   vector<bool>    *muonIsGlobalMuonPromptTight;
   vector<float>   *muonEoverP;
   vector<float>   *muonPtRatio;
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
   vector<bool>    *eleIsRA2;
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
   vector<bool>    *muonIsRA2_PF;
   vector<bool>    *muonIsGlobalMuon_PF;
   vector<bool>    *muonIsGlobalMuonPromptTight_PF;
   vector<float>   *muonEoverP_PF;
   vector<float>   *muonPtRatio_PF;
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
   vector<bool>    *eleIsRA2_PF;
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
   Int_t           SUSY_nb;
   Double_t        qScale;
   Double_t        ptHat;
   Double_t        mcWeight;
   vector<float>   *pdfWeights;
   vector<int>     *topDecayCode;
   Int_t           ZDecayMode;
   Int_t           flavorHistory;

   // List of branches
   TBranch        *b_runNumber;   //!
   TBranch        *b_lumiSection;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_bsx;   //!
   TBranch        *b_bsy;   //!
   TBranch        *b_bsz;   //!
   TBranch        *b_passTrigger;   //!
   TBranch        *b_hltPrescale;   //!
   TBranch        *b_pv_pass;   //!
   TBranch        *b_pv_isFake;   //!
   TBranch        *b_pv_z;   //!
   TBranch        *b_pv_rho;   //!
   TBranch        *b_pv_chi2;   //!
   TBranch        *b_pv_ndof;   //!
   TBranch        *b_muonPtDiff_PF;   //!
   TBranch        *b_passesBadPFMuonFilter;   //!
   TBranch        *b_passesInconsistentMuonPFCandidateFilter;   //!
   TBranch        *b_tightJetIndex_calo;   //!
   TBranch        *b_looseJetIndex_calo;   //!
   TBranch        *b_loosejetPt_calo;   //!
   TBranch        *b_loosejetPz_calo;   //!
   TBranch        *b_loosejetPtUncorr_calo;   //!
   TBranch        *b_loosejetEt_calo;   //!
   TBranch        *b_loosejetE_calo;   //!
   TBranch        *b_loosejetEta_calo;   //!
   TBranch        *b_loosejetPhi_calo;   //!
   TBranch        *b_loosejetPassLooseID_calo;   //!
   TBranch        *b_loosejetPassTightID_calo;   //!
   TBranch        *b_loosejetEnergyFracHadronic_calo;   //!
   TBranch        *b_loosejetJECUncPlus_calo;   //!
   TBranch        *b_loosejetJECUncMinus_calo;   //!
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
   TBranch        *b_tightJetIndex_PF;   //!
   TBranch        *b_looseJetIndex_PF;   //!
   TBranch        *b_loosejetPt_PF;   //!
   TBranch        *b_loosejetPz_PF;   //!
   TBranch        *b_loosejetPtUncorr_PF;   //!
   TBranch        *b_loosejetEt_PF;   //!
   TBranch        *b_loosejetE_PF;   //!
   TBranch        *b_loosejetEta_PF;   //!
   TBranch        *b_loosejetPhi_PF;   //!
   TBranch        *b_loosejetPassLooseID_PF;   //!
   TBranch        *b_loosejetPassTightID_PF;   //!
   TBranch        *b_loosejetEnergyFracHadronic_PF;   //!
   TBranch        *b_loosejetJECUncPlus_PF;   //!
   TBranch        *b_loosejetJECUncMinus_PF;   //!
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
   TBranch        *b_pf1MET;   //!
   TBranch        *b_pf1METphi;   //!
   TBranch        *b_pf1METsig;   //!
   TBranch        *b_pf1GenMET;   //!
   TBranch        *b_pf1GenMETphi;   //!
   TBranch        *b_trackPt;   //!
   TBranch        *b_trackEta;   //!
   TBranch        *b_trackPhi;   //!
   TBranch        *b_SumPtOverHT;   //!
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
   TBranch        *b_muonIsRA2;   //!
   TBranch        *b_muonIsGlobalMuon;   //!
   TBranch        *b_muonIsGlobalMuonPromptTight;   //!
   TBranch        *b_muonEoverP;   //!
   TBranch        *b_muonPtRatio;   //!
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
   TBranch        *b_eleIsRA2;   //!
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
   TBranch        *b_muonIsRA2_PF;   //!
   TBranch        *b_muonIsGlobalMuon_PF;   //!
   TBranch        *b_muonIsGlobalMuonPromptTight_PF;   //!
   TBranch        *b_muonEoverP_PF;   //!
   TBranch        *b_muonPtRatio_PF;   //!
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
   TBranch        *b_eleIsRA2_PF;   //!
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
   TBranch        *b_SUSY_nb;   //!
   TBranch        *b_qScale;   //!
   TBranch        *b_ptHat;   //!
   TBranch        *b_mcWeight;   //!
   TBranch        *b_pdfWeights;   //!
   TBranch        *b_topDecayCode;   //!
   TBranch        *b_ZDecayMode;   //!
   TBranch        *b_flavorHistory;   //!

   basicLoop(TTree *tree=0, TTree *infotree=0, TTree *ecaltree=0);    // ========================================== begin, end
   virtual ~basicLoop();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   // ========================================== begin
   //first the functions in basicLoop.C
   //these are "User" functions that loop over events.
   //the operation of the code should not depend on them
   virtual void     Loop(unsigned int dataindex=0);
   virtual void     exampleLoop();
   virtual void     screendump();
   virtual void     nbLoop();
   virtual void     ABCDtree(unsigned int dataindex=0);
   virtual void reducedTree(TString outputpath);
   virtual void triggerPlot();
   virtual void triggerPlotData();
   virtual void triggerTest();
   virtual void Nminus1plots();
   virtual void calculateTagProb();
   void cutflow(bool writeFiles=false);
   void cutflowPlotter();

   //below here are == functions in basicLoop.h ==
   //these are utilities that e.g. calculate useful quantities for a given event

   //some stuff that is used internally
   void fillTightJetInfo();
   void InitJets();
   void fillEcalVetoList(TTree* ecaltree);

   //performance timing
   void startTimer();
   void checkTimer(const Long64_t ndone, const Long64_t ntotal);
   void stopTimer(const Long64_t ntotal);

   //configuration options
   void printState(); //dump useful configuration info to the screen
   bool setCutScheme(CutScheme cutscheme);
   void setMETType(METType mettype);
   void setJetType(jetType jettype);
   void setLeptonType(leptonType leptontype);
   void setMETRange(METRange metrange);
   void setMETuncType(METuncType metunctype);
   void setDPType(dpType dptype);
   void setJESType(JESType jestype);
   void setJERType(JERType jertype);
   void setCleaningType(tailCleaningType cleanuptype);
   void setIgnoredCut(const TString cutTag);
   void resetIgnoredCut() ;
   void setRequiredCut(const TString cutTag);
   void resetRequiredCut() ;
   void setBCut(unsigned int nb);
   void setEleReq(int ne); //change Cut suffix to Req because it is an == requirement, not >=
   void setMuonReq(int nmu);
   TString getCutDescriptionString();
   TString getBCutDescriptionString();

   void setBTagEffType(BTagEffType btagefftype);

   //really special configuration options (for expert use)
   void specifyEvent(ULong64_t run, ULong64_t lumisection, ULong64_t event);
   bool eventIsSpecified();
   void setSpecialCutDescription(TString cutDesc) {specialCutDescription_=cutDesc;}
   void useRealDatasetNames(bool usereal) { realDatasetNames_=usereal;}

   //calculate useful stuff and get useful info
   TString getSampleName(TString inname="") ;
   double getWeight(Long64_t nentries);
   double getCrossSection(TString inname="") ;
   TString findInputName() ;

   int countGoodPV() ;

   float getMET(); //return MET determined by theMETType_
   float getMETphi(); //return MET determined by theMETType_
   float getGenMET(); //return MET determined by theMETType_
   float getGenMETphi(); //return MET determined by theMETType_
   int countGenBJets(float threshold);  //count gen b jets using loosejetFlavor
   float getJetInvisibleEnergyHT(); //could add a pT cut as an argument
   float getJetInvisibleEnergyMHT(); //could add a pT cut as an argument
   float getLargestJetPtRecoError(unsigned int maxjets);
   float getDeltaPhiMismeasuredMET(unsigned int maxjets);

   //for studying mc truth jet flavor (hiding the loop over jets from the user)
   void fillWithJetFlavor(TH1D* hh, double w, double threshold);
   void fillWithGenPDGId(TH1D* hh, double w, double threshold);

   //for systematics (return is a pair with METx, METy)
   std::pair<float, float> getJESAdjustedMETxy();
   std::pair<float, float> getJERAdjustedMETxy();
   std::pair<float, float> getJERAdjustedMHTxy();
   std::pair<float,float> getUnclusteredSmearedMETxy() ;
   float getJESExtraUnc(unsigned int ijet);

   bool cutRequired(TString cutTag) ;
   bool passCut(TString cutTag) ;

   //used for Sync1 and kBaseline0
   bool passMuVeto() ;
   bool passEleVeto() ;
   int countMu() ;
   int countEle() ;

   bool passTauVeto();
   bool passCleaning();
   bool passEcalDeadCellCleaning();
   bool passGreedyMuon();
   float findMaxMuonPtDiff();
   void findMostInconsistentPFMuon(float &mostInconsistentMuonPt, float &mostInconsistentMuonPtDiff);

   bool passSSVM(int i); //index is for loose jets

   //note that !(bad jet) != good jet
   //a bad jet is a high pT jet with too wide eta or failing jet ID
   //a good jet is high pT that has good eta and passes jet ID
   //many soft jets are neither good jets nor bad jets
   bool isBadJet(unsigned int ijet) ; //index on loose jet list
   unsigned nBadJets();

   bool passBadJetVeto(); //n bad jets == 0 or not

   float getLooseJetPt(unsigned int ijet); //includes JES
   float getLooseJetPx(unsigned int ijet);
   float getLooseJetPy(unsigned int ijet);
   float getLooseJetPtUncorr(unsigned int ijet); //includes JER
   bool isGoodJet_Sync1(unsigned int ijet);
   bool isGoodJet(unsigned int ijet); //index here is on the loose jet list
   bool isGoodJet10(unsigned int ijet); //index here is on the loose jet list
   bool isGoodJet30(unsigned int ijet); //index here is on the loose jet list
   bool isGoodJetMHT(unsigned int ijet); //index here is on the loose jet list
   bool isLooseJet_Sync1(unsigned int ijet); //looser pt cut

   float jetBTagEff(unsigned int ijet); //index here is on the loose jet list

   unsigned int nGoodJets_Sync1();
   unsigned int nGoodJets();
   float jetPtOfN(unsigned int n);
   float jetPhiOfN(unsigned int n);
   float jetEtaOfN(unsigned int n);
   float bjetPtOfN(unsigned int n);
   float bjetPhiOfN(unsigned int n);
   float bjetEtaOfN(unsigned int n);
   int countBJets_Sync1();
   int countBJets();
   bool passDeltaPhi_Sync1() ;
   float getHT_Sync1();
   float getHT();
   float getMHT();
   float getMHTphi();

   float getDeltaPhiTopTwoJets();

   bool passPV();

   bool passHLT(); bool printedHLT_; 
   TString lastTriggerPass_; //not really filled...but keep it in so that basicLoop.C compiles
   bool passHLT(const TString & triggerName);
   int getHLTPrescale(const TString & triggerName);

   void lookForPrescalePass();

   //for calculating the mass, cosHel of jet combinations
   void fillWTop(); //first is W mass, second it top mass
   double calc_mNj( std::vector<unsigned int> );
   double calc_mNj( unsigned int j1i, unsigned int j2i);
   double calc_mNj( unsigned int j1i, unsigned int j2i, unsigned int j3i);
   void calcCosHel( unsigned int j1i, unsigned int j2i, unsigned int j3i) ; //top, W cosHel

   //event shape variables from Luke
   void getSphericityJetMET(double & lambda1, double & lambda2, double & det,
			    const int jetmax, bool addMET);

   void getSphericityJetMET(float & lambda1, float & lambda2, float & det,
			    const int jetmax, bool addMET);

   double getDeltaPhiMPTMET();
   double getMinDeltaPhibMET() ;
   double getMinDeltaPhiMET(unsigned int maxjets) ;
   double getMinDeltaPhiMET30(unsigned int maxjets) ;
   double getMinDeltaPhiMET30_eta5(unsigned int maxjets) ;
   double getMinDeltaPhiMET30_eta5_noId(unsigned int maxjets) ;
   double getMaxDeltaPhiMET30(unsigned int maxjets) ;
   //   double getMinDeltaPhiMHT(unsigned int maxjets) ; //deprecated because MHT is now just another type of MET
   double getDeltaPhib1b2();
   double getUncorrectedHT(const double threshold);
   int getTopDecayCategory();
   double getMinDeltaPhi_bj(unsigned int bindex);
   //   double getOverallMinDeltaR_bj();
   //   double getMinDeltaR_bj(unsigned int bindex);
   double getDeltaPhi(double phi1, double phi2);
   // ========================================== end
};

#endif

#ifdef basicLoop_cxx
//====================== begin
//ecaltree is a tree of events failing the ECAL Dead Cell Filter
basicLoop::basicLoop(TTree *tree, TTree *infotree, TTree *ecaltree)
  :  theCutScheme_(kRA2),
     theMETType_(kMET),
     theMETRange_(kHigh),
     theJetType_(kCalo),
     theLeptonType_(kNormal),
     theDPType_(kminDP),
     theJESType_(kJES0),
     theJERType_(kJER0),
     theMETuncType_(kMETunc0),
     theCleaningType_(kNoCleaning),
     nBcut_(0),
     ne_(0),
     nmu_(0),
     isData_(false),
     realDatasetNames_(false),
     loadedEcalTree_(false),
     starttime_(0),
     specialCutDescription_(""),
     nbSSVM(0),
     WCosHel_(-99),
     topCosHel_(-99),
     bestWMass_(1e9),
     bestTopMass_(1e9),
     eleet1_(-1), muonpt1_(-1),
     printedHLT_(false),
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
   if (ecaltree != 0) fillEcalVetoList(ecaltree); else ecalVetoEvents_.clear();
   specifiedEvents_.clear();

   triggerList_.clear();
   triggerListAll_.clear();
   if (infotree!=0) {
     std::set<TString> triggersForCut;
     triggersForCut.insert("HLT_HT100U");
     //    triggersForCut.insert("HLT_HT120U");
     triggersForCut.insert("HLT_HT140U");
     //     triggersForCut.insert("HLT_HT150U");
     triggersForCut.insert("HLT_HT150U_v3");
     triggersForCut.insert("HLT_HT200");  //for old samples
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

   if (  (findInputName()).Contains("/V00-01-")  || (findInputName()).Contains("/V00-02-00") || (findInputName()).Contains("/V00-02-01")|| (findInputName()).Contains("/V00-02-02") ) {
     std::cout<<"Sorry, I am only compatible with V00-02-03 ntuples!"<<std::endl;
     assert(0);
   }
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

   WCosHel_=-99;
   topCosHel_=-99;
   bestWMass_=1e9;
   bestTopMass_=1e9;
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
   passTrigger = 0;
   hltPrescale = 0;
   pv_isFake = 0;
   pv_z = 0;
   pv_rho = 0;
   pv_chi2 = 0;
   pv_ndof = 0;
   muonPtDiff_PF = 0;
   tightJetIndex_calo = 0;
   looseJetIndex_calo = 0;
   loosejetPt_calo = 0;
   loosejetPz_calo = 0;
   loosejetPtUncorr_calo = 0;
   loosejetEt_calo = 0;
   loosejetE_calo = 0;
   loosejetEta_calo = 0;
   loosejetPhi_calo = 0;
   loosejetPassLooseID_calo = 0;
   loosejetPassTightID_calo = 0;
   loosejetEnergyFracHadronic_calo = 0;
   loosejetJECUncPlus_calo = 0;
   loosejetJECUncMinus_calo = 0;
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
   tightJetIndex_PF = 0;
   looseJetIndex_PF = 0;
   loosejetPt_PF = 0;
   loosejetPz_PF = 0;
   loosejetPtUncorr_PF = 0;
   loosejetEt_PF = 0;
   loosejetE_PF = 0;
   loosejetEta_PF = 0;
   loosejetPhi_PF = 0;
   loosejetPassLooseID_PF = 0;
   loosejetPassTightID_PF = 0;
   loosejetEnergyFracHadronic_PF = 0;
   loosejetJECUncPlus_PF = 0;
   loosejetJECUncMinus_PF = 0;
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
   muonIsRA2 = 0;
   muonIsGlobalMuon = 0;
   muonIsGlobalMuonPromptTight = 0;
   muonEoverP = 0;
   muonPtRatio = 0;
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
   eleIsRA2 = 0;
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
   muonIsRA2_PF = 0;
   muonIsGlobalMuon_PF = 0;
   muonIsGlobalMuonPromptTight_PF = 0;
   muonEoverP_PF = 0;
   muonPtRatio_PF = 0;
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
   eleIsRA2_PF = 0;
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
   pdfWeights = 0;
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
   fChain->SetBranchAddress("passTrigger", &passTrigger, &b_passTrigger);
   fChain->SetBranchAddress("hltPrescale", &hltPrescale, &b_hltPrescale);
   fChain->SetBranchAddress("pv_pass", &pv_pass, &b_pv_pass);
   fChain->SetBranchAddress("pv_isFake", &pv_isFake, &b_pv_isFake);
   fChain->SetBranchAddress("pv_z", &pv_z, &b_pv_z);
   fChain->SetBranchAddress("pv_rho", &pv_rho, &b_pv_rho);
   fChain->SetBranchAddress("pv_chi2", &pv_chi2, &b_pv_chi2);
   fChain->SetBranchAddress("pv_ndof", &pv_ndof, &b_pv_ndof);
   fChain->SetBranchAddress("muonPtDiff_PF", &muonPtDiff_PF, &b_muonPtDiff_PF);
   fChain->SetBranchAddress("passesBadPFMuonFilter", &passesBadPFMuonFilter, &b_passesBadPFMuonFilter);
   fChain->SetBranchAddress("passesInconsistentMuonPFCandidateFilter", &passesInconsistentMuonPFCandidateFilter, &b_passesInconsistentMuonPFCandidateFilter);
   fChain->SetBranchAddress("tightJetIndex_calo", &tightJetIndex_calo, &b_tightJetIndex_calo);
   fChain->SetBranchAddress("looseJetIndex_calo", &looseJetIndex_calo, &b_looseJetIndex_calo);
   fChain->SetBranchAddress("loosejetPt_calo", &loosejetPt_calo, &b_loosejetPt_calo);
   fChain->SetBranchAddress("loosejetPz_calo", &loosejetPz_calo, &b_loosejetPz_calo);
   fChain->SetBranchAddress("loosejetPtUncorr_calo", &loosejetPtUncorr_calo, &b_loosejetPtUncorr_calo);
   fChain->SetBranchAddress("loosejetEt_calo", &loosejetEt_calo, &b_loosejetEt_calo);
   fChain->SetBranchAddress("loosejetE_calo", &loosejetE_calo, &b_loosejetE_calo);
   fChain->SetBranchAddress("loosejetEta_calo", &loosejetEta_calo, &b_loosejetEta_calo);
   fChain->SetBranchAddress("loosejetPhi_calo", &loosejetPhi_calo, &b_loosejetPhi_calo);
   fChain->SetBranchAddress("loosejetPassLooseID_calo", &loosejetPassLooseID_calo, &b_loosejetPassLooseID_calo);
   fChain->SetBranchAddress("loosejetPassTightID_calo", &loosejetPassTightID_calo, &b_loosejetPassTightID_calo);
   fChain->SetBranchAddress("loosejetEnergyFracHadronic_calo", &loosejetEnergyFracHadronic_calo, &b_loosejetEnergyFracHadronic_calo);
   fChain->SetBranchAddress("loosejetJECUncPlus_calo", &loosejetJECUncPlus_calo, &b_loosejetJECUncPlus_calo);
   fChain->SetBranchAddress("loosejetJECUncMinus_calo", &loosejetJECUncMinus_calo, &b_loosejetJECUncMinus_calo);
   fChain->SetBranchAddress("loosejetFlavor_calo", &loosejetFlavor_calo, &b_loosejetFlavor_calo);
   fChain->SetBranchAddress("loosejetGenPt_calo", &loosejetGenPt_calo, &b_loosejetGenPt_calo);
   fChain->SetBranchAddress("loosejetGenEta_calo", &loosejetGenEta_calo, &b_loosejetGenEta_calo);
   fChain->SetBranchAddress("loosejetGenPhi_calo", &loosejetGenPhi_calo, &b_loosejetGenPhi_calo);
   fChain->SetBranchAddress("loosejetGenParticlePDGId_calo", &loosejetGenParticlePDGId_calo, &b_loosejetGenParticlePDGId_calo);
   fChain->SetBranchAddress("loosejetInvisibleEnergy_calo", &loosejetInvisibleEnergy_calo, &b_loosejetInvisibleEnergy_calo);
   fChain->SetBranchAddress("loosejetNTracks_calo", &loosejetNTracks_calo, &b_loosejetNTracks_calo);
   fChain->SetBranchAddress("loosejetNSV_calo", &loosejetNSV_calo, &b_loosejetNSV_calo);
   fChain->SetBranchAddress("loosejetSVUnWeightedMass_calo", &loosejetSVUnWeightedMass_calo, &b_loosejetSVUnWeightedMass_calo);
   fChain->SetBranchAddress("loosejetSVWeightedMass_calo", &loosejetSVWeightedMass_calo, &b_loosejetSVWeightedMass_calo);
   fChain->SetBranchAddress("loosejetSVUnWeightedLifetime_calo", &loosejetSVUnWeightedLifetime_calo, &b_loosejetSVUnWeightedLifetime_calo);
   fChain->SetBranchAddress("loosejetSVWeightedLifetime_calo", &loosejetSVWeightedLifetime_calo, &b_loosejetSVWeightedLifetime_calo);
   fChain->SetBranchAddress("loosejetSVUnWeightedCosTheta_calo", &loosejetSVUnWeightedCosTheta_calo, &b_loosejetSVUnWeightedCosTheta_calo);
   fChain->SetBranchAddress("loosejetSVWeightedCosTheta_calo", &loosejetSVWeightedCosTheta_calo, &b_loosejetSVWeightedCosTheta_calo);
   fChain->SetBranchAddress("loosejetBTagDisc_trackCountingHighPurBJetTags_calo", &loosejetBTagDisc_trackCountingHighPurBJetTags_calo, &b_loosejetBTagDisc_trackCountingHighPurBJetTags_calo);
   fChain->SetBranchAddress("loosejetBTagDisc_trackCountingHighEffBJetTags_calo", &loosejetBTagDisc_trackCountingHighEffBJetTags_calo, &b_loosejetBTagDisc_trackCountingHighEffBJetTags_calo);
   fChain->SetBranchAddress("loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags_calo", &loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags_calo, &b_loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags_calo);
   fChain->SetBranchAddress("loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags_calo", &loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags_calo, &b_loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags_calo);
   fChain->SetBranchAddress("loosejetBTagDisc_simpleSecondaryVertexBJetTags_calo", &loosejetBTagDisc_simpleSecondaryVertexBJetTags_calo, &b_loosejetBTagDisc_simpleSecondaryVertexBJetTags_calo);
   fChain->SetBranchAddress("tightJetIndex_PF", &tightJetIndex_PF, &b_tightJetIndex_PF);
   fChain->SetBranchAddress("looseJetIndex_PF", &looseJetIndex_PF, &b_looseJetIndex_PF);
   fChain->SetBranchAddress("loosejetPt_PF", &loosejetPt_PF, &b_loosejetPt_PF);
   fChain->SetBranchAddress("loosejetPz_PF", &loosejetPz_PF, &b_loosejetPz_PF);
   fChain->SetBranchAddress("loosejetPtUncorr_PF", &loosejetPtUncorr_PF, &b_loosejetPtUncorr_PF);
   fChain->SetBranchAddress("loosejetEt_PF", &loosejetEt_PF, &b_loosejetEt_PF);
   fChain->SetBranchAddress("loosejetE_PF", &loosejetE_PF, &b_loosejetE_PF);
   fChain->SetBranchAddress("loosejetEta_PF", &loosejetEta_PF, &b_loosejetEta_PF);
   fChain->SetBranchAddress("loosejetPhi_PF", &loosejetPhi_PF, &b_loosejetPhi_PF);
   fChain->SetBranchAddress("loosejetPassLooseID_PF", &loosejetPassLooseID_PF, &b_loosejetPassLooseID_PF);
   fChain->SetBranchAddress("loosejetPassTightID_PF", &loosejetPassTightID_PF, &b_loosejetPassTightID_PF);
   fChain->SetBranchAddress("loosejetEnergyFracHadronic_PF", &loosejetEnergyFracHadronic_PF, &b_loosejetEnergyFracHadronic_PF);
   fChain->SetBranchAddress("loosejetJECUncPlus_PF", &loosejetJECUncPlus_PF, &b_loosejetJECUncPlus_PF);
   fChain->SetBranchAddress("loosejetJECUncMinus_PF", &loosejetJECUncMinus_PF, &b_loosejetJECUncMinus_PF);
   fChain->SetBranchAddress("loosejetFlavor_PF", &loosejetFlavor_PF, &b_loosejetFlavor_PF);
   fChain->SetBranchAddress("loosejetGenPt_PF", &loosejetGenPt_PF, &b_loosejetGenPt_PF);
   fChain->SetBranchAddress("loosejetGenEta_PF", &loosejetGenEta_PF, &b_loosejetGenEta_PF);
   fChain->SetBranchAddress("loosejetGenPhi_PF", &loosejetGenPhi_PF, &b_loosejetGenPhi_PF);
   fChain->SetBranchAddress("loosejetGenParticlePDGId_PF", &loosejetGenParticlePDGId_PF, &b_loosejetGenParticlePDGId_PF);
   fChain->SetBranchAddress("loosejetInvisibleEnergy_PF", &loosejetInvisibleEnergy_PF, &b_loosejetInvisibleEnergy_PF);
   fChain->SetBranchAddress("loosejetNTracks_PF", &loosejetNTracks_PF, &b_loosejetNTracks_PF);
   fChain->SetBranchAddress("loosejetNSV_PF", &loosejetNSV_PF, &b_loosejetNSV_PF);
   fChain->SetBranchAddress("loosejetSVUnWeightedMass_PF", &loosejetSVUnWeightedMass_PF, &b_loosejetSVUnWeightedMass_PF);
   fChain->SetBranchAddress("loosejetSVWeightedMass_PF", &loosejetSVWeightedMass_PF, &b_loosejetSVWeightedMass_PF);
   fChain->SetBranchAddress("loosejetSVUnWeightedLifetime_PF", &loosejetSVUnWeightedLifetime_PF, &b_loosejetSVUnWeightedLifetime_PF);
   fChain->SetBranchAddress("loosejetSVWeightedLifetime_PF", &loosejetSVWeightedLifetime_PF, &b_loosejetSVWeightedLifetime_PF);
   fChain->SetBranchAddress("loosejetSVUnWeightedCosTheta_PF", &loosejetSVUnWeightedCosTheta_PF, &b_loosejetSVUnWeightedCosTheta_PF);
   fChain->SetBranchAddress("loosejetSVWeightedCosTheta_PF", &loosejetSVWeightedCosTheta_PF, &b_loosejetSVWeightedCosTheta_PF);
   fChain->SetBranchAddress("loosejetBTagDisc_trackCountingHighPurBJetTags_PF", &loosejetBTagDisc_trackCountingHighPurBJetTags_PF, &b_loosejetBTagDisc_trackCountingHighPurBJetTags_PF);
   fChain->SetBranchAddress("loosejetBTagDisc_trackCountingHighEffBJetTags_PF", &loosejetBTagDisc_trackCountingHighEffBJetTags_PF, &b_loosejetBTagDisc_trackCountingHighEffBJetTags_PF);
   fChain->SetBranchAddress("loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags_PF", &loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags_PF, &b_loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags_PF);
   fChain->SetBranchAddress("loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags_PF", &loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags_PF, &b_loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags_PF);
   fChain->SetBranchAddress("loosejetBTagDisc_simpleSecondaryVertexBJetTags_PF", &loosejetBTagDisc_simpleSecondaryVertexBJetTags_PF, &b_loosejetBTagDisc_simpleSecondaryVertexBJetTags_PF);
   fChain->SetBranchAddress("caloMET", &caloMET, &b_caloMET);
   fChain->SetBranchAddress("caloMETphi", &caloMETphi, &b_caloMETphi);
   fChain->SetBranchAddress("caloMETsig", &caloMETsig, &b_caloMETsig);
   fChain->SetBranchAddress("caloGenMET", &caloGenMET, &b_caloGenMET);
   fChain->SetBranchAddress("caloGenMETphi", &caloGenMETphi, &b_caloGenMETphi);
   fChain->SetBranchAddress("tcMET", &tcMET, &b_tcMET);
   fChain->SetBranchAddress("tcMETphi", &tcMETphi, &b_tcMETphi);
   fChain->SetBranchAddress("tcMETsig", &tcMETsig, &b_tcMETsig);
   fChain->SetBranchAddress("tcGenMET", &tcGenMET, &b_tcGenMET);
   fChain->SetBranchAddress("tcGenMETphi", &tcGenMETphi, &b_tcGenMETphi);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMETphi", &pfMETphi, &b_pfMETphi);
   fChain->SetBranchAddress("pfMETsig", &pfMETsig, &b_pfMETsig);
   fChain->SetBranchAddress("pfGenMET", &pfGenMET, &b_pfGenMET);
   fChain->SetBranchAddress("pfGenMETphi", &pfGenMETphi, &b_pfGenMETphi);
   fChain->SetBranchAddress("pf1MET", &pf1MET, &b_pf1MET);
   fChain->SetBranchAddress("pf1METphi", &pf1METphi, &b_pf1METphi);
   fChain->SetBranchAddress("pf1METsig", &pf1METsig, &b_pf1METsig);
   fChain->SetBranchAddress("pf1GenMET", &pf1GenMET, &b_pf1GenMET);
   fChain->SetBranchAddress("pf1GenMETphi", &pf1GenMETphi, &b_pf1GenMETphi);
   fChain->SetBranchAddress("trackPt", &trackPt, &b_trackPt);
   fChain->SetBranchAddress("trackEta", &trackEta, &b_trackEta);
   fChain->SetBranchAddress("trackPhi", &trackPhi, &b_trackPhi);
   fChain->SetBranchAddress("SumPtOverHT", &SumPtOverHT, &b_SumPtOverHT);
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
   fChain->SetBranchAddress("muonIsRA2", &muonIsRA2, &b_muonIsRA2);
   fChain->SetBranchAddress("muonIsGlobalMuon", &muonIsGlobalMuon, &b_muonIsGlobalMuon);
   fChain->SetBranchAddress("muonIsGlobalMuonPromptTight", &muonIsGlobalMuonPromptTight, &b_muonIsGlobalMuonPromptTight);
   fChain->SetBranchAddress("muonEoverP", &muonEoverP, &b_muonEoverP);
   fChain->SetBranchAddress("muonPtRatio", &muonPtRatio, &b_muonPtRatio);
   fChain->SetBranchAddress("muonPt", &muonPt, &b_muonPt);
   fChain->SetBranchAddress("muonEta", &muonEta, &b_muonEta);
   fChain->SetBranchAddress("muonPhi", &muonPhi, &b_muonPhi);
   fChain->SetBranchAddress("muonTrackIso", &muonTrackIso, &b_muonTrackIso);
   fChain->SetBranchAddress("muonEcalIso", &muonEcalIso, &b_muonEcalIso);
   fChain->SetBranchAddress("muonHcalIso", &muonHcalIso, &b_muonHcalIso);
   fChain->SetBranchAddress("muonChi2", &muonChi2, &b_muonChi2);
   fChain->SetBranchAddress("muonNdof", &muonNdof, &b_muonNdof);
   fChain->SetBranchAddress("muonNhits", &muonNhits, &b_muonNhits);
   fChain->SetBranchAddress("muonTrackd0", &muonTrackd0, &b_muonTrackd0);
   fChain->SetBranchAddress("muonTrackPhi", &muonTrackPhi, &b_muonTrackPhi);
   fChain->SetBranchAddress("muonPassID", &muonPassID, &b_muonPassID);
   fChain->SetBranchAddress("muonVtx_z", &muonVtx_z, &b_muonVtx_z);
   fChain->SetBranchAddress("eleIsRA2", &eleIsRA2, &b_eleIsRA2);
   fChain->SetBranchAddress("eleEt", &eleEt, &b_eleEt);
   fChain->SetBranchAddress("eleEta", &eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", &elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleTrackIso", &eleTrackIso, &b_eleTrackIso);
   fChain->SetBranchAddress("eleEcalIso", &eleEcalIso, &b_eleEcalIso);
   fChain->SetBranchAddress("eleHcalIso", &eleHcalIso, &b_eleHcalIso);
   fChain->SetBranchAddress("eledB", &eledB, &b_eledB);
   fChain->SetBranchAddress("eleVtx_z", &eleVtx_z, &b_eleVtx_z);
   fChain->SetBranchAddress("eleIDLoose", &eleIDLoose, &b_eleIDLoose);
   fChain->SetBranchAddress("eleIDRobustTight", &eleIDRobustTight, &b_eleIDRobustTight);
   fChain->SetBranchAddress("elePassID", &elePassID, &b_elePassID);
   fChain->SetBranchAddress("muonIsRA2_PF", &muonIsRA2_PF, &b_muonIsRA2_PF);
   fChain->SetBranchAddress("muonIsGlobalMuon_PF", &muonIsGlobalMuon_PF, &b_muonIsGlobalMuon_PF);
   fChain->SetBranchAddress("muonIsGlobalMuonPromptTight_PF", &muonIsGlobalMuonPromptTight_PF, &b_muonIsGlobalMuonPromptTight_PF);
   fChain->SetBranchAddress("muonEoverP_PF", &muonEoverP_PF, &b_muonEoverP_PF);
   fChain->SetBranchAddress("muonPtRatio_PF", &muonPtRatio_PF, &b_muonPtRatio_PF);
   fChain->SetBranchAddress("muonPt_PF", &muonPt_PF, &b_muonPt_PF);
   fChain->SetBranchAddress("muonEta_PF", &muonEta_PF, &b_muonEta_PF);
   fChain->SetBranchAddress("muonPhi_PF", &muonPhi_PF, &b_muonPhi_PF);
   fChain->SetBranchAddress("muonTrackIso_PF", &muonTrackIso_PF, &b_muonTrackIso_PF);
   fChain->SetBranchAddress("muonEcalIso_PF", &muonEcalIso_PF, &b_muonEcalIso_PF);
   fChain->SetBranchAddress("muonHcalIso_PF", &muonHcalIso_PF, &b_muonHcalIso_PF);
   fChain->SetBranchAddress("muonChi2_PF", &muonChi2_PF, &b_muonChi2_PF);
   fChain->SetBranchAddress("muonNdof_PF", &muonNdof_PF, &b_muonNdof_PF);
   fChain->SetBranchAddress("muonNhits_PF", &muonNhits_PF, &b_muonNhits_PF);
   fChain->SetBranchAddress("muonTrackd0_PF", &muonTrackd0_PF, &b_muonTrackd0_PF);
   fChain->SetBranchAddress("muonTrackPhi_PF", &muonTrackPhi_PF, &b_muonTrackPhi_PF);
   fChain->SetBranchAddress("muonPassID_PF", &muonPassID_PF, &b_muonPassID_PF);
   fChain->SetBranchAddress("muonVtx_z_PF", &muonVtx_z_PF, &b_muonVtx_z_PF);
   fChain->SetBranchAddress("eleIsRA2_PF", &eleIsRA2_PF, &b_eleIsRA2_PF);
   fChain->SetBranchAddress("eleEt_PF", &eleEt_PF, &b_eleEt_PF);
   fChain->SetBranchAddress("eleEta_PF", &eleEta_PF, &b_eleEta_PF);
   fChain->SetBranchAddress("elePhi_PF", &elePhi_PF, &b_elePhi_PF);
   fChain->SetBranchAddress("eleTrackIso_PF", &eleTrackIso_PF, &b_eleTrackIso_PF);
   fChain->SetBranchAddress("eleEcalIso_PF", &eleEcalIso_PF, &b_eleEcalIso_PF);
   fChain->SetBranchAddress("eleHcalIso_PF", &eleHcalIso_PF, &b_eleHcalIso_PF);
   fChain->SetBranchAddress("eledB_PF", &eledB_PF, &b_eledB_PF);
   fChain->SetBranchAddress("eleVtx_z_PF", &eleVtx_z_PF, &b_eleVtx_z_PF);
   fChain->SetBranchAddress("eleIDLoose_PF", &eleIDLoose_PF, &b_eleIDLoose_PF);
   fChain->SetBranchAddress("eleIDRobustTight_PF", &eleIDRobustTight_PF, &b_eleIDRobustTight_PF);
   fChain->SetBranchAddress("elePassID_PF", &elePassID_PF, &b_elePassID_PF);
   fChain->SetBranchAddress("SUSY_nb", &SUSY_nb, &b_SUSY_nb);
   fChain->SetBranchAddress("qScale", &qScale, &b_qScale);
   fChain->SetBranchAddress("ptHat", &ptHat, &b_ptHat);
   fChain->SetBranchAddress("mcWeight", &mcWeight, &b_mcWeight);
   fChain->SetBranchAddress("pdfWeights", &pdfWeights, &b_pdfWeights);
   fChain->SetBranchAddress("topDecayCode", &topDecayCode, &b_topDecayCode);
   fChain->SetBranchAddress("ZDecayMode", &ZDecayMode, &b_ZDecayMode);
   fChain->SetBranchAddress("flavorHistory", &flavorHistory, &b_flavorHistory);
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

bool basicLoop::setCutScheme(CutScheme cutscheme) {

  if (cutscheme == nCutSchemes) return false;
  theCutScheme_ = cutscheme;

  cutNames_.clear();
  cutTags_.clear();
  //  cutMap_.clear();

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
  //  for (unsigned int i=0;i<cutTags_.size(); i++) {
  //    cutMap_[cutTags_[i]] = i;
  //  }
  
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

    cutTags_.push_back("cut2SUSYb"); cutNames_[ cutTags_.back()] = "==2SUSYb";
    cutTags_.push_back("cut4SUSYb"); cutNames_[ cutTags_.back()] = "==4SUSYb";

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
    cutTags_.push_back("cutCleaning"); cutNames_[cutTags_.back()]="TailCleaning";

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

  //check if we are *requiring* this cut (special from the normal scheme)
  for (unsigned int i = 0; i< requiredCut_.size() ; i++) {
    if ( cutTag == requiredCut_.at(i) ) return true;
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

    else if (cutTag == "cut2SUSYb")  cutIsRequired = false;
    else if (cutTag == "cut4SUSYb")  cutIsRequired = false;

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
    else if (cutTag == "cutCleaning") cutIsRequired = theCleaningType_ != kNoCleaning;
    else assert(0);
  }
  else assert(0);

  return cutIsRequired;
}


bool basicLoop::passMuVeto() {
  //using nmu_ allows for inverting the veto
  return (countMu() == nmu_);
}

int basicLoop::countMu() {

  int ngoodmu=0;

  unsigned int nmu = 0;
  if (theLeptonType_ == kNormal) nmu=muonPt->size();  
  else if (theLeptonType_ ==kPFLeptons) nmu=muonPt_PF->size();
  else if (theLeptonType_ ==kPFLeptonsRA2) nmu=muonPt_PF->size();
  else {assert(0);}

  for ( unsigned int i = 0; i< nmu; i++) {
    float pt=0,eta=0,reliso=0;
    bool isglobal=false;    

    //load information depending on lepton type
    if (theLeptonType_ == kNormal) {
      pt = muonPt->at(i);
      eta = muonEta->at(i);
      reliso = (muonTrackIso->at(i) + muonHcalIso->at(i) + muonEcalIso->at(i))/muonPt->at(i);
      isglobal = muonIsGlobalMuon->at(i);
    }
    else if (theLeptonType_ ==kPFLeptons ||theLeptonType_ ==kPFLeptonsRA2) {
      pt = muonPt_PF->at(i);
      eta = muonEta_PF->at(i);
      reliso = (muonTrackIso_PF->at(i) + muonHcalIso_PF->at(i) + muonEcalIso_PF->at(i))/muonPt_PF->at(i);
      isglobal = muonIsGlobalMuon_PF->at(i);
    }
    else {assert(0);}

    //now make cuts
    if (theLeptonType_ == kPFLeptonsRA2) {
      if ( ! muonIsRA2_PF->at(i) ) continue;
    }
    else { //not bothering with the assert anymore!
      if ( !isglobal ) continue;
      if ( pt < 10 ) continue; 
      if (fabs(eta) > 2.5) continue;
      if ( reliso > 0.2) continue;
    }

    //once we reach here we've got a good muon in hand
    //FIXME adding this as a hack...i would like to do this more elegantly, but for now this will work
    if (ngoodmu==0) muonpt1_ = pt;

    ++ngoodmu;
  }
  
  return ngoodmu;
}

bool basicLoop::passTauVeto() {
  //TO DO
  return true;
}

bool basicLoop::passEleVeto() {
  //see comments in muon section

  //using ne_ allows for inverting the veto
  return (countEle() == ne_);
}

int basicLoop::countEle() {

  int ngoodele=0;

  unsigned int nele = 0;
  if (theLeptonType_ == kNormal) nele=eleEt->size();  
  else if (theLeptonType_ ==kPFLeptons) nele=eleEt_PF->size();
  else if (theLeptonType_ ==kPFLeptonsRA2) nele=eleEt_PF->size();
  else {assert(0);}

  for (unsigned int i=0; i < nele; i++) {
    float et=0,eta=0,reliso=0;
    
    //load information depending on lepton type
    if (theLeptonType_ == kNormal) {
      et = eleEt->at(i);
      eta = eleEta->at(i);
      reliso = (eleTrackIso->at(i) + eleHcalIso->at(i) + eleEcalIso->at(i))/eleEt->at(i);
    }
    else if (theLeptonType_ ==kPFLeptons || theLeptonType_ ==kPFLeptonsRA2) {
      et = eleEt_PF->at(i);
      eta = eleEta_PF->at(i);
      reliso = (eleTrackIso_PF->at(i) + eleHcalIso_PF->at(i) + eleEcalIso_PF->at(i))/eleEt_PF->at(i);
    }
    else {assert(0);}

    //make cuts    
    if (theLeptonType_ == kPFLeptonsRA2) {
      if ( ! eleIsRA2_PF->at(i) ) continue;
    }
    else {
      if ( et < 15 ) continue;
      if ( fabs(eta) > 2.5 ) continue;
      if ( reliso > 0.2) continue;
    }
    //if any electron passes all of these cuts, then it is good

    //FIXME adding this as a hack...i would like to do this more elegantly, but for now this will work
    if (ngoodele==0) eleet1_ = et;

    ++ngoodele;
  }

  return ngoodele;

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

//the standard code -- applies the simple OR of what is in triggerList_
bool basicLoop::passHLT() { 

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


/* modification of the above code to require that the triggers in the OR have prescale==1
bool basicLoop::passHLT() {

  if (!printedHLT_) {cout<<"WARNING -- using modified HLT filter!"<<endl; printedHLT_=true;}

  for (std::map<TString, int>::const_iterator itrig=triggerList_.begin(); 
       itrig!=triggerList_.end(); ++itrig) {

    if ( passTrigger->at( itrig->second) && (hltPrescale->at(itrig->second) == 1) ) {
      return true;
    }
  }
  return false;
}
*/

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
  if (cutTag=="cutMuVeto" && theCutScheme_==kRA2) return passMuVeto();
  else if (cutTag=="cutMuVeto" && theCutScheme_==kSync1) return passMuVeto();
  else if (cutTag=="cutMuVeto" && theCutScheme_==kBaseline0) return passMuVeto(); //passMuVetoDon()

  if (cutTag=="cutEleVeto" && theCutScheme_==kRA2) return passEleVeto();
  else if (cutTag=="cutEleVeto" && theCutScheme_==kSync1) return passEleVeto();
  else if (cutTag=="cutEleVeto" && theCutScheme_==kBaseline0) return passEleVeto(); //passEleVetoDon()

  if (cutTag=="cutTauVeto" && theCutScheme_==kBaseline0) return passTauVeto();

  if (cutTag=="cut3Jets" && theCutScheme_==kSync1) return (nGoodJets_Sync1() >= 3);

  if (theCutScheme_==kSync1 ) {
    if (cutTag=="cutJetPt1") return jetPtOfN(1)>100;
    if (cutTag=="cutJetPt2") return jetPtOfN(2)>50;
    if (cutTag=="cutJetPt3") return jetPtOfN(3)>50;;
  }

  if (cutTag=="cut3Jets" && theJetType_==kPF && theCutScheme_==kRA2) return (tightJetIndex->size() >= 3);

  if (cutTag=="cutHT" && theJetType_==kPF && theCutScheme_==kRA2) {
    cout<<"WARNING -- RA2 cut scheme needs reimplementation!"<<endl;
    //calculate HT
    double pfHT = 0;
    //    for (unsigned int i=0; i<tightJetIndex->size(); i++){
    //      pfHT+=jetPt.at(i);
    //    }
    return (pfHT>300.);
  }
  else if (cutTag=="cutHT" && theCutScheme_==kSync1) return getHT_Sync1() >300;

  //MET
  //it is now an anachronism that we treat MET and MHT as different cut flow steps
  //this block of code should now evaluate *the same* for both steps
  if (cutTag == "cutMET" || cutTag=="cutMHT") {
    float mymet = getMET();

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
    cout<<"WARNING -- kDeltaPhi is implemented in an old way (only calo jets)"<<endl;
    float phi_of_MET = getMETphi();
    
    // FIXME loose jet def'n has changed! this needs a careful update
    int nloosejets = loosejetPhi->size();
    //need to calculate DeltaPhi between jets and MET
    //for RA2 and MHT, this is done with *loose* jets!
    double dp0=0,dp1=0,dp2=0;
    if (nloosejets>0) {dp0=getDeltaPhi( loosejetPhi->at(0), phi_of_MET);} else {return false;}
    if (nloosejets>1) {dp1=getDeltaPhi( loosejetPhi->at(1), phi_of_MET);} else {return false;}
    if (nloosejets>2) {dp2=getDeltaPhi( loosejetPhi->at(2), phi_of_MET);} else {return false;}
    //here is the implementation of the DeltaPhi cuts
    if ( dp0 >0.3 && dp1 >0.5 && dp2 >0.3 ) { return true; } else {return false;}
  }
  else if (cutTag == "cutDeltaPhi" && //replace normal DeltaPhi cut with minDeltaPhi cut
	   ( theDPType_ == kminDP ) ) {
    
    if (theMETType_ ==kMET ||theMETType_==ktcMET ||theMETType_==kpfMET ||theMETType_ == kMHT) return ( getMinDeltaPhiMET(3) >= 0.3 );
    else {assert(0);}
  }
  else if (cutTag == "cutDeltaPhi" && //replace normal DeltaPhi cut with minDeltaPhi cut
	   ( theDPType_ == kminDPAll30 ) ) {
    
    if (theMETType_ ==kMET ||theMETType_==ktcMET ||theMETType_==kpfMET ||theMETType_ == kMHT) return ( getMinDeltaPhiMET30(99) >= 0.3 );
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
    if (cutTag == "cutEq1b") return nbSSVM == 1;
  }
  else if (theCutScheme_==kSync1) {
    if (cutTag == "cut1b" || cutTag == "cut2b" || cutTag == "cut3b" || cutTag=="cutEq1b") {
      int nb = countBJets_Sync1();
      if (cutTag == "cut1b") return nb >=1;
      if (cutTag == "cut2b") return nb >=2;
      if (cutTag == "cut3b") return nb >=3;
      if (cutTag == "cutEq1b") return nb == 1;
    }
  }

  if (cutTag == "cutCleaning") return passCleaning();

  if (cutTag == "cut2SUSYb") return (SUSY_nb == 2);
  if (cutTag == "cut4SUSYb") return (SUSY_nb == 4);

  //no longer storing cut results in ntuple!
  if (cutTag!="cutInclusive") {
    cout<<"[passCut] should not have reached this point! "<<cutTag<<endl;
    assert(0);
  }
  return true;
}

int basicLoop::countGoodPV() {

  //never run the code down here (except for special checks)

  int npass=0;
  for (unsigned int ipv = 0; ipv<pv_isFake->size(); ipv++) {
    if ( pv_isFake->at(ipv) ) continue;
    if ( fabs(pv_z->at(ipv)) > 24 ) continue;
    if ( fabs(pv_rho->at(ipv)) > 2 ) continue;
    if ( pv_ndof->at(ipv) <= 4 ) continue;

    ++npass;
  }

  return npass; 
}

bool basicLoop::passPV() {


  //i have now understood the discrepancies between the PVSelector and the hand-calculated result
  //the PV selector only looks at the first PV.

  //this is precomputed PVSelector result
  return pv_pass;

}

int basicLoop::countGenBJets(float threshold) {
  int nb=0;
  for (unsigned int i=0; i<loosejetGenPt->size(); ++i) {
    if ( (loosejetGenPt->at(i) > threshold)
	 && ( abs(loosejetFlavor->at(i)) == 5 ) ) ++nb;
  }
  return nb;
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

    if (getLooseJetPt(i) <30) continue;
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

  std::cout<<"getUncorrectedHT needs to be reimplemented! "<<threshold <<std::endl;
  return -1;

  /*
I have realized that if I want to do a crude trigger simulation, then I need to store every jet
with Uncorrect Pt>20 GeV, not just jets passing my 'loose' cuts

   vector<float>   *loosejetPtUncorr;
   vector<float>   *loosejetEtaUncorr;
   vector<float>   *loosejetPhiUncorr;*/

  /*
  double ht=0;
  for (unsigned int ijet=0; ijet<veryloosejetPtUncorr->size(); ijet++) {

    if ( fabs(veryloosejetEtaUncorr->at(ijet)) > 5 ) continue; //what is the eta range actually used by the HLT?

    double thispt = veryloosejetPtUncorr->at(ijet);
    if (thispt>=threshold)   ht+= thispt;
  }
  
   return ht;
  */
}

std::pair<float,float> basicLoop::getUnclusteredSmearedMETxy() {
  assert( theLeptonType_ == kPFLeptons || theLeptonType_ == kPFLeptonsRA2 );  //hard-code for PF leptons for now

  /*
twiki says to make sure we use JER bias-corrected uncorrected jet pT.
so we set that automatically in setMETuncType()
  */

  float myMET=-1;
  float myMETphi=-99;

  //adjust for global MET type
  if (theMETType_== kMET) {
    myMET= caloMET;
    myMETphi= caloMETphi;
  }
  else if (theMETType_ == ktcMET) {
    myMET= tcMET;
    myMETphi= tcMETphi;
  }
  else if (theMETType_ == kpfMET) {
    myMET= pfMET;
    myMETphi= pfMETphi;
  }
  //  else if (theMETType_ == kMHT) {  } //not implemented
  else {  assert(0);  }

  float myMETx = myMET * cos(myMETphi);
  float myMETy = myMET * sin(myMETphi);
  
  for (unsigned int ijet = 0; ijet< loosejetPt->size(); ++ijet) {
    float jetUx = getLooseJetPtUncorr(ijet) * cos(loosejetPhi->at(ijet));
    float jetUy = getLooseJetPtUncorr(ijet) * sin(loosejetPhi->at(ijet));

    myMETx += jetUx;
    myMETy += jetUy;
  }

  for (unsigned int imu = 0; imu< muonPt_PF->size(); ++imu) {
    float mux = muonPt_PF->at(imu) * cos(muonPhi_PF->at(imu));
    float muy = muonPt_PF->at(imu) * sin(muonPhi_PF->at(imu));

    myMETx += mux;
    myMETy += muy;
  }
  for (unsigned int iel = 0; iel< eleEt_PF->size(); ++iel) {
    float elx = eleEt_PF->at(iel) * cos(elePhi_PF->at(iel));
    float ely = eleEt_PF->at(iel) * sin(elePhi_PF->at(iel));

    myMETx += elx;
    myMETy += ely;
  }

  //now scale unclustered MET
  float factor=1;
  if (theMETuncType_ ==kMETuncDown) factor=0.9;
  else if (theMETuncType_ ==kMETuncUp) factor=1.1;
  else {assert(0);}
  myMETx *= factor;
  myMETy *= factor;

  //now repeat all of the loops but do -= instead of +=
  for (unsigned int ijet = 0; ijet< loosejetPt->size(); ++ijet) {
    float jetUx = getLooseJetPtUncorr(ijet) * cos(loosejetPhi->at(ijet));
    float jetUy = getLooseJetPtUncorr(ijet) * sin(loosejetPhi->at(ijet));

    myMETx -= jetUx;
    myMETy -= jetUy;
  }

  for (unsigned int imu = 0; imu< muonPt_PF->size(); ++imu) {
    float mux = muonPt_PF->at(imu) * cos(muonPhi_PF->at(imu));
    float muy = muonPt_PF->at(imu) * sin(muonPhi_PF->at(imu));

    myMETx -= mux;
    myMETy -= muy;
  }
  for (unsigned int iel = 0; iel< eleEt_PF->size(); ++iel) {
    float elx = eleEt_PF->at(iel) * cos(elePhi_PF->at(iel));
    float ely = eleEt_PF->at(iel) * sin(elePhi_PF->at(iel));

    myMETx -= elx;
    myMETy -= ely;
  }

  return make_pair(myMETx,myMETy);
}

float basicLoop::getJESExtraUnc(unsigned int ijet) {
  //add all of these pieces in quadrature

  float  unc = 0.015 * 0.015; //c_sw
  //          e_PU   JA  AvgPU    pT
  unc += pow(0.75 * 0.8 * 2.2 / loosejetPt->at(ijet) ,2);
  //b jet uncertainty
  if (abs(loosejetFlavor->at(ijet)) == 5) {
    if ( loosejetPt->at(ijet)>50 && loosejetPt->at(ijet)<200 && fabs(loosejetEta->at(ijet))<2 ) unc += 0.02*0.02;
    else unc += 0.03*0.03;
  }

  return sqrt(unc);
}

std::pair<float,float> basicLoop::getJESAdjustedMETxy() {

  if (theJESType_ == kJES0) {assert(0);}

  float myMET=-1;
  float myMETphi=-99;

  //adjust for global MET type
  if (theMETType_== kMET) {
    myMET= caloMET;
    myMETphi= caloMETphi;
  }
  else if (theMETType_ == ktcMET) {
    myMET= tcMET;
    myMETphi= tcMETphi;
  }
  else if (theMETType_ == kpfMET) {
    myMET= pfMET;
    myMETphi= pfMETphi;
  }
  //  else if (theMETType_ == kMHT) {  } //not implemented
  else {  assert(0);  }

  float myMETx = myMET * cos(myMETphi);
  float myMETy = myMET * sin(myMETphi);
  
  //loop over all jets
  for (unsigned int ijet=0; ijet<loosejetPt->size(); ++ijet) {
    float jetUx = loosejetPtUncorr->at(ijet) * cos(loosejetPhi->at(ijet));
    float jetUy = loosejetPtUncorr->at(ijet) * sin(loosejetPhi->at(ijet));

    myMETx += jetUx;
    myMETy += jetUy;
    float jes_factor=1;
    if (theJESType_ == kJESup) {
      float unc =   loosejetJECUncPlus->at(ijet);
      float cor = getJESExtraUnc(ijet);
      unc = sqrt(unc*unc + cor*cor);
      jes_factor += unc;
    }
    else if (theJESType_ == kJESdown) {
      float unc =  loosejetJECUncMinus->at(ijet);
      float cor = getJESExtraUnc(ijet);
      jes_factor -= sqrt(unc*unc + cor*cor);
    }
    else {assert(0);}
    jetUx *= jes_factor;
    jetUy *= jes_factor;

    myMETx -= jetUx;
    myMETy -= jetUy;
  }

  return make_pair(myMETx, myMETy);
}

std::pair<float,float> basicLoop::getJERAdjustedMETxy() {
  if (theJERType_ == kJER0) {assert(0);}

  float myMET=-1;
  float myMETphi=-99;

  //adjust for global MET type
  if (theMETType_== kMET) {
    myMET= caloMET;
    myMETphi= caloMETphi;
  }
  else if (theMETType_ == ktcMET) {
    myMET= tcMET;
    myMETphi= tcMETphi;
  }
  else if (theMETType_ == kpfMET) {
    myMET= pfMET;
    myMETphi= pfMETphi;
  }
  else if (theMETType_ == kMHT)  return getJERAdjustedMHTxy();
  else {  assert(0);  }

  float myMETx = myMET * cos(myMETphi);
  float myMETy = myMET * sin(myMETphi);
  
  //loop over all jets
  for (unsigned int ijet=0; ijet<loosejetPt->size(); ++ijet) {
    float genpt = loosejetGenPt->at(ijet);

    if (genpt < 15) continue;

    float jetUx = loosejetPtUncorr->at(ijet) * cos(loosejetPhi->at(ijet));
    float jetUy = loosejetPtUncorr->at(ijet) * sin(loosejetPhi->at(ijet));

    myMETx += jetUx;
    myMETy += jetUy;

    float recopt = loosejetPt->at(ijet);
    float factor = 0;
    if (theJERType_ == kJERbias )    factor = 0.1; //hard-coded factors from top twiki
    else if (theJERType_ == kJERup)  factor = 0.2;
    else {assert(0);}

    float  deltapt = (recopt - genpt) * factor;
    float frac = (recopt+deltapt)/recopt;
    float ptscale = frac>0 ? frac : 0;
        //

    jetUx *= ptscale;
    jetUy *= ptscale;

    myMETx -= jetUx;
    myMETy -= jetUy;
  }
  return make_pair(myMETx, myMETy);

}

float basicLoop::getMET() {
  float myMET=-1;

  if (theJESType_ == kJES0 && theJERType_ == kJER0 && theMETuncType_== kMETunc0) {
    //adjust for global MET type
    if (theMETType_== kMET) {
      myMET= caloMET;
    }
    else if (theMETType_ == ktcMET) {
      myMET= tcMET;
    }
    else if (theMETType_ == kpfMET) {
      myMET= pfMET;
    }
    else if (theMETType_ == kMHT) {
      myMET= getMHT();
    }
    else { 
      assert(0);
    }
  }
  else if (theJESType_ != kJES0 && theJERType_ == kJER0 && theMETuncType_== kMETunc0) {
    std::pair<float, float> metxy = getJESAdjustedMETxy();
    myMET = sqrt( metxy.first*metxy.first + metxy.second*metxy.second);
  }
  else if (theJESType_ == kJES0 && theJERType_ != kJER0 && theMETuncType_== kMETunc0) {
    std::pair<float, float> metxy = getJERAdjustedMETxy();
    myMET = sqrt( metxy.first*metxy.first + metxy.second*metxy.second);
  }
  else if (theJESType_ == kJES0 && theJERType_ == kJERbias && theMETuncType_!= kMETunc0) {
    std::pair<float, float> metxy = getUnclusteredSmearedMETxy();
    myMET = sqrt( metxy.first*metxy.first + metxy.second*metxy.second);
  }
  else {assert(0);}

  return myMET;
}

float basicLoop::getMETphi() {

  float myMETphi=-99;

  if (theJESType_ == kJES0 && theJERType_ == kJER0 && theMETuncType_== kMETunc0) {
    //adjust for global MET type
    if (theMETType_== kMET) {
      myMETphi = caloMETphi;
    }
    else if (theMETType_ == ktcMET) {
      myMETphi = tcMETphi;
    }
    else if (theMETType_ == kpfMET) {
      myMETphi = pfMETphi;
    }
    else if (theMETType_ == kMHT) {
      myMETphi = getMHTphi();
    }
    else {  assert(0); }
  }
  else if (theJESType_ != kJES0 && theJERType_ == kJER0 && theMETuncType_== kMETunc0) { //JES uncertainty
    std::pair<float, float> metxy = getJESAdjustedMETxy();
    myMETphi = atan2(metxy.second, metxy.first);
  }
  else if (theJESType_ == kJES0 && theJERType_ != kJER0 && theMETuncType_== kMETunc0) { //JER uncertainty
    std::pair<float, float> metxy = getJERAdjustedMETxy();
    myMETphi = atan2(metxy.second, metxy.first);
  }
  else if (theJESType_ == kJES0 && theJERType_ == kJERbias && theMETuncType_!= kMETunc0) {
    std::pair<float, float> metxy = getUnclusteredSmearedMETxy();
    myMETphi = atan2(metxy.second, metxy.first);
  }
  else {assert(0);}

  return myMETphi;
}

float basicLoop::getGenMET() {
  //for just one number i am happy to do this by value

  //adjust for global MET type
  if (theMETType_== kMET) {
    return caloGenMET;
  }
  else if (theMETType_ == ktcMET) {
    return tcGenMET;
  }
  else if (theMETType_ == kpfMET) {
    return pfGenMET;
  }
  else if (theMETType_ == kMHT) {
    assert(0); //not implemented
  }
  
  assert(0);

  return 0;
}

float basicLoop::getGenMETphi() {
  //for just one number i am happy to do this by value

  //adjust for global MET type
  if (theMETType_== kMET) {
    return caloGenMETphi;
  }
  else if (theMETType_ == ktcMET) {
    return tcGenMETphi;
  }
  else if (theMETType_ == kpfMET) {
    return pfGenMETphi;
  }
  else if (theMETType_ == kMHT) {
    assert(0); //not implemented
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
    cout<<"WARNING -- getMinDeltaPhiMET needs reimplementation for legacy cut schemes!"<<endl;
    /*
    unsigned int njets=  jetPhi.size();
    if (njets < maxjets) maxjets = njets;
    
    //get the minimum angle between the first n jets and MET
    for (unsigned int i=0; i< maxjets; i++) {
      
      double dp =  getDeltaPhi( jetPhi.at(i) , getMETphi());
      if (dp<mindp) mindp=dp;
      
    }
    */
  }

  return mindp;
}

double basicLoop::getMaxDeltaPhiMET30(unsigned int maxjets) {
  double maxdp=-99;

  if (theCutScheme_==kBaseline0) {
    unsigned int ngood=0;
    //get the maximum angle between the first n jets and MET
    for (unsigned int i=0; i< loosejetPhi->size(); i++) {

      if (isGoodJet30(i)) {
	++ngood;
	double dp =  getDeltaPhi( loosejetPhi->at(i) , getMETphi());
	if (dp>maxdp) maxdp=dp;
	if (ngood >= maxjets) break;
      }
    }
  }
  else {
    assert(0);
  }

  return maxdp;
}

double basicLoop::getMinDeltaPhiMET30(unsigned int maxjets) {

  double mindp=99;

  if (theCutScheme_==kBaseline0) {
    unsigned int ngood=0;
    //get the minimum angle between the first n jets and MET
    for (unsigned int i=0; i< loosejetPhi->size(); i++) {

      if (isGoodJet30(i)) {
	++ngood;
	double dp =  getDeltaPhi( loosejetPhi->at(i) , getMETphi());
	if (dp<mindp) mindp=dp;
	if (ngood >= maxjets) break;
      }
    }
  }
  else {
    assert(0);
  }

  return mindp;
}

double basicLoop::getMinDeltaPhiMET30_eta5(unsigned int maxjets) {

  double mindp=99;

  if (theCutScheme_==kBaseline0) {
    unsigned int ngood=0;
    //get the minimum angle between the first n jets and MET
    for (unsigned int i=0; i< loosejetPhi->size(); i++) {
      
      bool passJetId=false;
      if ( loosejetPassLooseID->at(i) ) passJetId=true;
      if (passJetId && getLooseJetPt(i)>30 && fabs(loosejetEta->at(i))<5.0){
	++ngood;
	double dp =  getDeltaPhi( loosejetPhi->at(i) , getMETphi());
	if (dp<mindp) mindp=dp;
	if (ngood >= maxjets) break;
      }
    }
  }
  else {
    assert(0);
  }
  
  return mindp;
}

double basicLoop::getMinDeltaPhiMET30_eta5_noId(unsigned int maxjets) {

  double mindp=99;

  if (theCutScheme_==kBaseline0) {
    unsigned int ngood=0;
    //get the minimum angle between the first n jets and MET
    for (unsigned int i=0; i< loosejetPhi->size(); i++) {
      
      bool passJetId=true;
      if (passJetId && getLooseJetPt(i)>30 && fabs(loosejetEta->at(i))<5.0){
	++ngood;
	double dp =  getDeltaPhi( loosejetPhi->at(i) , getMETphi());
	if (dp<mindp) mindp=dp;
	if (ngood >= maxjets) break;
      }
    }
  }
  else {
    assert(0);
  }
  
  return mindp;
}


double basicLoop::getMinDeltaPhibMET() {
  //get the minimum angle between a b jet and MET
  /* uses default tight jets and default MET */

  cout<<"WARNING -- getMinDeltaPhibMET() needs reimplemenation"<<endl;
  
  double mindp=99;
/*
  for (unsigned int i=0; i<jetPhi.size(); i++) {
    if (passBCut(i) ) {
      double dp =  getDeltaPhi( jetPhi.at(i), getMETphi());
      if (dp<mindp) mindp=dp;
    }
  }
*/
  return mindp;
}

double basicLoop::getDeltaPhib1b2() {
  //get the angle between the lead 2 b jets
  /* legacy code uses tight jets ; kBaseline0 uses loose jets and b cuts*/

  assert(theCutScheme_==kBaseline0);

  std::vector<float> phis;

  //careful...we're using tight jets for one case and loose jets for the other
  unsigned int maxj = loosejetPhi->size() ;

  for (unsigned int i=0; i< maxj; i++) {
    bool isGoodB = isGoodJet30(i) && passSSVM(i);
    
    if (isGoodB ) {
      
      float thisJetPhi = loosejetPhi->at(i) ;
      phis.push_back( thisJetPhi);

      if (phis.size() == 2) break;
    }
  }

  //make the code safe for events with less than 2 b tags
  if (phis.size() < 2) return -99;

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

double basicLoop::getDeltaPhi(double phi1, double phi2) {

  return acos(cos(phi1-phi2));
}

float basicLoop::getDeltaPhiTopTwoJets() {

  float phi1=0,phi2=0;
  int ngood=0;
  for (unsigned int i=0; i<loosejetPt->size(); i++) {
    if (isGoodJet(i)) {
      ngood++;
      if (ngood==1) phi1 = loosejetPhi->at(i);
      if (ngood==2) phi2 = loosejetPhi->at(i);
      else break;
    }
  }

  return getDeltaPhi(phi1,phi2);
}

double basicLoop::getMinDeltaPhi_bj(unsigned int bindex) {
  /*
this function assumes that bindex is of a b jet. it doesn't verify it
  */

  //using tight jets
  cout<<"WARNING -- getMinDeltaPhi_bj needs reimplementation for Baseline0"<<endl;
  double minDeltaPhi=99;
  /*

  double bphi = jetPhi.at(bindex);

  //loop over the jets
  for (unsigned int jindex=0; jindex<jetPhi.size(); jindex++) {
    if ( !passBCut(jindex) ) { //only look at non-b jets
      double dp = getDeltaPhi(bphi, jetPhi.at(jindex));
      if (dp < minDeltaPhi) minDeltaPhi = dp;
    }
  }
  */
  return minDeltaPhi;
}

//double basicLoop::getOverallMinDeltaR_bj() {
  /*
this code was written for a study that proved to be not-so-useful
still uses the 'tight' jets stored by the ntuple maker
  */

/*
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
*/

std::pair<float,float> basicLoop::getJERAdjustedMHTxy() {
  //no difference from MHT calculation -- just return both pieces
  float mhtx=0;
  float mhty=0;

  for (unsigned int i=0; i<loosejetPt->size(); i++) {
    if (isGoodJetMHT( i ) ) {

      mhtx -= getLooseJetPt(i) * cos(loosejetPhi->at(i));
      mhty -= getLooseJetPt(i) * sin(loosejetPhi->at(i));
    }
  }

  //this is an experiment! add taus in!
  //-->experiment successful! leave taus in.
  //FIXME hard-coded for PF
  for (unsigned int i=0; i<tauPt_PF->size(); i++) {
    //start by using same pT threshold as jets
    if ( tauPt_PF->at(i) > 30 && fabs(tauEta_PF->at(i))<5 ) {

      mhtx -= tauPt_PF->at(i) * cos(tauPhi_PF->at(i));
      mhty -= tauPt_PF->at(i) * sin(tauPhi_PF->at(i));
    }
  }
  
  return make_pair(mhtx,mhty);
}

float basicLoop::getMHT() {
  std::pair<float,float> mht=getJERAdjustedMHTxy();
  
  return sqrt(mht.first*mht.first + mht.second*mht.second);
}

float basicLoop::getMHTphi() {
  std::pair<float,float> mht=getJERAdjustedMHTxy();

  return atan2(mht.second,mht.first);
}

float basicLoop::getHT() {
  //use isGoodJet() to recalculate it for the specified jet type

  float ht=0;

  for (unsigned int i=0; i<loosejetPt->size(); i++) {
    if (isGoodJet( i ) ) ht+= getLooseJetPt(i);
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
    if (isGoodJet_Sync1( i ) ) ht+= getLooseJetPt(i);
  }

  return ht;
}

//the main accessor for jet pt
//automatically adjust for jes rescaling
//according to https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopLeptonPlusJets2010Systematics

//also adjust for JER according to that twiki
float basicLoop::getLooseJetPt( unsigned int ijet ) {

  //only allow one of JES or JER variations at a time

  //the normal case first
  if ( theJESType_ == kJES0 && theJERType_ == kJER0) return loosejetPt->at(ijet);
  //then JES variations
  else if (theJESType_ == kJESup && theJERType_ == kJER0) {
    float cor = getJESExtraUnc(ijet);
    float unc =  sqrt( loosejetJECUncPlus->at(ijet) * loosejetJECUncPlus->at(ijet) + cor*cor);
    return (1+unc) * loosejetPt->at(ijet) ;
  }
  else if (theJESType_ == kJESdown && theJERType_ == kJER0) {
    float cor = getJESExtraUnc(ijet);
    float unc =  sqrt( loosejetJECUncMinus->at(ijet) * loosejetJECUncMinus->at(ijet) + cor*cor);
    return (1-unc) *loosejetPt->at(ijet) ;
  }
  //now do JER variations
  else if (theJERType_ != kJER0 && theJESType_ == kJES0) {
    float genpt = loosejetGenPt->at(ijet);
    float recopt = loosejetPt->at(ijet);
    if (genpt >15 ) {
      float factor = 0;
      if (theJERType_ == kJERbias )    factor = 0.1; //hard-coded factors from top twiki
      else if (theJERType_ == kJERup)  factor = 0.2;
      else {assert(0);}

      float  deltapt = (recopt - genpt) * factor;
      float frac = (recopt+deltapt)/recopt;
      float ptscale = frac>0 ? frac : 0;
      recopt *= ptscale;
    }
    return recopt;
  }

  //if we haven't returned by now there must be something illegal going one!
  std::cout<<"Problem in getLooseJetPt with illegal combination of JES and JER settings!"<<std::endl;
  assert(0);

  return -1;
}

//apparently i need a method to get uncorrected jet pt adjusted by the JER correction...here it is
float basicLoop::getLooseJetPtUncorr( unsigned int ijet) {

  //the normal case first
  if ( theJESType_ == kJES0 && theJERType_ == kJER0) return loosejetPtUncorr->at(ijet);
  else if (theJESType_ != kJES0 ) {
    cout<<"Not implemented or desired!"<<endl;
  }
  else if (theJERType_ != kJER0 && theJESType_ == kJES0) {
    float genpt = loosejetGenPt->at(ijet);
    float recopt = loosejetPt->at(ijet);
    float uncorrpt = loosejetPtUncorr->at(ijet);
    if (genpt >15 ) {
      float factor = 0;
      if (theJERType_ == kJERbias )    factor = 0.1; //hard-coded factors from top twiki
      else if (theJERType_ == kJERup)  factor = 0.2;
      else {assert(0);}

      float  deltapt = (recopt - genpt) * factor;
      float frac = (recopt+deltapt)/recopt;
      float ptscale = frac>0 ? frac : 0;
      uncorrpt *= ptscale;
    }
    return uncorrpt;
  }

  std::cout<<"Problem in getLooseJetPtUncorr with illegal combination of JES and JER settings!"<<std::endl;
  assert(0);

  return -1;
}

//doJetID_ is a major kludge introduced in order to disable jetid for some tests
//in reality there should be a jetID enum (that also goes into the filename)
const bool doJetID_=true; //this should be *true* unless you're doing something special!
bool basicLoop::isGoodJet(unsigned int ijet) {

  if ( getLooseJetPt(ijet) <50) return false;
  if ( fabs(loosejetEta->at(ijet)) > 2.4) return false;
  if (doJetID_ && !(loosejetPassLooseID->at(ijet)) ) return false;

  return true;
}

bool basicLoop::isGoodJet30(unsigned int ijet) {

  if ( getLooseJetPt(ijet) <30) return false;
  if ( fabs(loosejetEta->at(ijet)) > 2.4) return false;
  if (doJetID_ && !(loosejetPassLooseID->at(ijet)) ) return false;

  return true;
}

//i should really just make the pT cut an argument of the function.
//but at this point i'm too nervous to change existing code
bool basicLoop::isGoodJet10(unsigned int ijet) {

  if ( getLooseJetPt(ijet) <10) return false;
  if ( fabs(loosejetEta->at(ijet)) > 2.4) return false;
  if (doJetID_ && !(loosejetPassLooseID->at(ijet)) ) return false;

  return true;
}

bool basicLoop::isGoodJetMHT(unsigned int ijet) {

  if ( getLooseJetPt(ijet) <30) return false;
  if ( fabs(loosejetEta->at(ijet)) > 5) return false;
  //no jet id for MHT

  return true;
}

bool basicLoop::isBadJet(unsigned int ijet) {
  //reverse the isGoodJet() cuts, except pT

  if ( getLooseJetPt(ijet) <50) return false;
  if ( fabs(loosejetEta->at(ijet)) > 5 || fabs(loosejetEta->at(ijet)) < 2.4) return false;
  if ( loosejetPassLooseID->at(ijet) ) return false;

  return true;
}


//get MC btag efficiency
float basicLoop::jetBTagEff(unsigned int ijet) {

  float btageff=0;
  float pt = getLooseJetPt(ijet);
  float eta = loosejetEta->at(ijet);
  int flavor = loosejetFlavor->at(ijet);

  if(isGoodJet30(ijet)){

    //double SF = 0.9;
    //double SF = 1;
    float SF_cen[5] = {0.948,0.883,0.856,0.848,0.91};   // 5 bins of pt 
    float SF_for[5] = {0.889,0.861,0.791,0.655,0.96};   // 5 bins of pt
    //double SF_unc = 1.15;                            
    //double SF_unc = 0.85;                            
    //float SF_unc = 1;
    float SFU_cen[5] = {1,1,1,1,1};
    float SFU_for[5] = {1,1,1,1,1};

    if (theBTagEffType_ == kBTagEffup) {
      SFU_cen[0] = 1.077;
      SFU_cen[1] = 1.042;
      SFU_cen[2] = 1.076;
      SFU_cen[3] = 1.034;
      SFU_cen[4] = 1.041;
      SFU_for[0] = 1.050;
      SFU_for[1] = 1.052;
      SFU_for[2] = 1.15;
      SFU_for[3] = 1.021;
      SFU_for[4] = 1.094;
    }
    else if (theBTagEffType_ == kBTagEffdown) {
      SFU_cen[0] = 0.967;
      SFU_cen[1] = 0.977;
      SFU_cen[2] = 0.958;
      SFU_cen[3] = 0.955;
      SFU_cen[4] = 0.800;
      SFU_for[0] = 0.962;
      SFU_for[1] = 0.976;
      SFU_for[2] = 0.936;
      SFU_for[3] = 0.960;
      SFU_for[4] = 0.967;     
    }
    


    if(fabs(eta)<1.4){ //"central" jets     
      if( abs(flavor) == 5){
	if( pt < 50 )        return SF_cen[0]*SFU_cen[0]*0.536914;
	else if ( pt < 75 )  return SF_cen[1]*SFU_cen[1]*0.623229;
	else if ( pt < 100 ) return SF_cen[2]*SFU_cen[2]*0.663002;
	else if ( pt < 150 ) return SF_cen[3]*SFU_cen[3]*0.671216;
	else                 return SF_cen[4]*SFU_cen[4]*0.623814;
      }
      else if (abs(flavor) == 4){
	if( pt < 50 ) return 0.123408;
	else if ( pt < 75 ) return 0.161444;
	else if ( pt < 100 ) return 0.181389;
	else if ( pt < 150 ) return 0.192003;
	else return 0.180417;
      }
      else if (abs(flavor) == 1 || abs(flavor) == 2
	       || abs(flavor) == 3 || abs(flavor) == 21){
	if( pt < 50 ) return 0.00829395;
	else if ( pt < 75 ) return 0.011271;
	else if ( pt < 100 ) return 0.013844;
	else if ( pt < 150 ) return 0.0166908;
	else return 0.0224991;
      }
    }
    else{//"forward" jets      
      if( abs(flavor) == 5){
	if( pt < 50 )        return SF_for[0]*SFU_for[0]*0.402681;
	else if ( pt < 75 )  return SF_for[1]*SFU_for[1]*0.560273;
	else if ( pt < 100 ) return SF_for[2]*SFU_for[2]*0.632681;
	else if ( pt < 150 ) return SF_for[3]*SFU_for[3]*0.662464;
	else                 return SF_for[4]*SFU_for[4]*0.66867;
      }
      else if (abs(flavor) == 4){
	if( pt < 50 ) return 0.0960317;
	else if ( pt < 75 ) return 0.144701;
	else if ( pt < 100 ) return 0.162798;
	else if ( pt < 150 ) return 0.190575;
	else return 0.186154;
      }
      else if (abs(flavor) == 1 || abs(flavor) == 2
	       || abs(flavor) == 3 || abs(flavor) == 21){
	if( pt < 50 ) return 0.00645874;
	else if ( pt < 75 ) return 0.0112103;
	else if ( pt < 100 ) return 0.0162004;
	else if ( pt < 150 ) return 0.0234471;
	else return 0.033612;
      }
    }


  }

  return btageff;
}






















void basicLoop::getSphericityJetMET(float & lambda1, float & lambda2, float & det,
				    const int jetmax, bool addMET) {
  double l1,l2,d;
  getSphericityJetMET(l1,l2,d,jetmax,addMET);

  lambda1 = l1;
  lambda2 = l2;
  det = d;

}

//event shape variables from Luke
void basicLoop::getSphericityJetMET(double & lambda1, double & lambda2, double & det,
				    const int jetmax, bool addMET) {


  TMatrixD top3Jets(2,2);
  double top3JetsScale = 0;
  
  unsigned int njets=  loosejetPt->size();
  int ngoodj=0;
  for (unsigned int i=0; i<njets ; i++) {
    if ( isGoodJet(i) ) {
      ++ngoodj;
      double phi = loosejetPhi->at(i);
      double eT = getLooseJetPt(i);
      double eX = eT*cos(phi);
      double eY = eT*sin(phi);
      if(ngoodj <= jetmax) {
	top3Jets[0][0] += eX*eX/eT;
	top3Jets[0][1] += eX*eY/eT;
	top3Jets[1][0] += eX*eY/eT;
	top3Jets[1][1] += eY*eY/eT;
	top3JetsScale += eT;
      }
      else break;
    }
  }

  if (addMET) {
    double phi = getMETphi();
    double eT = getMET();
    double eX = eT*cos(phi);
    double eY = eT*sin(phi);
  
    top3Jets[0][0] += eX*eX/eT;
    top3Jets[0][1] += eX*eY/eT;
    top3Jets[1][0] += eX*eY/eT;
    top3Jets[1][1] += eY*eY/eT;
    top3JetsScale += eT;
  }
  top3Jets*=1/top3JetsScale;
  TMatrixDEigen top3JetsEigen(top3Jets);
  lambda1 = top3JetsEigen.GetEigenValuesRe()[0];
  lambda2 = top3JetsEigen.GetEigenValuesRe()[1];
  det = top3Jets.Determinant();

}


float basicLoop::getJetInvisibleEnergyHT() {
  //scalar sum of the MC truth jet invisible energy

  //i am applying no jet cuts here!

  float totalInvisibleEnergy=0;
  for (unsigned int ij=0; ij< loosejetInvisibleEnergy->size(); ++ij ) {
    if ( loosejetInvisibleEnergy->at(ij) > 0) { //need to guard against unmatched jets
      totalInvisibleEnergy += loosejetInvisibleEnergy->at(ij);
    }
  }
  return totalInvisibleEnergy;
}

float basicLoop::getJetInvisibleEnergyMHT() {
  float x=0,y=0;

  for (unsigned int ij=0; ij< loosejetInvisibleEnergy->size(); ++ij ) {
    if ( loosejetInvisibleEnergy->at(ij) > 0) { //need to guard against unmatched jets
      x -= loosejetInvisibleEnergy->at(ij) * cos(loosejetGenPhi->at(ij));
      y -= loosejetInvisibleEnergy->at(ij) * sin(loosejetGenPhi->at(ij));
    }
  }
  return sqrt(x*x + y*y);
}

float basicLoop::getLargestJetPtRecoError(unsigned int maxjets) {
  float biggest=0;

  unsigned int loopmax = (loosejetPt->size() < maxjets) ? loosejetPt->size() : maxjets;
  //should i be applying jet cuts here?
  for (unsigned int ij=0; ij<loopmax; ++ij) {
    float residual = loosejetGenPt->at(ij) > 0 ? getLooseJetPt(ij) - loosejetGenPt->at(ij) : 0;
    if ( fabs(residual) > fabs(biggest) ) biggest = residual;
  }
  return biggest;
}

float basicLoop::getDeltaPhiMismeasuredMET(unsigned int maxjets) {

  float biggest=0;
  unsigned int ibiggest=0;

  unsigned int loopmax = (loosejetPt->size() < maxjets) ? loosejetPt->size() : maxjets;
  //should i be applying jet cuts here?
  for (unsigned int ij=0; ij<loopmax; ++ij) {
    float residual = loosejetGenPt->at(ij) > 0 ? getLooseJetPt(ij) - loosejetGenPt->at(ij) : 0;
    if ( fabs(residual) > fabs(biggest) ) { biggest = residual; ibiggest = ij; }
  }

  return getDeltaPhi(loosejetPhi->at(ibiggest),getMETphi());

}

bool basicLoop::isGoodJet_Sync1(unsigned int ijet) {

  if ( getLooseJetPt(ijet) <50) return false;
  
  if ( fabs(loosejetEta->at(ijet)) > 2.4) return false;
  
  if (loosejetEnergyFracHadronic->at(ijet) <0.1) return false;

  return true;
}

bool basicLoop::isLooseJet_Sync1(unsigned int ijet) {

  if ( getLooseJetPt(ijet) <30) return false;
  
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
      if (ngood==n) return  getLooseJetPt(i);
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

float basicLoop::bjetPtOfN(unsigned int n) {
  //updated to use Sync1 or Baseline0 cuts

  unsigned int ngood=0;
  for (unsigned int i=0; i<loosejetPt->size(); i++) {

    bool pass=false;
    if (theCutScheme_==kBaseline0) pass = isGoodJet30(i) && passSSVM(i);
    else {assert(0);}

    if (pass ) {
      ngood++;
      if (ngood==n) return  getLooseJetPt(i);
    }
  }
  return 0;
}

float basicLoop::bjetPhiOfN(unsigned int n) {
  //updated to use Sync1 or Baseline0 cuts

  unsigned int ngood=0;
  for (unsigned int i=0; i<loosejetPt->size(); i++) {

    bool pass=false;
    if (theCutScheme_==kBaseline0) pass = isGoodJet30(i) && passSSVM(i);
    else {assert(0);}

    if (pass ) {
      ngood++;
      if (ngood==n) return   loosejetPhi->at(i);
    }
  }
  return 0;
}

float basicLoop::bjetEtaOfN(unsigned int n) {
  //updated to use Sync1 or Baseline0 cuts

  unsigned int ngood=0;
  for (unsigned int i=0; i<loosejetPt->size(); i++) {

    bool pass=false;
    if (theCutScheme_==kBaseline0) pass = isGoodJet30(i) && passSSVM(i);
    else {assert(0);}

    if (pass ) {
      ngood++;
      if (ngood==n) return   loosejetEta->at(i);
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

unsigned int basicLoop::nBadJets() {
  
  unsigned int njets=0;
  for (unsigned int i=0; i<loosejetPt->size(); i++) {
    if (isBadJet(i) )   njets++;
  }
  return njets;
}

bool basicLoop::passBadJetVeto() {
  return (nBadJets() == 0);
}

bool basicLoop::passCleaning() {

  if (theCleaningType_ == kNoCleaning) return true;
  else if (theCleaningType_ == kMuonCleaning) {
    //these bools have now been simplified in the ntuple.
    //there are no longer extraneous _PF/not _PF copies
    return (passesBadPFMuonFilter && passesInconsistentMuonPFCandidateFilter);
  }
  else if (theCleaningType_ ==kMuonEcalCleaning) {
    bool passMuon = passesBadPFMuonFilter && passesInconsistentMuonPFCandidateFilter;
    bool passEcal = passEcalDeadCellCleaning();
    return passMuon && passEcal;
  }
  else {assert(0);}

  return false;
}

//look for greedy muons (subject to the limitations explained in the ntuple maker code)
bool basicLoop::passGreedyMuon() {

  for (unsigned int i=0; i<muonEoverP_PF->size(); ++i) {
    if (muonEoverP_PF->at(i) >= 1) return false;
  }

  return true;
}

//for studying the pT cut of the BadPFMuonFilter
float basicLoop::findMaxMuonPtDiff() {
  float max=0;

  for (unsigned int i=0; i<muonPtDiff_PF->size(); i++) {
    if (muonPtDiff_PF->at(i) > max) max = muonPtDiff_PF->at(i);
  }

  return max;
}

void basicLoop::findMostInconsistentPFMuon(float &mostInconsistentMuonPt, float &mostInconsistentMuonPtDiff) {
  mostInconsistentMuonPtDiff=0;
  mostInconsistentMuonPt=0;

  for (unsigned int i=0; i<muonPtRatio_PF->size(); i++) {
    float diff = fabs(muonPtRatio_PF->at(i) - 1);
    if (diff > mostInconsistentMuonPtDiff) {
      mostInconsistentMuonPt = muonPt_PF->at(i);
      mostInconsistentMuonPtDiff = muonPtRatio_PF->at(i);
    }
  }
  
}

bool basicLoop::passEcalDeadCellCleaning() {
  //the same logic as eventIsSpecificed(), just inverted

  if (ecalVetoEvents_.empty()) return true;

  jmt::eventID thisevent;
  thisevent.run = runNumber;
  thisevent.ls = lumiSection;
  thisevent.ev = eventNumber;

  if ( ecalVetoEvents_.find( thisevent) != ecalVetoEvents_.end()) {
    //    cout<<"failed BE filter: "<<runNumber<<" "<<lumiSection<<" "<<eventNumber<<endl;
    return false;
  }

  return true;
}


void basicLoop::fillWTop() {

  // cout<<" == event =="<<endl;

  double bestM2j=-9999;//, bestM2j_j1pt=0, bestM2j_j2pt=0;
  double bestM3j=-9999;//, bestM3j_j3pt=0;

  //adopting this code from Owen -- note the loop goes to the second to last jet only
  for (unsigned int j1i = 0; j1i < loosejetPt->size() -1; j1i++) {

    if ( isGoodJet30(j1i)) { //owen is using pT>30 cut

      //use exactly the same logic as Owen, to avoid bugs
      if (passSSVM(j1i) ) continue; //veto b jets

      //note how owen does the loop indexing here
      for (unsigned int j2i =j1i+1; j2i<loosejetPt->size(); j2i++) {
	if ( isGoodJet10(j2i)) { //owen is using a pT>10 cut here!

	  if (isGoodJet30(j2i) && passSSVM(j2i)) continue; //veto b jets with >30 gev

	  double m2j = calc_mNj(j1i,j2i);
	  if ( fabs(m2j- mW_) < fabs(bestM2j - mW_) ) {

	    bestM2j = m2j;
	    // bestM2j_j1pt = getLooseJetPt(j1i);
	    //bestM2j_j2pt = getLooseJetPt(j2i);

	    for ( unsigned int j3i=0; j3i<loosejetPt->size(); j3i++) {

	      if (j3i==j1i || j3i==j2i) continue;

	      if ( isGoodJet30(j3i) && passSSVM(j3i)) { //owen uses 30 GeV pT cut

		double m3j = calc_mNj(j1i,j2i,j3i);

		if ( fabs(m3j-mtop_) < fabs(bestM3j-mtop_) ) {
		  bestM3j=m3j;
		  calcCosHel(j1i,j2i,j3i); //update helicity angles
		  //bestM3j_j3pt = getLooseJetPt(j3i);
		  //		  cout<<"New best!"<<endl;
		}

	      } //is j3 good b jet
	    } //loop over j3

	  } // compare with W mass

	} //jet 2 is good
      } //j2i

    } //if jet is good

  } //j1i

  bestWMass_ = bestM2j;
  bestTopMass_ = bestM3j;

}


double basicLoop::calc_mNj( unsigned int j1i, unsigned int j2i) {
  std::vector<unsigned int> v;

  v.push_back(j1i);
  v.push_back(j2i);
  return calc_mNj(v);
}

double basicLoop::calc_mNj( unsigned int j1i, unsigned int j2i, unsigned int j3i) {
  std::vector<unsigned int> v;

  v.push_back(j1i);
  v.push_back(j2i);
  v.push_back(j3i);
  return calc_mNj(v);
}

double basicLoop::calc_mNj( std::vector<unsigned int> jNi ) {

//we could use an std::set which enforces that the elements are unique
  for (unsigned int i=0; i<jNi.size()-1; i++) {
    for (unsigned int j=i+1; j<jNi.size(); j++) {
      if (jNi.at(i) == jNi.at(j)) {
	cout<<"Problem in calc_mNj!"<<endl;
	return -1;
      }
    }
  }

  double sumE =0;
  double sumPx=0;
  double sumPy=0;
  double sumPz=0;

  for (unsigned int i=0; i<jNi.size(); i++)   {
    unsigned int j1i = jNi.at(i);
    sumE += loosejetE->at( j1i ); 

    sumPx += getLooseJetPt(j1i)*cos(loosejetPhi->at(j1i));
    sumPy += getLooseJetPt(j1i)*sin(loosejetPhi->at(j1i));
    sumPz += loosejetPz->at(j1i);
  }

   double sumP2 = sumPx*sumPx + sumPy*sumPy + sumPz*sumPz ;

   if ( (sumE*sumE) < sumP2 ) return -2. ;

   return sqrt( sumE*sumE - sumP2 );
}

float basicLoop::getLooseJetPx( unsigned int ijet ) {
  return getLooseJetPt(ijet) * cos(loosejetPhi->at(ijet));
}
float basicLoop::getLooseJetPy( unsigned int ijet ) {
  return getLooseJetPt(ijet) * sin(loosejetPhi->at(ijet));
}

//-- first two jets are from W.  third is b jet.
void basicLoop::calcCosHel( unsigned int j1i, unsigned int j2i, unsigned int j3i) {

   //bool verb(true) ;
   bool verb(false) ;
   if ( verb ) { printf( "\n" ) ; }

   if ( j1i == j2i || j1i == j3i || j2i == j3i) {
     WCosHel_ = -99;
     topCosHel_ = -99;
     return;
   }

   double   sumPx = getLooseJetPx(j1i) + getLooseJetPx(j2i) + getLooseJetPx(j3i);
   double   sumPy = getLooseJetPy(j1i) + getLooseJetPy(j2i) + getLooseJetPy(j3i);
   double   sumPz = loosejetPz->at(j1i) + loosejetPz->at(j2i) + loosejetPz->at(j3i); 

   //--- ignore quark masses.
   double j1Elab = sqrt( pow( getLooseJetPx(j1i), 2) + pow( getLooseJetPy(j1i), 2) + pow( loosejetPz->at(j1i), 2) ) ;
   double j2Elab = sqrt( pow( getLooseJetPx(j2i), 2) + pow( getLooseJetPy(j2i), 2) + pow( loosejetPz->at(j2i), 2) ) ;
   double j3Elab = sqrt( pow( getLooseJetPx(j3i), 2) + pow( getLooseJetPy(j3i), 2) + pow( loosejetPz->at(j3i), 2) ) ;

   double sumP2 = sumPx*sumPx + sumPy*sumPy + sumPz*sumPz ;
   double sumP = sqrt( sumP2 ) ;

   //--- assume the top mass when computing the energy.

   double sumE = sqrt( mtop_*mtop_ + sumP2 ) ;

   double betatop = sumP / sumE ;
   double gammatop = 1. / sqrt( 1. - betatop*betatop ) ;

   TLorentzVector j1p4lab( getLooseJetPx(j1i), getLooseJetPy(j1i), loosejetPz->at(j1i), j1Elab );
   TLorentzVector j2p4lab( getLooseJetPx(j2i), getLooseJetPy(j2i), loosejetPz->at(j2i), j2Elab );
   TLorentzVector j3p4lab( getLooseJetPx(j3i), getLooseJetPy(j3i), loosejetPz->at(j3i), j3Elab );

   if ( verb ) printf( " calc_wcoshel: betatop = %5.3f,  gammatop = %6.4f\n", betatop, gammatop ) ;

   double betatopx = sumPx / sumE ;
   double betatopy = sumPy / sumE ;
   double betatopz = sumPz / sumE ;

   TVector3 betatopvec( betatopx, betatopy, betatopz ) ;
   TVector3 negbetatopvec = -1.0 * betatopvec ;

   TLorentzVector j1p4trf( j1p4lab ) ;
   TLorentzVector j2p4trf( j2p4lab ) ;
   TLorentzVector j3p4trf( j3p4lab ) ;
   if (verb) { 
      printf(" calc_wcoshel: j1 tlorentzvector before boost:  %6.1f,  %6.1f,  %6.1f,   %6.1f\n",
          j1p4trf.Px(), j1p4trf.Py(), j1p4trf.Pz(), j1p4trf.E() ) ;
   }
   j1p4trf.Boost( negbetatopvec ) ;
   j2p4trf.Boost( negbetatopvec ) ;
   j3p4trf.Boost( negbetatopvec ) ;
   if (verb) { 
      printf(" calc_wcoshel: j1 tlorentzvector after  boost:  %6.1f,  %6.1f,  %6.1f,   %6.1f\n",
          j1p4trf.Px(), j1p4trf.Py(), j1p4trf.Pz(), j1p4trf.E() ) ;
   }


   topCosHel_ = j3p4trf.Vect().Dot( betatopvec ) / ( betatopvec.Mag() * j3p4trf.Vect().Mag() ) ;
   if ( verb ) { printf(" calc_wcoshel: top cos hel = %6.3f\n", topCosHel_ ) ; }

   //-- now boost j1 and j2 into W frame from top RF.

   TVector3 j1j2p3trf = j1p4trf.Vect() + j2p4trf.Vect() ;

   //-- Use W mass to compute energy
   double j1j2Etrf = sqrt( mW_*mW_ + j1j2p3trf.Mag2() );
   TLorentzVector j1j2p4trf( j1j2p3trf, j1j2Etrf ) ;
   TVector3 betawvec = j1j2p4trf.BoostVector() ;
   TVector3 negbetawvec = -1.0 * betawvec ;

   TLorentzVector j1p4wrf( j1p4trf ) ;
   TLorentzVector j2p4wrf( j2p4trf ) ;
   j1p4wrf.Boost( negbetawvec ) ;
   j2p4wrf.Boost( negbetawvec ) ;
   WCosHel_ = j1p4wrf.Vect().Dot( betawvec ) / ( betawvec.Mag() * j1p4wrf.Vect().Mag() ) ;
   if (verb) {
      printf(" calc_wcoshel: j1p4 in W rf: %6.1f,  %6.1f,  %6.1f,   %6.1f\n",
          j1p4wrf.Px(), j1p4wrf.Py(), j1p4wrf.Pz(), j1p4wrf.E() ) ;
      printf(" calc_wcoshel: j2p4 in W rf: %6.1f,  %6.1f,  %6.1f,   %6.1f\n",
          j2p4wrf.Px(), j2p4wrf.Py(), j2p4wrf.Pz(), j2p4wrf.E() ) ;
      TLorentzVector j1j2wrf = j1p4wrf + j2p4wrf ;
      printf(" calc_wcoshel: j1j2 mass in W rf: %7.2f\n", j1j2wrf.M() ) ;
      printf(" calc_wcoshel: w cos hel = = %6.3f\n", WCosHel_ ) ;
   }

   //   cout<<"cosHel calculator (t,W): "<<tcoshel<<" "<<wcoshel<<endl;
}

void basicLoop::fillWithJetFlavor(TH1D* hh, double w, double threshold) {

  //fill histogram hh with true jet flavor of all jets above pt threshold
  //use gen jets for threshold!
  for (unsigned int ij=0; ij<loosejetFlavor->size() ; ++ij) {
    if (loosejetGenPt->at(ij) > threshold) {
      hh->Fill( loosejetFlavor->at(ij), w);
    }
  }
}

void basicLoop::fillWithGenPDGId(TH1D* hh, double w, double threshold) {

  //fill histogram hh with gen pdg id of all jets above pt threshold
  //use gen jets for threshold!
  for (unsigned int ij=0; ij<loosejetGenParticlePDGId->size() ; ++ij) {
    if (loosejetGenPt->at(ij) > threshold) {
      hh->Fill( loosejetGenParticlePDGId->at(ij), w);
    }
  }
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

  /*
this function is now rather misnamed, since I have removed the tight jet variables

but we still use the nbSSVM, so it is critical to call this for every event!
  */

  if (theCutScheme_==kBaseline0) nbSSVM = countBJets();
  else {
    assert(0);
  }

}

//when new jet variables are added to the ntuple then need to update this!
void basicLoop::InitJets() {

  //i have never updated this for JPT jets, but I think it is irrelevant at the moment

  if (theJetType_==kCalo ) {
    tightJetIndex = tightJetIndex_calo;
    looseJetIndex = looseJetIndex_calo;
    loosejetPt = loosejetPt_calo;
    loosejetPz = loosejetPz_calo;
    loosejetPtUncorr = loosejetPtUncorr_calo;
    loosejetEt = loosejetEt_calo;
    loosejetE = loosejetE_calo;
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

    loosejetJECUncPlus = loosejetJECUncPlus_calo;
    loosejetJECUncMinus = loosejetJECUncMinus_calo;

  }
  else if (theJetType_==kPF) {
    tightJetIndex = tightJetIndex_PF;
    looseJetIndex = looseJetIndex_PF;
    loosejetPt = loosejetPt_PF;
    loosejetPz = loosejetPz_PF;
    loosejetPtUncorr = loosejetPtUncorr_PF;
    loosejetEt = loosejetEt_PF;
    loosejetE = loosejetE_PF;
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

    loosejetJECUncPlus = loosejetJECUncPlus_PF;
    loosejetJECUncMinus = loosejetJECUncMinus_PF;
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

  else if (inname.Contains("/QCD-Pt0to5-PythiaZ2/"))       return "PythiaQCD0";
  else if (inname.Contains("/QCD-Pt5to15-PythiaZ2/"))      return "PythiaQCD5";
  else if (inname.Contains("/QCD-Pt15to30-PythiaZ2/"))     return "PythiaQCD15";
  else if (inname.Contains("/QCD-Pt30to50-PythiaZ2/"))     return "PythiaQCD30";

  else if (inname.Contains("/QCD-Pt50to80-PythiaZ2/"))     return "PythiaQCD50";
  else if (inname.Contains("/QCD-Pt80to120-PythiaZ2/"))    return "PythiaQCD80";
  else if (inname.Contains("/QCD-Pt120to170-PythiaZ2/"))   return "PythiaQCD120";
  else if (inname.Contains("/QCD-Pt170to300-PythiaZ2/"))   return "PythiaQCD170";
  else if (inname.Contains("/QCD-Pt300to470-PythiaZ2/"))   return "PythiaQCD300";
  else if (inname.Contains("/QCD-Pt470to600-PythiaZ2/"))   return "PythiaQCD470";
  else if (inname.Contains("/QCD-Pt600to800-PythiaZ2/"))   return "PythiaQCD600";
  else if (inname.Contains("/QCD-Pt800to1000-PythiaZ2/"))  return "PythiaQCD800";
  else if (inname.Contains("/QCD-Pt1000to1400-PythiaZ2/")) return "PythiaQCD1000";
  else if (inname.Contains("/QCD-Pt1400to1800-PythiaZ2/")) return "PythiaQCD1400";
  else if (inname.Contains("/QCD-Pt1800toInf-PythiaZ2/"))  return "PythiaQCD1800";

  else if (inname.Contains("/QCD-Pt0to5-PythiaZ2-PU2010/"))       return "PythiaPUQCD0";
  else if (inname.Contains("/QCD-Pt5to15-PythiaZ2-PU2010/"))      return "PythiaPUQCD5";
  else if (inname.Contains("/QCD-Pt15to30-PythiaZ2-PU2010/"))     return "PythiaPUQCD15";
  else if (inname.Contains("/QCD-Pt30to50-PythiaZ2-PU2010/"))     return "PythiaPUQCD30";

  else if (inname.Contains("/QCD-Pt50to80-PythiaZ2-PU2010/"))     return "PythiaPUQCD50";
  else if (inname.Contains("/QCD-Pt80to120-PythiaZ2-PU2010/"))    return "PythiaPUQCD80";
  else if (inname.Contains("/QCD-Pt120to170-PythiaZ2-PU2010/"))   return "PythiaPUQCD120";
  else if (inname.Contains("/QCD-Pt170to300-PythiaZ2-PU2010/"))   return "PythiaPUQCD170";
  else if (inname.Contains("/QCD-Pt300to470-PythiaZ2-PU2010/"))   return "PythiaPUQCD300";
  else if (inname.Contains("/QCD-Pt470to600-PythiaZ2-PU2010/"))   return "PythiaPUQCD470";
  else if (inname.Contains("/QCD-Pt600to800-PythiaZ2-PU2010/"))   return "PythiaPUQCD600";
  else if (inname.Contains("/QCD-Pt800to1000-PythiaZ2-PU2010/"))  return "PythiaPUQCD800";
  else if (inname.Contains("/QCD-Pt1000to1400-PythiaZ2-PU2010/")) return "PythiaPUQCD1000";
  else if (inname.Contains("/QCD-Pt1400to1800-PythiaZ2-PU2010/")) return "PythiaPUQCD1400";
  else if (inname.Contains("/QCD-Pt1800toInf-PythiaZ2-PU2010/"))  return "PythiaPUQCD1800";

  else if (inname.Contains("/QCD-Pt15to3000-PythiaZ2-Flat-PU2010/"))  return "PythiaPUQCDFlat";


  else if (inname.Contains("/DATA/"))  {
    if (realDatasetNames_) {
      int lastslash=     inname.Last('/');
      TString wholepath=inname(0,lastslash);
      lastslash = wholepath.Last('/');
      return wholepath(lastslash+1,wholepath.Length());
    }
    else return "data";
  }
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

  else if (inname.Contains("/Zinvisible/"))                return 5760; //NNLO, taken from RA2 note //RA1 uses 5715
  
  else if (inname.Contains("/QCD-Pt0to5-PythiaZ2/"))       return 4.844e10;
  else if (inname.Contains("/QCD-Pt5to15-PythiaZ2/"))      return 3.675e10;
  else if (inname.Contains("/QCD-Pt15to30-PythiaZ2/"))     return 8.159e8;
  else if (inname.Contains("/QCD-Pt30to50-PythiaZ2/"))     return 5.312e7;

  else if (inname.Contains("/QCD-Pt50to80-PythiaZ2/"))     return 6.359e6;
  else if (inname.Contains("/QCD-Pt80to120-PythiaZ2/"))    return 7.843e5;
  else if (inname.Contains("/QCD-Pt120to170-PythiaZ2/"))   return 1.151e5;
  else if (inname.Contains("/QCD-Pt170to300-PythiaZ2/"))   return 2.426e4;
  else if (inname.Contains("/QCD-Pt300to470-PythiaZ2/"))   return 1.168e3;
  else if (inname.Contains("/QCD-Pt470to600-PythiaZ2/"))   return 7.022e1;
  else if (inname.Contains("/QCD-Pt600to800-PythiaZ2/"))   return 1.555e1;
  else if (inname.Contains("QCD-Pt800to1000-PythiaZ2/"))   return 1.844;
  else if (inname.Contains("QCD-Pt1000to1400-PythiaZ2/"))  return 3.321e-1;
  else if (inname.Contains("QCD-Pt1400to1800-PythiaZ2/"))  return 1.087e-2;
  else if (inname.Contains("QCD-Pt1800toInf-PythiaZ2/"))   return 3.575e-4;

  //these numbers are exactly the same as those above for the non-PU samples
  else if (inname.Contains("/QCD-Pt0to5-PythiaZ2-PU2010/"))       return 4.844e10;
  else if (inname.Contains("/QCD-Pt5to15-PythiaZ2-PU2010/"))      return 3.675e10;
  else if (inname.Contains("/QCD-Pt15to30-PythiaZ2-PU2010/"))     return 8.159e8;
  else if (inname.Contains("/QCD-Pt30to50-PythiaZ2-PU2010/"))     return 5.312e7;

  else if (inname.Contains("/QCD-Pt50to80-PythiaZ2-PU2010/"))     return 6.359e6;
  else if (inname.Contains("/QCD-Pt80to120-PythiaZ2-PU2010/"))    return 7.843e5;
  else if (inname.Contains("/QCD-Pt120to170-PythiaZ2-PU2010/"))   return 1.151e5;
  else if (inname.Contains("/QCD-Pt170to300-PythiaZ2-PU2010/"))   return 2.426e4;
  else if (inname.Contains("/QCD-Pt300to470-PythiaZ2-PU2010/"))   return 1.168e3;
  else if (inname.Contains("/QCD-Pt470to600-PythiaZ2-PU2010/"))   return 7.022e1;
  else if (inname.Contains("/QCD-Pt600to800-PythiaZ2-PU2010/"))   return 1.555e1;
  else if (inname.Contains("QCD-Pt800to1000-PythiaZ2-PU2010/"))   return 1.844;
  else if (inname.Contains("QCD-Pt1000to1400-PythiaZ2-PU2010/"))  return 3.321e-1;
  else if (inname.Contains("QCD-Pt1400to1800-PythiaZ2-PU2010/"))  return 1.087e-2;
  else if (inname.Contains("QCD-Pt1800toInf-PythiaZ2-PU2010/"))   return 3.575e-4;

  else if (inname.Contains("/QCD-Pt15to3000-PythiaZ2-Flat-PU2010/"))  return 2.213e+10;

  else if (inname.Contains("/DATA/"))                return -2;

  std::cout<<"Cannot find cross section for this sample!"<<std::endl;
  assert(0); 
  return -1;
  
}

double basicLoop::getWeight(Long64_t nentries) {
  if (isData_) return 1;

  const   double sigma = getCrossSection(findInputName());

  double w = lumi * sigma / double(nentries);

  if (findInputName().Contains("/QCD-Pt15to3000-PythiaZ2-Flat-PU2010/"))    w *= mcWeight;

  return  w;
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

void basicLoop::setJESType(JESType jestype) {
  theJESType_ = jestype;
}

void basicLoop::setJERType(JERType jertype) {
  theJERType_ = jertype;
}

void basicLoop::setMETuncType(METuncType metunctype) {
  theMETuncType_ = metunctype;
  if (theMETuncType_ != kMETunc0) {
    setJERType(kJERbias);
    cout<<"NOTE -- I am turning on the Jet Energy Resolution bias correction because you selected an Unclustered MET Correction!"<<endl;
  }
}

void basicLoop::setBTagEffType(BTagEffType btagefftype) {
  theBTagEffType_ = btagefftype;
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

void basicLoop::setCleaningType( tailCleaningType cleanuptype) {
  if (cleanuptype == kMuonEcalCleaning) assert(loadedEcalTree_);

  theCleaningType_ = cleanuptype;
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

void basicLoop::setRequiredCut(const TString cutTag) {

  bool ok=false;
  for (unsigned int i=0; i<cutTags_.size() ; i++) {
    if (cutTags_[i] == cutTag) {ok=true; break;}
  }
  if (!ok) {cout<<"Invalid cutIndex"<<endl; return;}

  requiredCut_.push_back(cutTag);

}

void basicLoop::resetIgnoredCut() {
  ignoredCut_.clear();
}

void basicLoop::resetRequiredCut() {
  requiredCut_.clear();
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

  //this implementation does not allow for a wildcard-type search (e.g. all of one lumi section)
  //can think about it....

  //this used to be 'return true'. Isn't that backwards?
  if (specifiedEvents_.empty()) return false;

  jmt::eventID thisevent;
  thisevent.run = runNumber;
  thisevent.ls = lumiSection;
  thisevent.ev = eventNumber;

  if ( specifiedEvents_.find( thisevent) != specifiedEvents_.end()) return true;

  return false;
}

void basicLoop::fillEcalVetoList(TTree* ecaltree) {

  ecalVetoEvents_.clear();

  //variables in tree
  Int_t run,ls,ev;
  std::vector<int> *cutFlowFlag=0;
  std::vector<std::string> *cutFlowStr =0;

  ecaltree->SetBranchAddress("run",&run);
  ecaltree->SetBranchAddress("event",&ev);
  ecaltree->SetBranchAddress("lumi",&ls);
  ecaltree->SetBranchAddress("cutFlowFlag",&cutFlowFlag);
  ecaltree->SetBranchAddress("cutFlowStr",&cutFlowStr);

  /*
note that this is way more complexity than we actually need. I wrote out to the tree
only the failed events, so we don't really need to check the value in the tree for each event.
But this structure is more general.
  */

  //loop over tree
  for (Long64_t i = 0; i < ecaltree->GetEntries(); ++i ) {
    ecaltree->GetEvent(i);
    for (unsigned int j = 0 ; j< cutFlowStr->size(); ++j) {
      if (TString(cutFlowStr->at(j).c_str()) == "BE") {
	if ( cutFlowFlag->at(j) == 1) { //event failed BE filter
	  jmt::eventID theEvent;
	  theEvent.run = run;
	  theEvent.ls=ls;
	  theEvent.ev = ev;
	  ecalVetoEvents_.insert(theEvent);
	}
	else {
	  cout<<"unexpected result in fillEcalVetoList! event passed BE filter!"<<endl;
	}
      }
    }
  }

  cout<<"Added "<<ecalVetoEvents_.size()<< " events to the ECAL veto list!"<<endl;
  loadedEcalTree_=true;
}

TString basicLoop::getCutDescriptionString() {

  TString cuts = CutSchemeNames_[theCutScheme_];
  if (specialCutDescription_ != "") {
    cuts += ".";
    cuts += specialCutDescription_;
    cuts += ".";
  }
  else {
    cuts += "_";
  }
  cuts += jetTypeNames_[theJetType_];
  cuts += "_";
  if (theJESType_ != kJES0) {
    cuts+= jesTypeNames_[theJESType_];
    cuts+= "_";
  }
  if (theJERType_ != kJER0) {
    cuts+= jerTypeNames_[theJERType_];
    cuts+= "_";
  }
  cuts += METTypeNames_[theMETType_];
  cuts += METRangeNames_[theMETRange_];
  cuts += "_";
  if (theMETuncType_ != kMETunc0) {
    cuts+= unclusteredMetUncNames_[theMETuncType_];
    cuts+= "_";
  }
  cuts += leptonTypeNames_[theLeptonType_];
  cuts += ne_; 
  cuts += "e";
  cuts += nmu_; 
  cuts += "mu";
  cuts += "_";
  cuts+= dpTypeNames_[theDPType_];
  if (theCleaningType_ != kNoCleaning) {
    cuts+="_";
    cuts+=tailCleaningNames_[theCleaningType_];
  }
  for (unsigned int icut=0; icut<ignoredCut_.size() ; icut++) {
    cuts+="_No";
    //it would be more robust to use .find() instead of []
    cuts += jmt::fortranize( cutNames_[ignoredCut_.at(icut)]);
  }
  for (unsigned int icut=0; icut<requiredCut_.size() ; icut++) {
    cuts+="_With";
    //it would be more robust to use .find() instead of []
    cuts += jmt::fortranize( cutNames_[requiredCut_.at(icut)]);
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
  cout<<"JES uncertainty set:  "<<jesTypeNames_[theJESType_]<<endl;
  cout<<"JER scaling set:      "<<jerTypeNames_[theJERType_]<<endl;
  cout<<"lepton type set to:   "<<leptonTypeNames_[theLeptonType_]<<endl;
  cout<<"Requiring exactly     "<<ne_<<" electrons"<<endl;
  cout<<"Requiring exactly     "<<nmu_<<" muons"<<endl;
  cout<<"MET type set to:      "<<METTypeNames_[theMETType_]<<endl;
  cout<<"MET range set to:     "<<METRangeNames_[theMETRange_]<<endl;
  cout<<"Unclustered energy:   "<<unclusteredMetUncNames_[theMETuncType_]<<endl;
  cout<<"DeltaPhi type set to: "<<dpTypeNames_[theDPType_]<<endl;
  cout<<"Requiring at least    "<<nBcut_<<" b tags"<<endl;
  cout<<"Tail cleanup set to:  "<<tailCleaningNames_[theCleaningType_]<<endl;
  for (unsigned int i = 0; i< ignoredCut_.size() ; i++) {
    cout<<"Will ignore cut:    "<<cutNames_[ignoredCut_.at(i)]<<endl;
  }
  for (unsigned int i = 0; i< requiredCut_.size() ; i++) {
    cout<<"Will require cut:   "<<cutNames_[requiredCut_.at(i)]<<endl;
  }
  
}

void basicLoop::startTimer() {
  starttime_ = new TDatime();
}

void basicLoop::checkTimer(const Long64_t ndone, const Long64_t ntotal) {
  if (ndone==0) return;
  double fracdone = double(ndone)/double(ntotal);

  TDatime timenow;
  UInt_t elapsed= timenow.Convert() - starttime_->Convert();

  cout<<int(100*fracdone) << " ["<<(1-fracdone)*double(elapsed)/fracdone<<" seconds remaining]"<<endl;

}

void basicLoop::stopTimer(const Long64_t ntotal) {
  TDatime stoptime; //default ctor is for current time
  UInt_t elapsed= stoptime.Convert() - starttime_->Convert();
  cout<<"events / time = "<<ntotal<<" / "<<elapsed<<" = "<<double(ntotal)/double(elapsed)<<" Hz"<<endl;

  delete starttime_;
  starttime_=0;
}

// ========================================== end



#endif // #ifdef basicLoop_cxx
