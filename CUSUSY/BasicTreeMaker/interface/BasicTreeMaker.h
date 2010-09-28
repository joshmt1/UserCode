// -*- C++ -*-
//
// Package:    BasicTreeMaker
// Class:      BasicTreeMaker
// 

//see cc file for more description

/*
to do list:
probably for synchronization purposes, we need even more info in the trees
so that cut flow can be customized:

eg:
-- lepton ID info
-- PV info
*/

//
// Original Author:  Joshua Thompson,6 R-029,+41227678914,
//         Created:  Thu Jul  8 16:33:08 CEST 2010
// $Id: BasicTreeMaker.h,v 1.1 2010/09/28 07:50:28 joshmt Exp $
//
//

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h" //test

//ROOT includes
#include "TTree.h"
#include "TH1.h"

//
// class declaration
//

class BasicTreeMaker : public edm::EDAnalyzer {
public:
  explicit BasicTreeMaker(const edm::ParameterSet&);
  ~BasicTreeMaker();
  
  
private:
  virtual void beginJob() ;
  virtual void beginRun(const edm::Run & iRun, const edm::EventSetup & iSetup);

  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  void findHLTProcessName(const edm::Event & iEvent);

  //  bool passPV(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillJetInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillLeptonInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillTrackInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillMCInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillTriggerInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillBTagDBInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  int  findSUSYMaternity( const reco::Candidate & cand );
  int  findTopDecayMode( const reco::Candidate & cand );

  //a clever piece of code stolen from Freya
  template <class C>
  struct IndexSorter {
    IndexSorter (const C& values, bool decreasing = true) : 
      values_(values), decrease_(decreasing) {}
    std::vector<int> operator() () const {
      std::vector<int> result;
      result.reserve(values_.size());
      for ( size_t i=0; i<values_.size(); ++i )  result.push_back(i);
      sort(result.begin(),result.end(),*this);
      return result;
    }
    bool operator() (int a, int b) {
      if ( decrease_ )
	return values_[a]>values_[b];
      else
	return values_[a]<values_[b];
    }
    const C& values_;
    bool decrease_;
  };

  // ----------member data ---------------------------
  edm::CPUTimer* thetimer_;
  TH1D* Heventcount_;
  TTree* tree_;
  TTree* infotree_;

  bool isMC_;
  bool doPrescale_;

  std::vector<std::string> btagAlgorithmNames_;
  std::vector<std::string> triggersOfInterest_;

  //  SUSYEventSelector* susyCutFlow_;

  PVSelector                           pvSelector_;

  //configuration strings (largely copied from don's and freya's code)
  //  edm::InputTag triggerLabel_;
  edm::InputTag jetLabel_;
  edm::InputTag caloMetLabel_;
  edm::InputTag tcMetLabel_;

  edm::InputTag eleLabel_;
  edm::InputTag muoLabel_;

  std::string susyTrigger_;

  JetIDSelectionFunctor                jetIdLoose_;
  MuonVPlusJetsIDSelectionFunctor      muonId_;
  ElectronVPlusJetsIDSelectionFunctor  electronId_;

  HLTConfigProvider hltConfig_; //test
  std::string processName_; // process name of (HLT) process for which to get HLT configuration

  //cut values
  int minJets_; //min number of Tight jets

  double jetPtMin_ ;
  double jetEtaMax_;  
  double loosejetPtMin_ ;
  double loosejetEtaMax_;  

  double muPtMin_  ;
  double muEtaMax_ ;
  double eleEtMin_ ;
  double eleEtaMax_;

  double mhtMin_;
  double metMin_;
 
  //bookkeeping
  bool jetInfoFilled_;
  bool leptonInfoFilled_;
  bool trackInfoFilled_;

  // ====== define variables for the tree ======
  std::vector<std::string> cutNames;
  //  std::vector<bool> cutResultsDon;
  std::vector<bool> cutResults;

  std::vector<bool> passTrigger;
  std::vector<unsigned int> hltPrescale;

  //tight jet info
  std::vector<int> looseJetIndex; //map from tight jet list to loose jet list
  std::vector<float> jetPt;
  std::vector<float> jetEta;
  std::vector<float> jetPhi;
  std::vector<int> jetFlavor;
  std::map < std::string, std::vector<float> > jetBTagDisc; 

  //this structure duplicates info between the tight and loose lists
  //i am starting to fix this using the looseJetIndex.
  //for convenience I will keep some duplication for now

  //loose jet info
  std::vector<int> verylooseJetIndex; //map from loose jet to very loose jet list
  std::vector<float> loosejetPt;
  std::vector<float> loosejetEt;
  std::vector<float> loosejetEta;
  std::vector<float> loosejetPhi;
  std::vector<int> loosejetFlavor;
  std::vector<int> loosejetGenParticlePDGId;
  std::vector<float> loosejetInvisibleEnergy;
  std::map < std::string, std::vector<float> > loosejetBTagDisc;

  //realized that I don't want even the 'loose' offline cuts for these jets
  std::vector<float> veryloosejetPtUncorr; //uncorrected jet info
  std::vector<float> veryloosejetEtaUncorr;
  std::vector<float> veryloosejetPhiUncorr;

  int nbSSVM; //FIXME this is still in need of a more sophisticated approach (for now this works...)

  float HT;
  float MHT;
  float MHTphi;
  float DeltaPhi_JetMHT1;
  float DeltaPhi_JetMHT2;
  float DeltaPhi_JetMHT3;
  //MET info
  float MET;
  float METphi;
  float tcMET;
  float tcMETphi;

  //track info
  std::vector<float> trackPt;
  std::vector<float> trackEta;
  std::vector<float> trackPhi;

  //MC info
  int SUSY_nb;
  float qScale;
  std::vector<int> topDecayCode;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
