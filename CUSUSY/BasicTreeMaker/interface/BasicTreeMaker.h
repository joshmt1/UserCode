// -*- C++ -*-
//
// Package:    BasicTreeMaker
// Class:      BasicTreeMaker
// 

//see cc file for more description

//
// Original Author:  Joshua Thompson,6 R-029,+41227678914,
//         Created:  Thu Jul  8 16:33:08 CEST 2010
// $Id: BasicTreeMaker.h,v 1.6 2010/10/21 22:13:36 joshmt Exp $
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

  void resetTreeVariables() ;

  void findHLTProcessName(const edm::Event & iEvent);
  void fillShortNames();

  void fillPVInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillJetInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup, unsigned int jetIndex);
  void fillMetInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillLeptonInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup, unsigned int il);
  void fillTrackInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillMCInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillTriggerInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillBTagDBInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  int  findSUSYMaternity( const reco::Candidate & cand );
  int  findTopDecayMode( const reco::Candidate & cand );

  bool passJetId(const pat::Jet & jet);

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
  edm::InputTag pvLabel_;
  //  edm::InputTag jetLabel_;

  std::vector<std::string> jetAlgorithmNames_; //the real collection names
  std::vector<std::string> jetAlgorithmTags_; //id for the ntuple branches

  std::vector<std::string> metAlgorithmNames_; //the real collection names
  std::vector<std::string> metAlgorithmTags_; //id for the ntuple branches

  std::vector<std::string> eleAlgorithmNames_; //the real collection names
  std::vector<std::string> muonAlgorithmNames_; //the real collection names

  JetIDSelectionFunctor                jetIdLoose_;
  JetIDSelectionFunctor                jetIdTight_;
  PFJetIDSelectionFunctor              PFjetIdLoose_;
  PFJetIDSelectionFunctor              PFjetIdTight_;
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
  //  bool jetInfoFilled_;
  //  bool leptonInfoFilled_;
  bool trackInfoFilled_;

  // ====== define variables for the tree ======
  ULong64_t runNumber;
  ULong64_t eventNumber;
  ULong64_t lumiSection;

  std::vector<std::string> cutNames;
  std::vector<bool> cutResults;

  std::vector<bool> passTrigger;
  std::vector<unsigned int> hltPrescale;
  int SUSYtriggerIndex;

  //primary vertex info
  std::vector<bool> pv_isFake;
  std::vector<float>  pv_z;
  std::vector<float>  pv_ndof;
  std::vector<float>  pv_chi2;
  std::vector<float>  pv_rho;

  //beamspot info
  float bsx;
  float bsy;
  float bsz;

  std::map< std::string, std::vector<bool> > muonIsGlobalMuonPromptTight;
  //muon info for all *AllGlobalMuons*
  std::map< std::string, std::vector<float> > muonPt;
  std::map< std::string, std::vector<float> > muonEta;
  std::map< std::string, std::vector<float> > muonTrackIso;
  std::map< std::string, std::vector<float> > muonEcalIso;
  std::map< std::string, std::vector<float> > muonHcalIso;

  std::map< std::string, std::vector<float> > muonChi2;
  std::map< std::string, std::vector<float> > muonNdof;
  std::map< std::string, std::vector<float> > muonNhits;
  std::map< std::string, std::vector<float> > muonTrackd0;
  std::map< std::string, std::vector<float> > muonTrackPhi;
  std::map< std::string, std::vector<bool> > muonPassID;

  std::map< std::string, std::vector<float> > muonVtx_z;

  std::map< std::string, std::vector<float> > muonHcalVeto;
  std::map< std::string, std::vector<float> > muonEcalVeto;

  std::map< std::string, int > nMuons; //good muons (passing pt and eta cuts)

  //electron info for *all* electrons (!)
  std::map< std::string, std::vector<float> > eleEt;
  std::map< std::string, std::vector<float> > eleEta;
  std::map< std::string, std::vector<float> > eleTrackIso;
  std::map< std::string, std::vector<float> > eleEcalIso;
  std::map< std::string, std::vector<float> > eleHcalIso;

  std::map< std::string, std::vector<float> > eledB;

  std::map< std::string, std::vector<float> > eleVtx_z;

  std::map< std::string, std::vector<float> > eleIDLoose;
  std::map< std::string, std::vector<float> > eleIDRobustTight;
  std::map< std::string, std::vector<bool> > elePassID;

  std::map< std::string, int > nElectrons; //good electrons

  /*
all jet quantities must because a map<string, whatever> where the
string is the jetAlgorithmTag
  */

  //tight jet info
  std::map< std::string,  std::vector<int> > tightJetIndex; //map from tight jet list to loose jet list

  //loose jet info (no longer appying any jet id here!)
  std::map< std::string,  std::vector<int> > looseJetIndex; //map from loose jet to very loose jet list
  std::map< std::string,  std::vector<float> > loosejetPt;
  std::map< std::string,  std::vector<float> > loosejetEt;
  std::map< std::string,  std::vector<float> > loosejetEta;
  std::map< std::string,  std::vector<float> > loosejetPhi;
  std::map< std::string,  std::vector<bool> > loosejetPassLooseID;
  std::map< std::string,  std::vector<bool> > loosejetPassTightID;
  std::map< std::string,  std::vector<float> > loosejetEnergyFracHadronic;

  std::map< std::string,  std::vector<int> > loosejetFlavor;
  std::map< std::string,  std::vector<int> > loosejetGenParticlePDGId;
  std::map< std::string,  std::vector<float> > loosejetInvisibleEnergy;
  std::map< std::string,  std::map < std::string, std::vector<float> > > loosejetBTagDisc;

  std::map< std::string,  std::vector<float> > badjetPt;
  std::map< std::string,  std::vector<float> > badjetEta;
  std::map< std::string,  std::vector<float> > badjetPhi;

  //here we store just caloJets (uncorrected)
  //realized that I don't want even the 'loose' offline cuts for these jets
  std::vector<float> veryloosejetPtUncorr; //uncorrected jet info
  std::vector<float> veryloosejetEtaUncorr;
  std::vector<float> veryloosejetPhiUncorr;

  int nbSSVM;

  float HT;
  float MHT;
  float MHTphi;
  float DeltaPhi_JetMHT1;
  float DeltaPhi_JetMHT2;
  float DeltaPhi_JetMHT3;
  //MET info
  std::map< std::string, float> MET;
  std::map< std::string, float> METphi;
  std::map< std::string, float> METsig; //met significance

  //track info
  std::vector<float> trackPt;
  std::vector<float> trackEta;
  std::vector<float> trackPhi;

  //MC info
  int SUSY_nb;
  float qScale;
  std::vector<int> topDecayCode;
  int flavorHistory;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
