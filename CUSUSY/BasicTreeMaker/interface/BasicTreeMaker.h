/// -*- C++ -*-
//
// Package:    BasicTreeMaker
// Class:      BasicTreeMaker
// 

//see cc file for more description

//
// Original Author:  Joshua Thompson,6 R-029,+41227678914,
//         Created:  Thu Jul  8 16:33:08 CEST 2010
// $Id: BasicTreeMaker.h,v 1.18 2011/03/03 16:37:20 joshmt Exp $
//
//

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h" //test

//Sezen & Harrison's code for PDF uncertainties
#include "CUSUSY/BasicTreeMaker/src/PDFFactory.h"

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

  void fillPVInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillJetInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup, unsigned int jetIndex);
  void fillMetInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillLeptonInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup, unsigned int il);
  void fillTauInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup, unsigned int il);
  void fillTrackInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillMCInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillTriggerInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillBTagDBInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  int  findSUSYMaternity( const reco::Candidate & cand );
  int  findTopDecayMode( const reco::Candidate & cand );

  bool passJetId(const pat::Jet & jet);

  //Code from RA2 for filtering fake MHT with muons:

  bool badPFMuonFilter(const edm::Event& iEvent, edm::InputTag pfCandSource, edm::InputTag muonSource, double maxPtDiff = 100, bool doPtDiff = true, bool doPJCut = false, bool debug = false);
  bool inconsistentMuonPFCandidateFilter(const edm::Event& iEvent, edm::InputTag muonSource, double ptMin = 100, double maxPTDiff = 0.1, bool verbose = false);

  //End Code from RA2 for filtering fake MHT with muons

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

  // Instantiate the PDFFactory class:
  PDFFactory pdffactory_;

  bool isMC_;
  bool doPrescale_;

  std::vector<std::string> btagAlgorithmNames_;
  std::vector<std::string> tauidAlgorithmNames_;
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
  std::vector<std::string> tauAlgorithmNames_; //the real collection names

  std::vector<std::string> eleSelected_; //name of a list of selected electrons
  std::vector<std::string> muonSelected_; //name of a list of selected muons

  edm::InputTag pfCandSrc_;

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

  //bookkeeping
  //  bool jetInfoFilled_;
  //  bool leptonInfoFilled_;
  bool trackInfoFilled_;

  //for matching of particles
  float tolerance_;

  // ====== define variables for the tree ======
  ULong64_t runNumber;
  ULong64_t eventNumber;
  ULong64_t lumiSection;

  std::vector<bool> passTrigger;
  std::vector<unsigned int> hltPrescale;

  //primary vertex info
  bool pv_pass;
  std::vector<bool> pv_isFake;
  std::vector<float>  pv_z;
  std::vector<float>  pv_ndof;
  std::vector<float>  pv_chi2;
  std::vector<float>  pv_rho;

  //beamspot info
  float bsx;
  float bsy;
  float bsz;

  std::vector<float>  pdfWeights; //for PDF uncertainties

  std::map< std::string, std::vector<bool> > muonIsRA2;
  std::map< std::string, std::vector<bool> > muonIsGlobalMuon;
  std::map< std::string, std::vector<bool> > muonIsGlobalMuonPromptTight;
  std::map< std::string, std::vector<float> > muonPt;
  std::map< std::string, std::vector<float> > muonEta;
  std::map< std::string, std::vector<float> > muonPhi;
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

 
  //beginning RA2 fake MHT muon reducing filter info
  std::map< std::string, bool > passesBadPFMuonFilter;
  std::map< std::string, bool > passesInconsistentMuonPFCandidateFilter;
  //ending RA2 fake MHT muon reducing filter info


  //electron info for *all* electrons (!)
  std::map< std::string, std::vector<bool> > eleIsRA2;
  std::map< std::string, std::vector<float> > eleEt;
  //  std::map< std::string, std::vector<float> > elePt; //confirmed that pT and eT are the same!
  std::map< std::string, std::vector<float> > eleEta;
  std::map< std::string, std::vector<float> > elePhi;
  std::map< std::string, std::vector<float> > eleTrackIso;
  std::map< std::string, std::vector<float> > eleEcalIso;
  std::map< std::string, std::vector<float> > eleHcalIso;

  std::map< std::string, std::vector<float> > eledB;

  std::map< std::string, std::vector<float> > eleVtx_z;

  std::map< std::string, std::vector<float> > eleIDLoose;
  std::map< std::string, std::vector<float> > eleIDRobustTight;
  std::map< std::string, std::vector<bool> > elePassID;

  //tau info all *all* taus
  std::map< std::string, std::vector<float> > tauPt;
  std::map< std::string, std::vector<float> > tauEta;
  std::map< std::string, std::vector<float> > tauPhi;
  std::map< std::string, std::vector<float> > tauTaNC;
  std::map< std::string,  std::map < std::string, std::vector<bool> > > tauID;
  //againstElectron againstMuon byIsolation byTaNC byTaNCfrHalfPercent byTaNCfrQuarterPercent 

  /*
all jet quantities must because a map<string, whatever> where the
string is the jetAlgorithmTag
  */

  //tight jet info
  std::map< std::string,  std::vector<int> > tightJetIndex; //map from tight jet list to loose jet list

  //loose jet info (no longer appying any jet id here!)
  std::map< std::string,  std::vector<int> > looseJetIndex; //map from loose jet to very loose jet list
  std::map< std::string,  std::vector<float> > loosejetPt;
  std::map< std::string,  std::vector<float> > loosejetPtUncorr; //uncorrected
  std::map< std::string,  std::vector<float> > loosejetEt;
  std::map< std::string,  std::vector<float> > loosejetEta;
  std::map< std::string,  std::vector<float> > loosejetPhi;

  std::map< std::string,  std::vector<float> > loosejetPz;
  std::map< std::string,  std::vector<float> > loosejetE;

  std::map< std::string,  std::vector<bool> > loosejetPassLooseID;
  std::map< std::string,  std::vector<bool> > loosejetPassTightID;
  std::map< std::string,  std::vector<float> > loosejetEnergyFracHadronic;
  std::map< std::string,  std::vector<int> > loosejetNTracks;
  std::map< std::string,  std::vector<int> > loosejetNSV;
  std::map< std::string,  std::vector<float> > loosejetSVUnWeightedMass;
  std::map< std::string,  std::vector<float> > loosejetSVWeightedMass;
  std::map< std::string,  std::vector<float> > loosejetSVUnWeightedLifetime;
  std::map< std::string,  std::vector<float> > loosejetSVWeightedLifetime;
  std::map< std::string,  std::vector<float> > loosejetSVUnWeightedCosTheta;
  std::map< std::string,  std::vector<float> > loosejetSVWeightedCosTheta;
  std::map< std::string,  std::vector<float> > loosejetJECUncPlus;
  std::map< std::string,  std::vector<float> > loosejetJECUncMinus;

  std::map< std::string,  std::vector<int> > loosejetFlavor;
  std::map< std::string,  std::vector<int> > loosejetGenPt;
  std::map< std::string,  std::vector<int> > loosejetGenPhi;
  std::map< std::string,  std::vector<int> > loosejetGenEta;
  std::map< std::string,  std::vector<int> > loosejetGenParticlePDGId;
  std::map< std::string,  std::vector<float> > loosejetInvisibleEnergy;
  std::map< std::string,  std::map < std::string, std::vector<float> > > loosejetBTagDisc;


  //MET info
  std::map< std::string, float> MET;
  std::map< std::string, float> METphi;
  std::map< std::string, float> METsig; //met significance
  //mc met info
  std::map< std::string, float> GenMET;
  std::map< std::string, float> GenMETphi;

  //track info
  std::vector<float> trackPt;
  std::vector<float> trackEta;
  std::vector<float> trackPhi;

  //MC info
  int SUSY_nb;
  double qScale;
  double mcWeight;
  std::vector<int> topDecayCode;
  int ZDecayMode;
  int flavorHistory;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
