// -*- C++ -*-
//
// Package:    BasicTreeMaker
// Class:      BasicTreeMaker
// 
/**\class BasicTreeMaker BasicTreeMaker.cc CUSUSY/BasicTreeMaker/src/BasicTreeMaker.cc

 Description: The usual code for making a very simple tree from PATtuples
Created for studies of inclusive hadronic SUSY searches with b tags

 Implementation:
Uses STL vectors for arrays of data. This seems to work fine for analysis in bare ROOT
using a class created with MakeClass.

Developed and tested with CMSSW_3_6_2
(Now moving to CMSSW_3_6_3)

*/
//
// Original Author:  Joshua Thompson,6 R-029,+41227678914,
//         Created:  Thu Jul  8 16:33:08 CEST 2010
// $Id: BasicTreeMaker.cc,v 1.4 2010/09/10 16:16:28 joshmt Exp $
//
//


// system include files
#include <memory>
#include <iostream>
#include <vector>
#include <string>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//added by jmt
#include "FWCore/Utilities/interface/CPUTimer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"

//amazingly, this won't build if these two includes are not in this order!
#include "PhysicsTools/SelectorUtils/interface/EventSelector.h"
#include "PhysicsTools/SelectorUtils/interface/PVSelector.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "CUSUSY/Selection/interface/SUSYEventSelector.h"

//Root includes
#include "TFile.h"
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
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  //  bool passPV(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillJetInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillLeptonInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillTrackInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillMCInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillTriggerInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);

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

  std::vector<std::string> btagAlgorithmNames_;
  std::vector<std::string> triggersOfInterest_;

  SUSYEventSelector* susyCutFlow_;

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
  std::vector<bool> cutResultsDon;
  std::vector<bool> cutResults;

  std::vector<bool> passTrigger;

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
  std::vector<float> loosejetPt;
  std::vector<float> loosejetEta;
  std::vector<float> loosejetPhi;
  std::vector<int> loosejetFlavor;
  std::vector<int> loosejetGenParticlePDGId;
  std::vector<float> loosejetInvisibleEnergy;
  std::map < std::string, std::vector<float> > loosejetBTagDisc;

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

//
// constructors and destructor
//
BasicTreeMaker::BasicTreeMaker(const edm::ParameterSet& iConfig) :
  thetimer_(new edm::CPUTimer()),
  Heventcount_(0),
  tree_(0),
  infotree_(0),

  btagAlgorithmNames_(iConfig.getParameter<std::vector<std::string> >("btagAlgorithms")),
  triggersOfInterest_(iConfig.getParameter<std::vector<std::string> >("triggersOfInterest")),

  susyCutFlow_( new SUSYEventSelector( iConfig.getParameter<edm::ParameterSet>("susyBJetsSelection" ))),

  pvSelector_      (iConfig.getParameter<edm::ParameterSet>("pvSelector") ),

  //  triggerLabel_    (iConfig.getParameter<edm::InputTag>("triggerTag")),
  jetLabel_        (iConfig.getParameter<edm::InputTag>("jetTag")),
  
  caloMetLabel_    (iConfig.getParameter<edm::InputTag>("caloMetTag")),
  tcMetLabel_      (iConfig.getParameter<edm::InputTag>("tcMetTag")),

  eleLabel_        (iConfig.getParameter<edm::InputTag>("eleTag")),
  muoLabel_        (iConfig.getParameter<edm::InputTag>("muoTag")),

  //  susyTrigger_     (iConfig.getParameter<std::string>("susyTrigger")),

  jetIdLoose_      (iConfig.getParameter<edm::ParameterSet>("jetIdLoose") ),
  muonId_          (iConfig.getParameter<edm::ParameterSet>("muonId") ),
  electronId_      (iConfig.getParameter<edm::ParameterSet>("electronId") ),

  minJets_         (iConfig.getParameter<int>("minJets")),

  jetPtMin_        (iConfig.getParameter<double>("jetPtMin")), 
  jetEtaMax_       (iConfig.getParameter<double>("jetEtaMax")),
  loosejetPtMin_   (iConfig.getParameter<double>("loosejetPtMin")), 
  loosejetEtaMax_  (iConfig.getParameter<double>("loosejetEtaMax")),

  muPtMin_         (iConfig.getParameter<double>("muPtMin")), 
  muEtaMax_        (iConfig.getParameter<double>("muEtaMax")), 
  eleEtMin_        (iConfig.getParameter<double>("eleEtMin")), 
  eleEtaMax_       (iConfig.getParameter<double>("eleEtaMax")), 

  mhtMin_          (iConfig.getParameter<double>("mhtMin")),
  metMin_          (iConfig.getParameter<double>("metMin")),

  jetInfoFilled_(false),
  leptonInfoFilled_(false),
  trackInfoFilled_(false)
{

  edm::Service<TFileService> fs;
  Heventcount_ = fs->make<TH1D>( "Heventcount"  , "events processed", 1,  0, 1 );
  tree_ = new TTree("tree","tree");
  infotree_ = new TTree("infotree","infotree");
  //branch creation done in beginJob instead of here

}


BasicTreeMaker::~BasicTreeMaker()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  delete thetimer_;
 
}


//
// member functions
//

//
void
BasicTreeMaker::fillTriggerInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  using namespace edm;

  bool passTrig=false;
  
  //some code from Josh Bendavid to automagically get the trigger tag (in MC)
  //https://hypernews.cern.ch/HyperNews/CMS/get/physTools/1791/1/1/1/1/1/2.html

  //i feel like i should not have to do this for every event (to save time)
  //but I can't put it in beginJob (I think) because there is no edm::Event at that point
  //-- could implement a global so that it only runs once
  Handle<trigger::TriggerEvent> triggerEventHLT;
  iEvent.getByLabel("hltTriggerSummaryAOD", triggerEventHLT);
  
  InputTag trigResultsTag("TriggerResults","",triggerEventHLT.provenance()->processName()); 
  //std::cout<<"Will use trigger tag = "<<triggerEventHLT.provenance()->processName()<<std::endl;

  // get HLT trigger information
  Handle<TriggerResults> HLTResults;
  iEvent.getByLabel(trigResultsTag, HLTResults);
  
  //based on code copied from Don
  if (HLTResults.isValid()) {
    
    const edm::TriggerNames & triggerNames = iEvent.triggerNames(*HLTResults);
    unsigned int trigger_size = HLTResults->size();

    bool foundAGoodTrigger=false;
    //loop over the list of triggers
    for (unsigned int itrig = 0; itrig< triggersOfInterest_.size(); itrig++) {
      
      unsigned int trigger_position = triggerNames.triggerIndex(triggersOfInterest_.at(itrig));
      
      if (trigger_position < trigger_size) {
	bool passed=HLTResults->accept(trigger_position);
	passTrigger.push_back( passed );
	if (!foundAGoodTrigger) { //this must be the first valid trigger
	  //	  std::cout<<"Will store in cutflow the results for: "<<triggersOfInterest_.at(itrig)<<std::endl; //debug
	  foundAGoodTrigger=true;
	  passTrig = passed;
	}
      }
      else {	//	std::cout << "WARNING -- Trigger position not found" << std::endl;
	passTrigger.push_back(false);
      }
      
    }

    //verbose code for looking at the trigger menu
    if (false) {
      for (unsigned int i = 0; i<trigger_size; i++) {
	std::cout<<triggerNames.triggerName(i)<<"\t"<<HLTResults->accept(i)<<std::endl;
      }
    }
  }
  else {
    std::cout << "WARNING -- TriggerResults missing from this event." << std::endl;
  }
  
  cutResults.push_back(passTrig); //store main trigger result in cut flow
}

//return the number of SUSY mothers
int
BasicTreeMaker::findSUSYMaternity( const reco::Candidate & cand ) {

  int numSUSY=0;

  int numOfMothers =  cand.numberOfMothers();

  if (numOfMothers >=1 ) {
    //    if (numOfMothers >1)    std::cout<<"multiple mothers!"<<std::endl;
    for (int j = 0; j<numOfMothers ; j++) {
      int motherId = abs(cand.mother(j)->pdgId()); //get rid of minus signs!
      if ( (motherId>= 1000001 && motherId <=1000037) || (motherId>= 2000001 && motherId <=2000015)) {numSUSY++;}
      //the case where the b comes from a gluon
      //don't need to do anything in this case, but we do need to avoid entering the else{} block
      else if (motherId == 21) {/*std::cout<<"Found gluon splitting!"<<std::endl;*/}
      else { numSUSY+= findSUSYMaternity( *cand.mother(j) ); }
    }
  }
  
  return numSUSY;
}

int
BasicTreeMaker::findTopDecayMode( const reco::Candidate & cand) {
  //tested using a TTbarJets madgraph sample
  const int numOfDaughters = cand.numberOfDaughters();
  bool foundb=false, foundq=false,founde=false,foundm=false,foundt=false;
  bool found_tau_had=false,found_tau_elec=false,found_tau_muon=false;
  for (int i=0; i<numOfDaughters; ++i) {
    //std::cout<<cand.daughter(i)->pdgId()<<"\t"<<cand.daughter(i)->status()<<std::endl;

    const reco::Candidate *daughter = cand.daughter(i); //should we check if it is null?
    if ( abs(daughter->pdgId()) == 5) foundb=true;
    else if (abs(daughter->pdgId()) == 24) { //this is the W
      for (size_t j=0; j < daughter->numberOfDaughters(); ++j) {
	//seems that the daughters of the W will be status 3 decay products and a status 2 copy of the W
	if (daughter->daughter(j)->status() == 3) {
	  int wdau = abs(daughter->daughter(j)->pdgId());
	  if (wdau <= 4) foundq = true;
	  else if (wdau == 11) founde = true;
	  else if (wdau == 13) foundm = true;
	  else if (wdau == 15) {
	    foundt = true;
	    const reco::Candidate *mytau = daughter->daughter(j);
	    for (size_t k=0; k<mytau->numberOfDaughters(); k++) {
	      //	      std::cout<<mytau->daughter(k)->pdgId()<<"\t"<<mytau->daughter(k)->status()<<std::endl;
	      const reco::Candidate *mytau2 = mytau->daughter(k);
	      if (abs(mytau2->pdgId()) == 15) {
		
		for (size_t l=0; l<mytau2->numberOfDaughters(); l++) {
		  int taudau = abs(mytau2->daughter(l)->pdgId());
		  if (taudau == 11) found_tau_elec=true;
		  else if (taudau == 13) found_tau_muon=true;
		  else if (taudau ==213 || taudau==111 || taudau==113 ||taudau==211||taudau==130
			   || taudau==223 ||taudau==321 ||taudau ==323 ||taudau ==310||taudau==221) found_tau_had=true;
		  else if (taudau == 14 || taudau==16 || taudau==12) { } //do nothing for neutrinos
		  else {
		    std::cout<<"WARNING -- Unknown tau daughter found = "<<taudau<<std::endl; //debug info
		  }
		}
	      }
	      else if (abs(mytau2->pdgId()) ==22) { //photon
		//do nothing
	      }
	      else std::cout<<"WARNING -- tau daughter is not as expected! "<<mytau2->pdgId()<<"\t"<<mytau2->status()<<std::endl;
	    }
	  }
	  //skip neutrinos completely....
	}
      }
    }
  }

  if ( !foundb ) return 0; //this shouldn't happen

  //i should really put in a check that only one of these bools is set!
  if ( foundq ) return 1; //hadronic W decay
  if ( founde ) return 2; //W->e
  if ( foundm ) return 3; //W->mu
  if ( foundt ) { //W->tau
    if ((found_tau_had && found_tau_elec) || (found_tau_had && found_tau_muon) ||(found_tau_elec&&found_tau_muon))
      std::cout<<"WARNING -- something weird happened in tau decay detection!"<<std::endl;
    if (found_tau_had) return 4;   
    if (found_tau_elec) return 5;   
    if (found_tau_muon) return 6;
    //    std::cout<<"Returning tau decay status 7!"<<std::endl; //debug info
    return 7; //not yet categorized (probably hadronic)
  }
  //shouldn't get here
  return -1;
}

void
BasicTreeMaker::fillMCInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // === get the process qscale (pthat) ===

  edm::Handle<GenEventInfoProduct> genEventInfo;
  iEvent.getByLabel("generator", genEventInfo);

  qScale=genEventInfo->qScale();

  // === get some info about the decay structure in the event ===

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles",genParticles); //TODO make this configurable

  for (size_t k = 0 ; k<genParticles->size(); k++) {
    const reco::Candidate & TCand = (*genParticles)[ k ];
    if (abs(TCand.pdgId())== 5 && TCand.status()==3) { //find b quark
      int nSUSY=findSUSYMaternity(TCand);
      //if we have one b, and it has nSUSY greater than 1, we still want to just count the b once
      if (  nSUSY>0 ) SUSY_nb++;
    }
    else if (abs(TCand.pdgId())== 6 && TCand.status()==3) { //find t quark
      int topcode = findTopDecayMode(TCand);
      topDecayCode.push_back(topcode);
    }

  }


}

void
BasicTreeMaker::fillTrackInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  if (trackInfoFilled_) return;

  edm::Handle<reco::TrackCollection> generalTracks;
  iEvent.getByLabel("generalTracks", generalTracks); //TODO make this name configurable in py

  //no need to bother with sorting tracks
  for (reco::TrackCollection::const_iterator trit=generalTracks->begin(); trit!=generalTracks->end() ; trit++) {

    if (trit->pt() > 5 && fabs(trit->eta()) < 5) { //cuts TODO make them configurable
      trackPt.push_back(trit->pt());
      trackEta.push_back(trit->eta());
      trackPhi.push_back(trit->phi());
    }
  }

  trackInfoFilled_=true;
}

void
BasicTreeMaker::fillLeptonInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  if (leptonInfoFilled_) return;

  //to get the cut flow in the right (arbitrary) order, need to do muons first

  int nmuon=0;

  edm::Handle<edm::View<pat::Muon> > muonHandle;
  iEvent.getByLabel(muoLabel_,muonHandle);
  const edm::View<pat::Muon> & muons = *muonHandle;
  for (edm::View<pat::Muon>::const_iterator imuon = muons.begin(); imuon!=muons.end(); ++imuon) {
    
    if ( !imuon->muonID("GlobalMuonPromptTight")) continue;
    // Muon veto
    bool passMuon = muonId_(*imuon,iEvent);
    if (  imuon->pt() > muPtMin_ && fabs(imuon->eta()) < muEtaMax_ && passMuon )       nmuon++;
  }
  cutResults.push_back( nmuon == 0 );


  edm::Handle<edm::View<pat::Electron> > electronHandle;
  iEvent.getByLabel(eleLabel_,electronHandle);
  const edm::View<pat::Electron> & electrons = *electronHandle;

  int nelectron=0;
  for (edm::View<pat::Electron>::const_iterator ielectron = electrons.begin(); ielectron!=electrons.end(); ++ielectron) {
    if ( ielectron->et() > eleEtMin_ && fabs(ielectron->eta()) < eleEtaMax_ && 
	 electronId_(*ielectron) &&
	 ielectron->electronID( "eidLoose" ) > 0  ) {
      nelectron++;
    }
    
  }
  cutResults.push_back( nelectron == 0 );
  
  leptonInfoFilled_=true;
}

void
BasicTreeMaker::fillJetInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //  std::cout<<" == fillJetInfo =="<<std::endl; //debug

  //we're gonna do MET in here too

  if (jetInfoFilled_) return;

  edm::Handle<edm::View<pat::Jet> > jetHandle;
  iEvent.getByLabel(jetLabel_,jetHandle);
  const edm::View<pat::Jet> & jets = *jetHandle;

  //freya sort using et. don sorts with pt. stick with pt
  std::vector<float> jetPtVector;
  for(edm::View<pat::Jet>::const_iterator ijet = jets.begin(); ijet != jets.end(); ++ijet) {
    jetPtVector.push_back( ijet->pt() );
  }

  std::vector<int> sortedJetIndices = IndexSorter< std::vector<float> >(jetPtVector,true)();
  
  //variables for use in the jet loop
  pat::strbitset ret1 = jetIdLoose_.getBitTemplate();
  double MHTx = 0, MHTy = 0;
  //now loop over sorted jets
  for (size_t ii = 0 ; ii< sortedJetIndices.size(); ++ii) {
    //get the index of the ii'th jet
    int jj = sortedJetIndices[ii];
    //now pull out that jet
    const pat::Jet & jet = jets[jj];
    
    //MHT is calculated using loose jet cuts
    //HT is calculated using regular jet cuts
    bool passJetID = false;
    if ( jet.isCaloJet() ) passJetID = jetIdLoose_(jet, ret1);
    
    
    // === make cuts ===
    bool passLooseCuts = (jet.pt() > loosejetPtMin_) && (fabs(jet.eta())<loosejetEtaMax_) && passJetID;
    bool passTightCuts = (jet.pt() > jetPtMin_) && (fabs(jet.eta())<jetEtaMax_) && passJetID;
    
    //let's enforce that the tight cuts are always a subset of the loose cuts
    if (passTightCuts && !passLooseCuts) assert(0);
    
    //now fill ntuple!
    if (passLooseCuts) {
      MHTx -= jet.px();
      MHTy -= jet.py();
      loosejetPt.push_back( jet.pt() );
      loosejetEta.push_back(jet.eta());
      loosejetPhi.push_back(jet.phi());
      loosejetFlavor.push_back( jet.partonFlavour() );

      if ( jet.genJet() != 0 && jet.genParticle() != 0) {
	
	loosejetGenParticlePDGId.push_back( jet.genParticle()->pdgId() );
	loosejetInvisibleEnergy.push_back( jet.genJet()->invisibleEnergy());
	//	std::cout<<jet.partonFlavour()<<"\t"<<jet.genParticle()->pdgId()<<"\t"<<jet.genJet()->invisibleEnergy()<<std::endl;
	//std::cout<<"\t\t"<<jet.genJet()->getGenConstituents().size()<<std::endl;
	
	//We could add a function here to drill down through the jet constituents
	//and sum up the missing *transverse* energy belonging to e.g. neutrinos
	//i'm not going to implement this now
	
	// 	for (unsigned int ijet=0; ijet<jet.genJet()->getGenConstituents().size(); ijet++) {
	// 	  std::cout<<"\t\t"<<jet.genJet()->getGenConstituents().at(ijet)->pdgId()
	// 		   <<"\t"<<jet.genJet()->getGenConstituents().at(ijet)->energy()
	// 		   <<"\t"<<jet.genJet()->getGenConstituents().at(ijet)->et()<<std::endl;
	// 	}
      }
      else {
	//	std::cout<<jet.partonFlavour()<<"\tNo gen match!"<<std::endl;
	loosejetGenParticlePDGId.push_back( -99 );
	loosejetInvisibleEnergy.push_back( -99 );
      }
      
      
      for (unsigned int ib=0; ib<btagAlgorithmNames_.size(); ib++) {
	loosejetBTagDisc[btagAlgorithmNames_[ib]].push_back( jet.bDiscriminator(btagAlgorithmNames_[ib]) );
      }
    }
    if (passTightCuts) {
      //allow us to reference things only stored for loose jets
      looseJetIndex.push_back( loosejetPt.size() - 1);
      
      HT += jet.pt();
      jetPt.push_back( jet.pt() );
      jetEta.push_back(jet.eta());
      jetPhi.push_back(jet.phi());
      jetFlavor.push_back( jet.partonFlavour() );
      
      float bdisc = 0;
      bool foundSSVM=false;
      for (unsigned int ib=0; ib<btagAlgorithmNames_.size(); ib++) {
	bdisc = jet.bDiscriminator(btagAlgorithmNames_[ib]);
	jetBTagDisc[btagAlgorithmNames_[ib]].push_back( bdisc );
	//look especially for SSV
	//when the tagger doesn't exist, it seems to return -1000. So this 'if' should fire only once per jet
	//but to be safe we'll use this foundSSVM variable
	if ( ((btagAlgorithmNames_[ib] == "simpleSecondaryVertexBJetTags") && (bdisc > -999))
	     || ((btagAlgorithmNames_[ib] == "simpleSecondaryVertexHighEffBJetTags") && (bdisc >-999))) {
	  //	  std::cout<<"Found SSV algorithm! jet = "<<ii<<std::endl; //debug
	  if (bdisc>=1.74) foundSSVM=true;  //hard-coded SSV medium cut
	}
      }
      if (foundSSVM) nbSSVM++; //hard-coded SSV medium cut
      
    } //end of tight cuts block
  } //end of loop over jets

  MHT = sqrt( MHTx*MHTx + MHTy*MHTy);
  MHTphi = atan2(MHTy,MHTx);
  int nloosejets = (int) loosejetPhi.size();
  if (nloosejets >0)  DeltaPhi_JetMHT1 = fabs(reco::deltaPhi( MHTphi, loosejetPhi.at(0)));
  if (nloosejets >1)  DeltaPhi_JetMHT2 = fabs(reco::deltaPhi( MHTphi, loosejetPhi.at(1)));
  if (nloosejets >2)  DeltaPhi_JetMHT3 = fabs(reco::deltaPhi( MHTphi, loosejetPhi.at(2)));

  //fill cut flow info
  // == are there at least 3 good jets?
  int ngoodjets = (int) jetPt.size();
  cutResults.push_back( ngoodjets >= minJets_ );
  // == check pt of lead 3 jets
  if (ngoodjets >0)  cutResults.push_back( jetPt.at(0) >= 50); else cutResults.push_back(false); //FIXME ... 50 should come from python
  if (ngoodjets >1)  cutResults.push_back( jetPt.at(1) >= 50); else cutResults.push_back(false); //FIXME ... 50 should come from python
  if (ngoodjets >2)  cutResults.push_back( jetPt.at(2) >= 50); else cutResults.push_back(false); //FIXME ... 50 should come from python

  //next is the HT cut //FIXME hardcoded to RA2 value
  cutResults.push_back(HT >= 300);

  //next step in cut flow is MET
  edm::Handle<edm::View<pat::MET> > metHandle;
  iEvent.getByLabel(caloMetLabel_,metHandle);
  MET = metHandle->front().et();
  METphi = metHandle->front().phi();

  edm::Handle<edm::View<pat::MET> > tcmetHandle; //i dunno if I really need a new one of these; can't hurt
  iEvent.getByLabel(tcMetLabel_,tcmetHandle);
  tcMET = tcmetHandle->front().et();
  tcMETphi = tcmetHandle->front().phi();

  cutResults.push_back(MET >= metMin_); //note that we are cutting on caloMET
  //the cut will have to be redone in order to use tcMET

  cutResults.push_back(MHT >= mhtMin_);

  jetInfoFilled_=true;
}

// bool
// BasicTreeMaker::passPV(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// {
//   //Primary vertex
//   edm::Handle<std::vector<reco::Vertex> > pvHandle;
//   iEvent.getByLabel(pvLabel_,pvHandle);

//   //TODO
// }

// ------------ method called to for each event  ------------
void
BasicTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  
  //  cout<<" ==== "<<endl;
  Heventcount_->Fill(0.5);
  
  pat::strbitset  ret =susyCutFlow_->getBitTemplate();

  //reset tree variables
  ret.set(false);
  cutResults.clear();
  cutResultsDon.clear();
  passTrigger.clear();
  jetInfoFilled_=false;
  leptonInfoFilled_=false;
  trackInfoFilled_=false;
  looseJetIndex.clear();
  jetPt.clear();
  jetEta.clear();
  jetPhi.clear();
  jetFlavor.clear();

  loosejetPt.clear();
  loosejetEta.clear();
  loosejetPhi.clear();
  loosejetFlavor.clear();
  loosejetGenParticlePDGId.clear();
  loosejetInvisibleEnergy.clear();

  for (unsigned int ib=0; ib<btagAlgorithmNames_.size(); ib++) {
    jetBTagDisc[btagAlgorithmNames_[ib]].clear();
    loosejetBTagDisc[btagAlgorithmNames_[ib]].clear();
  }

  nbSSVM=0;
  HT=0;
  MHT=0;
  MET=0;
  METphi=0;
  MHTphi=0;
  tcMET=0;
  tcMETphi=0;
  DeltaPhi_JetMHT1 = -99;
  DeltaPhi_JetMHT2 = -99;
  DeltaPhi_JetMHT3 = -99;
  trackPt.clear();
  trackEta.clear();
  trackPhi.clear();
  SUSY_nb=0;
  qScale=0;
  topDecayCode.clear();
  //start analyzin'

  //run don's selector and copy to tree
  (*susyCutFlow_)(iEvent, ret);
  for (vector<string>::const_iterator iname = cutNames.begin(); iname!=cutNames.end() ; ++iname) {
    cutResultsDon.push_back( ret[*iname]);

    //i'm not taking the trigger info from here anymore either....this code is rather silly now
    if (*iname == string("Inclusive") ) {
      cutResults.push_back( ret[*iname]);
      //cout<<"Found cut name = "<<*iname<<endl;
    }
  }

  //trigger
  fillTriggerInfo(iEvent,iSetup);

  //primary vertex
  //  cutResults.push_back( passPV(iEvent,iSetup)); //TODO
  cutResults.push_back( pvSelector_( iEvent ) );

  //fills jets cuts, as well as MET, and the values of DeltaPhi
  fillJetInfo( iEvent,iSetup);

  //fill muon, electron
  fillLeptonInfo(iEvent,iSetup);

  //the next step in the cut flow is DeltaPhi
  //FIXME remove hard-coding
  bool dphi1cut = DeltaPhi_JetMHT1 > 0.3;//dphi1Min_;
  bool dphi2cut = DeltaPhi_JetMHT2 > 0.5;//dphi2Min_;
  bool dphi3cut = DeltaPhi_JetMHT3 > 0.3;//dphi3Min_;
  cutResults.push_back( dphi1cut && dphi2cut && dphi3cut);

  //finish the cut flow by counting b jets
  cutResults.push_back(nbSSVM>=1);
  cutResults.push_back(nbSSVM>=2);
  cutResults.push_back(nbSSVM>=3);



  fillTrackInfo(iEvent,iSetup);
  fillMCInfo(iEvent,iSetup);

  tree_->Fill();

// #ifdef THIS_IS_AN_EVENT_EXAMPLE
//    Handle<ExampleData> pIn;
//    iEvent.getByLabel("example",pIn);
// #endif
   
// #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
//    ESHandle<SetupData> pSetup;
//    iSetup.get<SetupRecord>().get(pSetup);
// #endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
BasicTreeMaker::beginJob()
{
  //nfail_=0;

  //== create tree branches ==

  //tree that is filled only once per job
  infotree_->Branch("SUSY_cutNames",&cutNames);
  infotree_->Branch("btagAlgorithms",&btagAlgorithmNames_);
  infotree_->Branch("triggerList",&triggersOfInterest_);

  //tree that is filled for each event
  tree_->Branch("SUSY_cutResults",&cutResultsDon);
  tree_->Branch("cutResults",&cutResults);
  tree_->Branch("passTrigger",&passTrigger); //stores results for triggers listed in triggerList
  //jet branches
  tree_->Branch("looseJetIndex",&looseJetIndex);
  tree_->Branch("jetPt",&jetPt);
  tree_->Branch("jetEta",&jetEta);
  tree_->Branch("jetPhi",&jetPhi);
  tree_->Branch("jetFlavor",&jetFlavor);
  //ROOT does not allow us to put a map of vectors into a tree
  //so we have to go with this approach
  for (unsigned int ib=0; ib<btagAlgorithmNames_.size(); ib++) {
    std::string bname = "jetBTagDisc_";
    bname += btagAlgorithmNames_[ib];
    tree_->Branch(bname.c_str(),&jetBTagDisc[btagAlgorithmNames_[ib]]);
  }

  //loose jet branches
  tree_->Branch("loosejetPt",&loosejetPt);
  tree_->Branch("loosejetEta",&loosejetEta);
  tree_->Branch("loosejetPhi",&loosejetPhi);
  tree_->Branch("loosejetFlavor",&loosejetFlavor);
  tree_->Branch("loosejetGenParticlePDGId",&loosejetGenParticlePDGId);
  tree_->Branch("loosejetInvisibleEnergy",&loosejetInvisibleEnergy);
  //ROOT does not allow us to put a map of vectors into a tree
  //so we have to go with this approach
  for (unsigned int ib=0; ib<btagAlgorithmNames_.size(); ib++) {
    std::string bname = "loosejetBTagDisc_";
    bname += btagAlgorithmNames_[ib];
    tree_->Branch(bname.c_str(),&loosejetBTagDisc[btagAlgorithmNames_[ib]]);
  }

  tree_->Branch("nbSSVM",&nbSSVM,"nbSSVM/I");
  tree_->Branch("HT",&HT,"HT/F");
  tree_->Branch("MHT",&MHT,"MHT/F");
  tree_->Branch("MHTphi",&MHTphi,"MHTphi/F");
  tree_->Branch("DeltaPhi_JetMHT1",&DeltaPhi_JetMHT1,"DeltaPhi_JetMHT1/F");
  tree_->Branch("DeltaPhi_JetMHT2",&DeltaPhi_JetMHT2,"DeltaPhi_JetMHT2/F");
  tree_->Branch("DeltaPhi_JetMHT3",&DeltaPhi_JetMHT3,"DeltaPhi_JetMHT3/F");
  tree_->Branch("MET",&MET,"MET/F");
  tree_->Branch("METphi",&METphi,"METphi/F");
  tree_->Branch("tcMET",&tcMET,"tcMET/F");
  tree_->Branch("tcMETphi",&tcMETphi,"tcMETphi/F");
  tree_->Branch("trackPt",&trackPt);
  tree_->Branch("trackEta",&trackEta);
  tree_->Branch("trackPhi",&trackPhi);

  tree_->Branch("SUSY_nb",&SUSY_nb,"SUSY_nb/I");
  tree_->Branch("qScale",&qScale,"qScale/F");
  tree_->Branch("topDecayCode",&topDecayCode);

  //copied from Don's code
  //someday i need to find a better way to do it
  cutNames.push_back( "Inclusive"        );
  cutNames.push_back( "Trigger"          );
  cutNames.push_back( "PV"               );
  cutNames.push_back( ">=3 Jets"         );
  cutNames.push_back( "First Jet"        );
  cutNames.push_back( "Second Jet"       );
  cutNames.push_back( "Third Jet"        );
  cutNames.push_back( "HT Cut"           );
  cutNames.push_back( "MET Cut"          );
  cutNames.push_back( "MHT Cut"          );
  cutNames.push_back( "Muon Veto"        );
  cutNames.push_back( "Electron Veto"    );
  cutNames.push_back( "DeltaPhi Cut"     );
  cutNames.push_back( ">=1 B-Tagged Jet" );
  cutNames.push_back( ">=2 B-Tagged Jets");
  cutNames.push_back( ">=3 B-Tagged Jets");

  thetimer_->start();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BasicTreeMaker::endJob() {
  using namespace std;

  thetimer_->stop();

  infotree_->Fill();

  susyCutFlow_->print(std::cout);

  //  cout<<"nfail = "<<nfail_<<endl;

  cout<<"Job took [CPU/Wall] seconds: "<<thetimer_->cpuTime()<<" / "<<thetimer_->realTime()<<endl;
  cout<<"events/sec = "<<Heventcount_->GetEntries()<<" / "<<thetimer_->realTime()<<" = "<<Heventcount_->GetEntries() /thetimer_->realTime()<<endl;  

}

//define this as a plug-in
DEFINE_FWK_MODULE(BasicTreeMaker);
