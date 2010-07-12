// -*- C++ -*-
//
// Package:    BasicTreeMaker
// Class:      BasicTreeMaker
// 
/**\class BasicTreeMaker BasicTreeMaker.cc CUSUSY/BasicTreeMaker/src/BasicTreeMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Joshua Thompson,6 R-029,+41227678914,
//         Created:  Thu Jul  8 16:33:08 CEST 2010
// $Id$
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

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

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

  int  findSUSYMaternity( const reco::Candidate & cand );

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
  SUSYEventSelector* susyCutFlow_;

  PVSelector                           pvSelector_;

  //configuration strings (largely copied from don's and freya's code)
  edm::InputTag jetLabel_;
  edm::InputTag eleLabel_;
  edm::InputTag muoLabel_;

  //  edm::InputTag               metLabel_; //TODO
  JetIDSelectionFunctor                jetIdLoose_;
  MuonVPlusJetsIDSelectionFunctor      muonId_;
  ElectronVPlusJetsIDSelectionFunctor  electronId_;


  //cut values
  double jetPtMin_ ;
  double jetEtaMax_;  
  double loosejetPtMin_ ;
  double loosejetEtaMax_;  
 
  double muPtMin_  ;
  double muEtaMax_ ;
  double eleEtMin_ ;
  double eleEtaMax_;
 
  //bookkeeping
  bool jetInfoFilled_;
  bool leptonInfoFilled_;
  bool trackInfoFilled_;

  //define variables for the tree
  std::vector<std::string> cutNames;
  std::vector<bool> cutResultsDon;
  std::vector<bool> cutResults;
  //jet info
  std::vector<float> jetPt;
  std::vector<float> jetEta;
  std::vector<float> jetPhi;
  std::vector<float> jetBTagSSV; //this is probably not the most elegant or general soln, but let's just do it
  int nbSSVM;
  float HT;
  float MHT;
  float MHTphi;
  float DeltaPhi_JetMHT1;
  float DeltaPhi_JetMHT2;
  float DeltaPhi_JetMHT3;
  //MET info
  float MET;
  float METphi;
  //NB -- as long as we only care about 'good' (tighter) jets, then 
  //we can calculate the DeltaPhi between jets and MET, as well as minDeltaPhi, later

  //track info
  std::vector<float> trackPt;
  std::vector<float> trackEta;
  std::vector<float> trackPhi;

  //MC info
  int SUSY_nb;

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
  susyCutFlow_( new SUSYEventSelector( iConfig.getParameter<edm::ParameterSet>("susyBJetsSelection" ))),

  pvSelector_      (iConfig.getParameter<edm::ParameterSet>("pvSelector") ),

  jetLabel_        (iConfig.getParameter<edm::InputTag>("jetTag")),
  eleLabel_        (iConfig.getParameter<edm::InputTag>("eleTag")),
  muoLabel_        (iConfig.getParameter<edm::InputTag>("muoTag")),

  jetIdLoose_      (iConfig.getParameter<edm::ParameterSet>("jetIdLoose") ),
  muonId_          (iConfig.getParameter<edm::ParameterSet>("muonId") ),
  electronId_      (iConfig.getParameter<edm::ParameterSet>("electronId") ),

  jetPtMin_        (iConfig.getParameter<double>("jetPtMin")), 
  jetEtaMax_       (iConfig.getParameter<double>("jetEtaMax")),
  loosejetPtMin_   (iConfig.getParameter<double>("loosejetPtMin")), 
  loosejetEtaMax_  (iConfig.getParameter<double>("loosejetEtaMax")),

  muPtMin_         (iConfig.getParameter<double>("muPtMin")), 
  muEtaMax_        (iConfig.getParameter<double>("muEtaMax")), 
  eleEtMin_        (iConfig.getParameter<double>("eleEtMin")), 
  eleEtaMax_       (iConfig.getParameter<double>("eleEtaMax")), 

  jetInfoFilled_(false),
  leptonInfoFilled_(false),
  trackInfoFilled_(false)
{

  edm::Service<TFileService> fs;
  Heventcount_ = fs->make<TH1D>( "Heventcount"  , "events processed", 1,  0, 1 );
  tree_ = new TTree("tree","tree");
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

void
BasicTreeMaker::fillMCInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles",genParticles); //TODO make this configurable

  for (size_t k = 0 ; k<genParticles->size(); k++) {
    const reco::Candidate & TCand = (*genParticles)[ k ];
    if (abs(TCand.pdgId())== 5 && TCand.status()==3) { //find b quark
      int nSUSY=findSUSYMaternity(TCand);
      //if we have one b, and it has nSUSY greater than 1, we still want to just count the b once
      if (  nSUSY>0 ) SUSY_nb++;
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
  std::vector<float> loosejetPhi;
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

    if (passLooseCuts) {
      MHTx -= jet.px();
      MHTy -= jet.py();
    }
    if (passTightCuts) {
      //now fill ntuple!
      HT += jet.pt();
      jetPt.push_back( jet.pt() );
      jetEta.push_back(jet.eta());
      jetPhi.push_back(jet.phi());
      float bdisc = jet.bDiscriminator("simpleSecondaryVertexBJetTags");
      jetBTagSSV.push_back( bdisc );
      if (bdisc>=1.74) nbSSVM++;
    }
    if (passLooseCuts || passTightCuts) {
      loosejetPhi.push_back(jet.phi());
    }

  }

  MHT = sqrt( MHTx*MHTx + MHTy*MHTy);
  MHTphi = atan2(MHTy,MHTx);
  int nloosejets = (int) loosejetPhi.size();
  if (nloosejets >0)  DeltaPhi_JetMHT1 = fabs(reco::deltaPhi( MHTphi, loosejetPhi.at(0)));
  if (nloosejets >1)  DeltaPhi_JetMHT2 = fabs(reco::deltaPhi( MHTphi, loosejetPhi.at(1)));
  if (nloosejets >2)  DeltaPhi_JetMHT3 = fabs(reco::deltaPhi( MHTphi, loosejetPhi.at(2)));

  //fill cut flow info
  // == are there at least 3 good jets?
  int ngoodjets = (int) jetPt.size();
  cutResults.push_back( ngoodjets >= 3 ); //FIXME this value (and those below) should be passed in from python
  // == check pt of lead 3 jets
  if (ngoodjets >0)  cutResults.push_back( jetPt.at(0) >= 50); else cutResults.push_back(false);
  if (ngoodjets >1)  cutResults.push_back( jetPt.at(1) >= 50); else cutResults.push_back(false);
  if (ngoodjets >2)  cutResults.push_back( jetPt.at(2) >= 50); else cutResults.push_back(false);

  //next is the HT cut //FIXME hardcoded to RA2 value
  cutResults.push_back(HT >= 300);

  //next step in cut flow is MET
  edm::Handle<edm::View<pat::MET> > metHandle;
  //  iEvent.getByLabel(metLabel_,metHandle);
  iEvent.getByLabel("patMETsAK5Calo",metHandle); //FIXME dear god, this is awful
  MET = metHandle->front().et();
  METphi = metHandle->front().phi();
  cutResults.push_back(MET >= 150); //FIXME

  cutResults.push_back(MHT >= 150); //FIXME

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
  jetInfoFilled_=false;
  leptonInfoFilled_=false;
  trackInfoFilled_=false;
  jetPt.clear();
  jetEta.clear();
  jetPhi.clear();
  jetBTagSSV.clear();
  nbSSVM=0;
  HT=0;
  MHT=0;
  MET=0;
  METphi=0;
  MHTphi=0;
  DeltaPhi_JetMHT1 = -99;
  DeltaPhi_JetMHT2 = -99;
  DeltaPhi_JetMHT3 = -99;
  trackPt.clear();
  trackEta.clear();
  trackPhi.clear();
  SUSY_nb=0;
  //start analyzin'

  //run don's selector and copy to tree
  (*susyCutFlow_)(iEvent, ret);
  for (vector<string>::const_iterator iname = cutNames.begin(); iname!=cutNames.end() ; ++iname) {
    cutResultsDon.push_back( ret[*iname]);

    //we trust the trigger info...and unfortunately nothing else
    if (*iname == string("Inclusive") || *iname == string("Trigger")) {
      cutResults.push_back( ret[*iname]);
      //cout<<"Found cut name = "<<*iname<<endl;
    }
  }

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

  /* i realize that this cross check is not going to work unless i add a lot of intelligence
  if ( cutResultsDon.size() != cutResults.size() ) {
    nfail_++;
    std::cout<<"cutResults size = "<<cutResults.size()
	     <<" ; cutResultsDon size = "<<cutResultsDon.size()<<std::endl;
  }
  for (size_t icut=0; icut<cutResultsDon.size(); icut++) {
    if (!cutResultsDon.at(icut)) break;
    if (cutResultsDon.at(icut) != cutResults.at(icut)) {
      nfail_++;
      cout<<"Fail: "<<cutNames.at(icut)<<" "<<cutResultsDon.at(icut)
	  <<cutResults.at(icut)<<endl;
    }
  }
  */

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

  tree_->Branch("SUSY_cutNames",&cutNames);
  tree_->Branch("SUSY_cutResults",&cutResultsDon);
  tree_->Branch("cutResults",&cutResults);
  //jet branches
  tree_->Branch("jetPt",&jetPt);
  tree_->Branch("jetEta",&jetEta);
  tree_->Branch("jetPhi",&jetPhi);
  tree_->Branch("jetBTagSSV",&jetBTagSSV);
  tree_->Branch("nbSSVM",&nbSSVM,"nbSSVM/I");
  tree_->Branch("HT",&HT,"HT/F");
  tree_->Branch("MHT",&MHT,"MHT/F");
  tree_->Branch("MHTphi",&MHTphi,"MHTphi/F");
  tree_->Branch("DeltaPhi_JetMHT1",&DeltaPhi_JetMHT1,"DeltaPhi_JetMHT1/F");
  tree_->Branch("DeltaPhi_JetMHT2",&DeltaPhi_JetMHT2,"DeltaPhi_JetMHT2/F");
  tree_->Branch("DeltaPhi_JetMHT3",&DeltaPhi_JetMHT3,"DeltaPhi_JetMHT3/F");
  tree_->Branch("MET",&MET,"MET/F");
  tree_->Branch("METphi",&METphi,"METphi/F");
  tree_->Branch("trackPt",&trackPt);
  tree_->Branch("trackEta",&trackEta);
  tree_->Branch("trackPhi",&trackPhi);

  tree_->Branch("SUSY_nb",&SUSY_nb,"SUSY_nb/I");

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

  susyCutFlow_->print(std::cout);

  //  cout<<"nfail = "<<nfail_<<endl;

  cout<<"Job took [CPU/Wall] seconds: "<<thetimer_->cpuTime()<<" / "<<thetimer_->realTime()<<endl;
  cout<<"events/sec = "<<Heventcount_->GetEntries()<<" / "<<thetimer_->realTime()<<" = "<<Heventcount_->GetEntries() /thetimer_->realTime()<<endl;  

}

//define this as a plug-in
DEFINE_FWK_MODULE(BasicTreeMaker);
