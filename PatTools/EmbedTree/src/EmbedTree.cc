// -*- C++ -*-
//
// Package:    EmbedTree
// Class:      EmbedTree
// 
/**\class EmbedTree EmbedTree.cc PatTools/EmbedTree/src/EmbedTree.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Joshua Thompson
//         Created:  Tue Nov  6 10:04:04 CST 2012
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//my includes
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/METReco/interface/MET.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/Math/interface/deltaPhi.h"

//ROOT includes
#include "TTree.h"

//
// class declaration
//

class EmbedTree : public edm::EDAnalyzer {
   public:
      explicit EmbedTree(const edm::ParameterSet&);
      ~EmbedTree();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  void fillTtbarDecayCode(const edm::Event&, const edm::EventSetup&);
  void fillIsoTracks(const edm::Event&, const edm::EventSetup&);
  void fillLeptons(const edm::Event&, const edm::EventSetup&);
  void fillPreembedInfo(const edm::Event&, const edm::EventSetup&);

  void clearTreeVariables();
  int translateTopCode2PdgId(const int code);


      // ----------member data ---------------------------

  // -- config info
  edm::InputTag jetSrc_, metSrc_;
  edm::InputTag ttbarDecaySrc_;
  edm::InputTag muonSrc_,electronSrc_;
  edm::InputTag mindrSrc_;
  edm::InputTag originalMuonSrc_;

  // ntuple stuff

  TTree* tree_;

  //jet counters
  int njets70;
  int njets50;
  int njets30;

  //TODO b tag info [can't until that is working in embedding]
  //TODO top tagging, HepTopTagger
  //TODO MT2
  //TODO mT(b,MET)

  //for embedded events only, quantities for the original muon
  float minDRmuonJet;
  float seedmuonpt,seedmuoneta,seedmuonphi;

  //MET, jet info
  float DeltaPhiMetJet1;
  float DeltaPhiMetJet2;
  float DeltaPhiMetJet3;

  //  float DeltaPhiMetIsotrack1; //can't do this until I have phi of isotrack
  float DeltaPhiMetMuon1,DeltaPhiMetElectron1;

  //MET
  float MET;
  float METphi;

  //isolated tracks
  float isotrackpt1;//,isotrackphi1, isotracketa1; (no dice here...would need to modify producer)

 //leptons
  float electronpt1,electroneta1,electronphi1;
  float muonpt1,muoneta1,muonphi1;


  //mc truth, for unembedded MC only
  int ttbarDecayCode;
  float genLeptonPt1,genLeptonPt2;
  float genLeptonEta1,genLeptonEta2;
  int genLeptonFlavor1,genLeptonFlavor2; //provides more info that the (legacy) ttbarDecayCode

};

void EmbedTree::clearTreeVariables() {

  //to be incremented -- must be zero
  njets70=0;
  njets50=0;
  njets30=0;

  //generic bogus values
  ttbarDecayCode=-1;
  genLeptonPt1=-1; genLeptonPt2=-1;
  genLeptonEta1=-99; genLeptonEta2=-99;
  genLeptonFlavor1=-1; genLeptonFlavor2=-1; 

  MET=-99;
  METphi=-99;

  isotrackpt1=-99; //isotrackphi1=-99; isotracketa1=-99;

  electronpt1=-99; electroneta1=-99; electronphi1=-99;
  muonpt1=-99; muoneta1=-99; muonphi1=-99;

  DeltaPhiMetJet1=-99; DeltaPhiMetJet2=-99; DeltaPhiMetJet3=-99;
  //  DeltaPhiMetIsotrack1=-99;
  DeltaPhiMetMuon1=-99;
  DeltaPhiMetElectron1=-99;
  minDRmuonJet=-99;
  seedmuonpt=-99; seedmuoneta=-99; seedmuonphi=-99;

}

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
EmbedTree::EmbedTree(const edm::ParameterSet& iConfig) :
  jetSrc_(iConfig.getParameter<edm::InputTag>("JetSource")),
  metSrc_(iConfig.getParameter<edm::InputTag>("MetSource")),
  ttbarDecaySrc_(iConfig.getParameter<edm::InputTag>("TtbarDecaySource")),
  muonSrc_(iConfig.getParameter<edm::InputTag>("MuonSource")),
  electronSrc_(iConfig.getParameter<edm::InputTag>("ElectronSource")),
  mindrSrc_(iConfig.getParameter<edm::InputTag>("MinDRSource")),
  originalMuonSrc_(iConfig.getParameter<edm::InputTag>("OriginalMuonSource")),
  tree_(0)
{
   //now do what ever initialization is needed

   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  //branch creation for tree done in beginJob instead of here
  tree_=  fs->make<TTree>("reducedTree","reducedTree");


}


EmbedTree::~EmbedTree()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
EmbedTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  clearTreeVariables();

  //important that MET is filled first
  edm::Handle<edm::View<reco::MET> > mets;
  iEvent.getByLabel(metSrc_, mets);
  MET = mets->front().pt();
  METphi = mets->front().phi();

  //do not use pat::Jet
  edm::Handle<edm::View<reco::Jet> > jets;
  iEvent.getByLabel(jetSrc_,jets);

  for (edm::View<reco::Jet>::const_iterator jet = jets->begin(); jet != jets->end(); ++jet) {

    if ( fabs(jet->eta()) < 2.4 ) {
      float pt = jet->pt();
      if ( pt >= 30) {
	++njets30;
	if (njets30     ==1) DeltaPhiMetJet1 = std::abs(reco::deltaPhi(jet->phi(),METphi));
	else if (njets30==2) DeltaPhiMetJet2 = std::abs(reco::deltaPhi(jet->phi(),METphi));
	else if (njets30==3) DeltaPhiMetJet3 = std::abs(reco::deltaPhi(jet->phi(),METphi));
      }
      if (pt>=50) 	  ++njets50;
      if (pt>=70) 	  ++njets70;
    }
  }

  fillTtbarDecayCode(iEvent, iSetup);

  fillIsoTracks(iEvent,iSetup);

  fillLeptons(iEvent,iSetup);

  fillPreembedInfo(iEvent,iSetup);

  //fill tree
  tree_->Fill();
}

void EmbedTree::fillPreembedInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  if ( std::string(mindrSrc_.label()).empty())   minDRmuonJet = -1;
  else {
    edm::Handle<float > mdr;
    iEvent.getByLabel(mindrSrc_.label(),"minDeltaRToJet",mdr);
    minDRmuonJet = *mdr;
  }

  if (std::string(originalMuonSrc_.label()).empty()) {
    //don't really need to do anything. the variables are reset in clearTreeVariables()
  }
  else {
    edm::Handle<edm::View<reco::PFCandidate> > mulist;
    iEvent.getByLabel(originalMuonSrc_.label(),"pfMu",mulist);

    assert( mulist->size() == 1);
    for (edm::View<reco::PFCandidate>::const_iterator themu = mulist->begin(); themu != mulist->end(); ++themu) {
      seedmuonpt = themu->pt();
      seedmuoneta = themu->eta();
      seedmuonphi = themu->phi();
    }
  }

}

void EmbedTree::fillIsoTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //for now the only thing we store is the pT of the leading isolated track

  edm::Handle<std::vector<float> > cand_dzpv;
  iEvent.getByLabel("trackIsolationMaker","pfcandsdzpv",cand_dzpv);

  edm::Handle<std::vector<float> > cand_pt;
  iEvent.getByLabel("trackIsolationMaker","pfcandspt",cand_pt);

  edm::Handle<std::vector<float> > cand_trkiso;
  iEvent.getByLabel("trackIsolationMaker","pfcandstrkiso",cand_trkiso);

  edm::Handle<std::vector<int> > cand_chg;
  iEvent.getByLabel("trackIsolationMaker","pfcandschg",cand_chg);

  size_t n = cand_dzpv->size();
  for ( size_t ii = 0; ii<n; ii++) {
    
    if ( cand_chg->at(ii) != 0  && (fabs(cand_dzpv->at(ii)) < 0.05)) { //charged and compatible with PV

      //No pT cut! always store lead iso track
      //if (cand_pt->at(ii) > 10) { //pT >10
	
	if (cand_trkiso->at(ii)/cand_pt->at(ii) <0.10) { 
	  if ( cand_pt->at(ii) > isotrackpt1) {
	    isotrackpt1 = cand_pt->at(ii);
	  }
	}
	//      }
    }
  }

}

int EmbedTree::translateTopCode2PdgId(const int code) {


  //undo my old scheme and return the well-known pdg id

  if (code==2) return 11; //e
  else if (code==3) return 13;//mu
  else if (code>=4) return 15;//tau
  
  return 0;

}

void EmbedTree::fillTtbarDecayCode(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  if ( std::string(ttbarDecaySrc_.label()).empty()) ttbarDecayCode = -1;
  else {
    edm::Handle<int > code;
    iEvent.getByLabel(ttbarDecaySrc_.label(),"ttbarDecayCode",code);
    ttbarDecayCode = *code;

    //also get the generator-level lepton pT and eta
    edm::Handle<std::vector<float> > genleppt;
    edm::Handle<std::vector<float> > genlepeta;
    edm::Handle<std::vector<int> > genlepflavor;
    iEvent.getByLabel(ttbarDecaySrc_.label(),"lepGenPt",genleppt);
    iEvent.getByLabel(ttbarDecaySrc_.label(),"lepGenEta",genlepeta);
    iEvent.getByLabel(ttbarDecaySrc_.label(),"topDecayCode",genlepflavor);

    //sort by pT, so that the highest pT lepton will be in the *1 variables
    int i=0; int j=1;
    if ( (*genleppt)[0] > (*genleppt)[1] ) {
      i=0;
      j=1;
    }
    else {
      i=1;
      j=0;
    }

    genLeptonPt1 = (*genleppt)[i];
    genLeptonPt2 = (*genleppt)[j];

    genLeptonEta1 = (*genlepeta)[i];
    genLeptonEta2 = (*genlepeta)[j];

    genLeptonFlavor1 = translateTopCode2PdgId((*genlepflavor)[i]);
    genLeptonFlavor2 = translateTopCode2PdgId((*genlepflavor)[j]);
  }

}


void EmbedTree::fillLeptons(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //fill info on highest pt leptons

  //i'll bet a templated function would serve well here, but let's keep it simple for now

  edm::Handle<edm::View<reco::GsfElectron> > electrons;
  iEvent.getByLabel(electronSrc_, electrons);
  for (edm::View<reco::GsfElectron>::const_iterator ie = electrons->begin(); ie != electrons->end(); ++ie) {
    if (   ie->et() > electronpt1) {
      electronpt1 = ie->et();
      electronphi1 = ie->phi();
      electroneta1 = ie->eta();
      DeltaPhiMetElectron1 = std::abs(reco::deltaPhi(electronphi1,METphi));
    }
  }


  edm::Handle<edm::View<reco::Muon> > muons;
  iEvent.getByLabel(muonSrc_, muons);
  for (edm::View<reco::Muon >::const_iterator im = muons->begin(); im != muons->end(); ++im) {
    if (   im->pt() > muonpt1) {
      muonpt1 = im->pt();
      muonphi1 = im->phi();
      muoneta1 = im->eta();
      DeltaPhiMetMuon1 = std::abs(reco::deltaPhi(muonphi1,METphi));
    }
  }


}


// ------------ method called once each job just before starting event loop  ------------

void 
EmbedTree::beginJob()
{

  tree_->Branch("njets70",&njets70,"njets70/I");
  tree_->Branch("njets50",&njets50,"njets50/I");
  tree_->Branch("njets30",&njets30,"njets30/I");

  tree_->Branch("MET",&MET,"MET/F");
  tree_->Branch("METphi",&METphi,"METphi/F");

  tree_->Branch("DeltaPhiMetJet1",&DeltaPhiMetJet1,"DeltaPhiMetJet1/F");
  tree_->Branch("DeltaPhiMetJet2",&DeltaPhiMetJet2,"DeltaPhiMetJet2/F");
  tree_->Branch("DeltaPhiMetJet3",&DeltaPhiMetJet3,"DeltaPhiMetJet3/F");

  tree_->Branch("DeltaPhiMetMuon1",&DeltaPhiMetMuon1,"DeltaPhiMetMuon1/F");
  tree_->Branch("DeltaPhiMetElectron1",&DeltaPhiMetElectron1,"DeltaPhiMetElectron1/F");

  tree_->Branch("isotrackpt1",&isotrackpt1,"isotrackpt1/F");

  tree_->Branch("electronpt1",&electronpt1,"electronpt1/F");
  tree_->Branch("electronphi1",&electronphi1,"electronphi1/F");
  tree_->Branch("electroneta1",&electroneta1,"electroneta1/F");

  tree_->Branch("muonpt1",&muonpt1,"muonpt1/F");
  tree_->Branch("muonphi1",&muonphi1,"muonphi1/F");
  tree_->Branch("muoneta1",&muoneta1,"muoneta1/F");

  tree_->Branch("seedmuonpt",&seedmuonpt,"seedmuonpt/F");
  tree_->Branch("seedmuonphi",&seedmuonphi,"seedmuonphi/F");
  tree_->Branch("seedmuoneta",&seedmuoneta,"seedmuoneta/F");

  tree_->Branch("minDRmuonJet",&minDRmuonJet,"minDRmuonJet/F");


  tree_->Branch("ttbarDecayCode",&ttbarDecayCode,"ttbarDecayCode/I");
  tree_->Branch("genLeptonFlavor1",&genLeptonFlavor1,"genLeptonFlavor1/I");
  tree_->Branch("genLeptonPt1",&genLeptonPt1,"genLeptonPt1/F");
  tree_->Branch("genLeptonEta1",&genLeptonEta1,"genLeptonEta1/F");

  tree_->Branch("genLeptonFlavor2",&genLeptonFlavor2,"genLeptonFlavor2/I");
  tree_->Branch("genLeptonPt2",&genLeptonPt2,"genLeptonPt2/F");
  tree_->Branch("genLeptonEta2",&genLeptonEta2,"genLeptonEta2/F");


}

// ------------ method called once each job just after ending the event loop  ------------
void 
EmbedTree::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
EmbedTree::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
EmbedTree::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
EmbedTree::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
EmbedTree::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EmbedTree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EmbedTree);
