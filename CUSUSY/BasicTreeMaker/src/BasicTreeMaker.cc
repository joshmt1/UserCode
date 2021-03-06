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

Originally developed and tested with CMSSW_3_6_2
Have also used in in 384.
Latest incarnation is for 386 and 387

  Recipe is kept here:
https://wiki.lepp.cornell.edu/lepp/bin/view/CMS/JMTBasicNtuples
*/
//
// Original Author:  Joshua Thompson,6 R-029,+41227678914,
//         Created:  Thu Jul  8 16:33:08 CEST 2010
// $Id: BasicTreeMaker.cc,v 1.35 2011/03/25 22:25:45 joshmt Exp $
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
#include "FWCore/Framework/interface/ESHandle.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

//amazingly, this won't build if these two includes are not in this order!
#include "PhysicsTools/SelectorUtils/interface/EventSelector.h"
#include "PhysicsTools/SelectorUtils/interface/PVSelector.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "RecoBTag/PerformanceDB/interface/BtagPerformance.h"
#include "RecoBTag/Records/interface/BTagPerformanceRecord.h"

#include "RecoBTag/SecondaryVertex/interface/TrackKinematics.h"

#include "CondFormats/PhysicsToolsObjects/interface/BinningPointByMap.h"

//JEC uncertainty
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

//test
//#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

//pulled in from Don's code
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "PhysicsTools/SelectorUtils/interface/ElectronVPlusJetsIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/MuonVPlusJetsIDSelectionFunctor.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

//test
//#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"


#include "CUSUSY/BasicTreeMaker/interface/BasicTreeMaker.h"


//
// constructors and destructor
//
BasicTreeMaker::BasicTreeMaker(const edm::ParameterSet& iConfig) :
  thetimer_(new edm::CPUTimer()),
  Heventcount_(0),
  tree_(0),
  infotree_(0),

  pdffactory_("cteq66.LHgrid", 44),

  isMC_(iConfig.getParameter<bool>("MCflag")),
  doPrescale_(true),

  btagAlgorithmNames_(iConfig.getParameter<std::vector<std::string> >("btagAlgorithms")),
  tauidAlgorithmNames_(iConfig.getParameter<std::vector<std::string> >("tauidAlgorithms")),
  triggersOfInterest_(iConfig.getParameter<std::vector<std::string> >("triggersOfInterest")),

  pvSelector_      (iConfig.getParameter<edm::ParameterSet>("pvSelector") ),

  pvLabel_         (iConfig.getParameter<edm::InputTag>("pvTag")),

  jetAlgorithmNames_(iConfig.getParameter<std::vector<std::string> >("jetAlgorithms")),
  jetAlgorithmTags_(iConfig.getParameter<std::vector<std::string> >("jetNames")),
  
  metAlgorithmNames_(iConfig.getParameter<std::vector<std::string> >("metAlgorithms")),
  metAlgorithmTags_(iConfig.getParameter<std::vector<std::string> >("metNames")),

  eleAlgorithmNames_(iConfig.getParameter<std::vector<std::string> >("eleAlgorithms")),
  muonAlgorithmNames_(iConfig.getParameter<std::vector<std::string> >("muonAlgorithms")),
  tauAlgorithmNames_(iConfig.getParameter<std::vector<std::string> >("tauAlgorithms")),

  eleSelected_(iConfig.getParameter<std::vector<std::string> >("eleSelected")),
  muonSelected_(iConfig.getParameter<std::vector<std::string> >("muonSelected")),

  pfCandSrc_       (iConfig.getParameter<edm::InputTag>("PFCandSource")),

  jetIdLoose_      (iConfig.getParameter<edm::ParameterSet>("jetIdLoose") ),
  jetIdTight_      (iConfig.getParameter<edm::ParameterSet>("jetIdTight") ),
  PFjetIdLoose_    (iConfig.getParameter<edm::ParameterSet>("pfjetIdLoose") ),
  PFjetIdTight_    (iConfig.getParameter<edm::ParameterSet>("pfjetIdTight") ),
  muonId_          (iConfig.getParameter<edm::ParameterSet>("muonId") ),
  electronId_      (iConfig.getParameter<edm::ParameterSet>("electronId") ),

  //need something for hlt Config here?
  processName_(""),

  minJets_         (iConfig.getParameter<int>("minJets")),

  jetPtMin_        (iConfig.getParameter<double>("jetPtMin")), 
  jetEtaMax_       (iConfig.getParameter<double>("jetEtaMax")),
  loosejetPtMin_   (iConfig.getParameter<double>("loosejetPtMin")), 
  loosejetEtaMax_  (iConfig.getParameter<double>("loosejetEtaMax")),

  tolerance_(1e-4)
{

  assert(  jetAlgorithmNames_.size() == jetAlgorithmTags_.size());

  assert(  metAlgorithmNames_.size() == metAlgorithmTags_.size());

  //the fill lepton info method is not very flexible!
  assert(  eleAlgorithmNames_.size() == muonAlgorithmNames_.size() );
  assert( eleAlgorithmNames_.size() == eleSelected_.size() );
  assert( muonAlgorithmNames_.size() == muonSelected_.size() );

  edm::Service<TFileService> fs;
  Heventcount_ = fs->make<TH1D>( "Heventcount"  , "events processed", 1,  0, 1 );
  tree_ = new TTree("tree","tree");
  infotree_ = new TTree("infotree","infotree");
  //branch creation done in beginJob instead of here

}


BasicTreeMaker::~BasicTreeMaker()
{
  delete thetimer_;
}


//
// member functions
//

void BasicTreeMaker::fillPVInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<std::vector<reco::Vertex> > vtxs;
  iEvent.getByLabel(pvLabel_, vtxs);

  for ( unsigned int ipv = 0; ipv<vtxs->size() ; ipv++) {
    reco::Vertex const & pv = vtxs->at(ipv);

    pv_isFake.push_back(pv.isFake());
    
    pv_z.push_back( pv.z());
    pv_ndof.push_back(pv.ndof());
    pv_chi2.push_back(pv.chi2());
    pv_rho.push_back( pv.position().Rho());
    
  }
  
}

//
void
BasicTreeMaker::fillTriggerInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  using namespace edm;

  
  findHLTProcessName(iEvent); //fill processName_
  InputTag trigResultsTag("TriggerResults","",processName_); 

  // get HLT trigger information
  Handle<TriggerResults> HLTResults;
  iEvent.getByLabel(trigResultsTag, HLTResults);
  
  //based on code copied from Don
  if (HLTResults.isValid()) {
    
    const edm::TriggerNames & triggerNames = iEvent.triggerNames(*HLTResults);
    unsigned int trigger_size = HLTResults->size();

    //loop over the list of triggers
    for (unsigned int itrig = 0; itrig< triggersOfInterest_.size(); itrig++) {
      
      unsigned int trigger_position = triggerNames.triggerIndex(triggersOfInterest_.at(itrig));
      unsigned int prescale=1;
      if (trigger_position < trigger_size) {
	bool passed=HLTResults->accept(trigger_position);
	passTrigger.push_back( passed );

        if (!isMC_ && doPrescale_) { //in data, get the prescale
	  prescale =  hltConfig_.prescaleValue(iEvent,iSetup,triggersOfInterest_.at(itrig));
	  //  std::cout<< triggersOfInterest_.at(itrig) <<"Prescale = "<<prescale<<std::endl;
	}
	hltPrescale.push_back( prescale);
      }
      else {	//	std::cout << "WARNING -- Trigger position not found" << std::endl;
	passTrigger.push_back(false);
	hltPrescale.push_back( 0);
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
			   || taudau==223 ||taudau==321 ||taudau ==323 ||taudau ==310||taudau==221 ||taudau==20213) found_tau_had=true;
		  else if (taudau == 14 || taudau==16 || taudau==12) { } //do nothing for neutrinos
		  else if (taudau == 22 || taudau==24 ) { } //do nothing for photons (and W?)
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

  //std::cout<<" == fillMCInfo =="<<std::endl;

  // === get the process qscale (pthat) ===
  edm::Handle<GenEventInfoProduct> genEventInfo;
  iEvent.getByLabel("generator", genEventInfo);

  qScale=genEventInfo->qScale();
  ptHat = genEventInfo->hasBinningValues() ? genEventInfo->binningValues()[0] : 0;
  mcWeight = genEventInfo->weight();


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
    else if (abs(TCand.pdgId())== 23 && TCand.status()==3) { //find Z
      //      std::cout<<"Found a status 3 Z"<<std::endl;
      for (unsigned int iz=0; iz<TCand.numberOfDaughters(); ++iz) {
	//	std::cout<<TCand.daughter(iz)->pdgId()<<"\t"<<TCand.daughter(iz)->status()<<std::endl;
	if ( TCand.daughter(iz)->status()==3 ) {
	  int Zdau = 	abs(TCand.daughter(iz)->pdgId());
	  if  ( Zdau >= 11 && Zdau <= 18 ) {
	    ZDecayMode=Zdau;
	    break;
	  }
	}
      }
    } //end Z block
  }

  // === fill flavor history info ===
  edm::Handle<unsigned int> path;
  iEvent.getByLabel("flavorHistoryFilter", path);
  
  flavorHistory=(int) *path;
  //std::cout<<"flavorHistory = "<<flavorHistory<<std::endl;
  
  // === fill PDF info ===
  bool pdfInfook = pdffactory_.extractInfo(iEvent);
  if (!pdfInfook) std::cout << "WARNING -- PDF information not available!" << std::endl;
  // Get the PDF weights and print them::
  std::vector<double>& pdfweights = pdffactory_.weights();
  // std::cout << "PDF weights for event : " << iEvent.id().event() << std::endl;
  for (unsigned int i=0; i<pdfweights.size(); i++) {
    //     cout << "Weight " << i << ":\t" << pdfweights[i] << endl;
    pdfWeights.push_back(pdfweights[i]); //copy into the tree array
  }
  //cout << endl;

  // === get mSugra scan information (if available) ===
  //in reality the handle can probably be reused. but who cares....
  edm::Handle<double> mSugraHandle1;
  iEvent.getByLabel("susyScantanbeta", mSugraHandle1);

  edm::Handle<double> mSugraHandle2;
  iEvent.getByLabel("susyScanA0", mSugraHandle2);

  edm::Handle<double> mSugraHandle3;
  iEvent.getByLabel("susyScanCrossSection", mSugraHandle3);

  edm::Handle<double> mSugraHandle4;
  iEvent.getByLabel("susyScanM0", mSugraHandle4);

  edm::Handle<double> mSugraHandle5;
  iEvent.getByLabel("susyScanM12", mSugraHandle5);

  edm::Handle<double> mSugraHandle6;
  iEvent.getByLabel("susyScanMu", mSugraHandle6);

  edm::Handle<double> mSugraHandle7;
  iEvent.getByLabel("susyScanRun", mSugraHandle7);

  //some of these numbers are probably really integers,
  //so i'm tempted to cast them as such here
  //but i don't want to tempt fate

  if (mSugraHandle1.isValid())  susy_tanBeta = *mSugraHandle1;
  if (mSugraHandle2.isValid())  susy_A0      = *mSugraHandle2;
  if (mSugraHandle3.isValid())  susy_crossSection=*mSugraHandle3;
  if (mSugraHandle4.isValid())  susy_m0      = *mSugraHandle4;
  if (mSugraHandle5.isValid())  susy_m12     = *mSugraHandle5;
  if (mSugraHandle6.isValid())  susy_mu      = *mSugraHandle6;
  if (mSugraHandle7.isValid())  susy_run     = *mSugraHandle7;

}

void
BasicTreeMaker::fillTrackInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<reco::TrackCollection> generalTracks;
  iEvent.getByLabel("generalTracks", generalTracks); //TODO make this name configurable in py

  //also need the main PV
  edm::Handle<std::vector<reco::Vertex> > vtxs;
  iEvent.getByLabel(pvLabel_, vtxs);
  const reco::Vertex * vtx = &((*vtxs)[0]);

  double sumpt=0;
  //no need to bother with sorting tracks
  for (reco::TrackCollection::const_iterator trit=generalTracks->begin(); trit!=generalTracks->end() ; trit++) {
    //for general ntuple storage
    if (trit->pt() > 5 && fabs(trit->eta()) < 5) { //cuts TODO make them configurable
      trackPt.push_back(trit->pt());
      trackEta.push_back(trit->eta());
      trackPhi.push_back(trit->phi());
    }

    //for the tracking failure filter ...not so good that i'm hard-coding this stuff...
    if (fabs(trit->dz(vtx->position())) > 1) continue;
    if (fabs(trit->dxy(vtx->position())) > 0.2) continue;
    sumpt += trit->pt();
  }
  SumPtOverHT = pfht_>0  ? sumpt/pfht_ : -1; //for tracking failure detection. can apply cut offline

}

void
BasicTreeMaker::fillTauInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup, unsigned int il) {
  std::string tauTag=tauAlgorithmNames_[il];

  //  std::cout<<" == "<<tauTag<<" =="<<std::endl;

  edm::Handle<edm::View<pat::Tau> > tauHandle;
  iEvent.getByLabel(tauTag,tauHandle);
  const edm::View<pat::Tau> & taus = *tauHandle;
  for (edm::View<pat::Tau>::const_iterator itau = taus.begin(); itau!=taus.end(); ++itau) {
    //    std::cout<<itau->pt()<<"\t"<<itau->eta()<<"\t"<<itau->phi()<<std::endl;

    tauPt[tauTag].push_back( itau->pt() );
    tauEta[tauTag].push_back( itau->eta());
    tauPhi[tauTag].push_back( itau->phi());
    tauTaNC[tauTag].push_back(itau->tauID("byTaNC")); //hard-code this one because it is the only float (rest are bool)
    for (unsigned int ialg=0; ialg<tauidAlgorithmNames_.size(); ialg++) {
      tauID[tauTag][tauidAlgorithmNames_[ialg]].push_back( itau->tauID(tauidAlgorithmNames_[ialg]) );
      //      std::cout<<tauidAlgorithmNames_[ialg]<<" = "<<itau->tauID(tauidAlgorithmNames_[ialg])<<std::endl;
    }

  }
  
}

  //Code from RA2 for filtering fake MHT with muons:

void BasicTreeMaker::badPFMuonFilter(const edm::Event& iEvent, edm::InputTag pfCandSource, edm::InputTag muonSource, double maxPtDiff, bool doPtDiff, bool doPJCut, bool debug){

  passesBadPFMuonFilter=true;

  //std::cout<<"[badPFMuonFilter] tag: "<<muonSource<<std::endl;

  edm::Handle<reco::PFCandidateCollection> pfCands;
  iEvent.getByLabel(pfCandSource, pfCands);
  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByLabel(muonSource, muons);
  reco::PFCandidateCollection pfCandidates_ = *pfCands;
//  double metBefore = 0.;
//  double metAfter = 0.;
//  if (doPJCut) {
//    bool fixed = postMuonCleaning(&pfCandidates_,metBefore,metAfter);
//    if ( fixed ) return false;
//  }
  if (doPtDiff) {
    for (reco::PFCandidateCollection::const_iterator c = pfCands->begin(); c != pfCands->end(); ++c) {
      if (std::abs(c->pdgId()) != 13) continue;
      if (!muons.failedToGet() && muons.isValid()) { // if these are pat muons
	for (edm::View<pat::Muon>::const_iterator m = muons->begin(); m != muons->end(); ++m) {
	  float dr = reco::deltaR(c->eta(), c->phi(), m->eta(), m->phi());
	  if (m->originalObjectRef().id() == c->muonRef().id() &&
	      m->originalObjectRef().key() == c->muonRef().key()) {
	    if (debug) std::cout << "pf-reco muon match (mindr/pf pt/reco pt): " << dr << " " << c->pt() << " " << m->pt() << std::endl;
	    if (c->pt() - m->pt() > maxPtDiff) passesBadPFMuonFilter=false;
	    if (c->pt() - m->pt() > 25) muonPtDiff_PF.push_back(c->pt() - m->pt());
	  }
	}
      }
      //else { // assume reco::Muon ref from pfcandidate will be present
      //	float dr = reco::deltaR(c->eta(), c->phi(), c->muonRef()->eta(), c->muonRef()->phi());
      //	if (debug) std::cout << "pf-reco muon match (mindr/pf pt/reco pt): " << dr << " " << c->pt() << " " << c->muonRef()->pt() << std::endl;
      //	if (c->pt() - c->muonRef()->pt() > maxPtDiff) return false;
      //}
    }
  }
  // if no bad match found, return true

}

void BasicTreeMaker::inconsistentMuonPFCandidateFilter(const edm::Event& iEvent, edm::InputTag muonSource, double ptMin, double maxPTDiff, bool verbose){
  using namespace std;  
  using namespace edm;

  // cout<<"[inconsistentMuonPFCandidateFilter] tag: "<<muonSource<<endl;


  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByLabel(muonSource, muons);
  
  passesInconsistentMuonPFCandidateFilter=true;

  if (!muons.failedToGet() && muons.isValid()) { // if these are pat muons
    for (edm::View<pat::Muon>::const_iterator muon = muons->begin(); muon != muons->end(); ++muon) {
      if ( muon->pt() < ptMin )                          continue;       
      if (  muon->isTrackerMuon()
	    && muon->isGlobalMuon()
	    && fabs(muon->innerTrack()->pt()/muon->globalTrack()->pt() - 1) <= maxPTDiff   )  continue; 
      passesInconsistentMuonPFCandidateFilter = false;
    }
  }

}

//more or less copied and pasted from
//http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/RecoParticleFlow/PFProducer/src/PFMuonAlgo.cc?revision=1.20&view=markup
//except that I'm changing it to use pat::Muon instead of a 'muonRef'
//This is used only for my half-way implementation of the greedy muon filter
bool
BasicTreeMaker::isIsolatedMuon( const pat::Muon& muonRef ){
  // Isolated Muons which are missed by standard cuts are nearly always global+tracker
  if ( !muonRef.isGlobalMuon() ) return false;

  // If it's not a tracker muon, only take it if there are valid muon hits

  reco::TrackRef standAloneMu = muonRef.standAloneMuon();

  if ( !muonRef.isTrackerMuon() ){
    if(standAloneMu->hitPattern().numberOfValidMuonDTHits() == 0 &&
       standAloneMu->hitPattern().numberOfValidMuonCSCHits() ==0) return false;
  }
  
  // for isolation, take the smallest pt available to reject fakes

  reco::TrackRef combinedMu = muonRef.combinedMuon();
  double smallestMuPt = combinedMu->pt();
  
  if(standAloneMu->pt()<smallestMuPt) smallestMuPt = standAloneMu->pt();
  
  if(muonRef.isTrackerMuon())
    {
      reco::TrackRef trackerMu = muonRef.track();
      if(trackerMu->pt() < smallestMuPt) smallestMuPt= trackerMu->pt();
    }
     
  double sumPtR03 = muonRef.isolationR03().sumPt;
  double emEtR03 = muonRef.isolationR03().emEt;
  double hadEtR03 = muonRef.isolationR03().hadEt;
  
  double relIso = (sumPtR03 + emEtR03 + hadEtR03)/smallestMuPt;
  //double relIso = sumPtR03/muonRef.pt();

  /*
  // protection against fake tracks -- only needed for tracker only, which isn't found in signal sample
  if(muonRef.isTrackerMuon() && !muonRef.isGlobalMuon() && !muonRef.isStandAloneMuon())
    {
      reco::TrackRef trackerMu = muonRef.track();
      
      if(trackerMu->ptError()/trackerMu->pt() > 0.20) return false;
    }
  //std::cout<<" relIso "<<relIso<<std::endl;
  */

  if (relIso<0.1) return true;
  else return false;
}


void
BasicTreeMaker::fillLeptonInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup, unsigned int il) {

  const bool eleDebug=false;
  const bool muonDebug=false;

  std::string muTag=muonAlgorithmNames_[il];
  std::string eTag=eleAlgorithmNames_[il];
  
  if (muonDebug || eleDebug)  std::cout<<muTag<<"\t"<<eTag <<std::endl;
  
  //Added by Luke -- Trying to get functions for the muons that are implemented by RA2 to get rid of fake MHT
  //note -- these used to be called inside the muon loop, but that doesn't seem right at all
  //that means they get called 0 times if there are 0 muons, but then the bool never got filled.
  //in short, it was wrong. now call them exactly once per event!

  //for  badPFMuonFilter we want to call it only for the non-PF muon collection
  if (muTag.find("PF")==std::string::npos) badPFMuonFilter(iEvent, pfCandSrc_,edm::InputTag(muTag)) ;
  //for  the inconsistent muon we want to call it for the PF muon collection
  if (muTag.find("PF")!=std::string::npos) inconsistentMuonPFCandidateFilter(iEvent, edm::InputTag(muTag)) ;

  //to get the cut flow in the right (arbitrary) order, need to do muons first
  edm::Handle<edm::View<pat::Muon> > muonHandle;
  iEvent.getByLabel(muTag,muonHandle);
  const edm::View<pat::Muon> & muons = *muonHandle;

  //get the RA2 muon list
  edm::Handle<edm::View<pat::Muon> > muonHandleRA2;
  iEvent.getByLabel(muonSelected_[il],muonHandleRA2);
  const edm::View<pat::Muon> & muonsRA2 = *muonHandleRA2;

  for (edm::View<pat::Muon>::const_iterator imuon = muons.begin(); imuon!=muons.end(); ++imuon) {

    //check if this muon is also on the RA2 list
    bool foundOverlap=false;
    for (edm::View<pat::Muon>::const_iterator iRA2muon = muonsRA2.begin(); iRA2muon!=muonsRA2.end(); ++iRA2muon) {
      //ok, for pat::Muon this method worked. for non-PF pat::Electron is seems to work
      //but for PF pat::Electron it did not work.
      //this makes me nervous, so I go with a simple pT,eta,phi match
//       if ( ((reco::RecoCandidate*) &(*imuon))-> overlap(  *iRA2muon ) ) {
// 	foundOverlap=true;
// 	//	std::cout<<"matching muons with pT = "<<imuon->pt()<<" "<<iRA2muon->pt()<<std::endl;
//       }
      if ( ( fabs(imuon->pt() - iRA2muon->pt()) < tolerance_) 
	   && (fabs(imuon->eta() - iRA2muon->eta()) < tolerance_) 
	   && (fabs(imuon->phi() - iRA2muon->phi()) < tolerance_)) {
 	foundOverlap=true;
	//	std::cout<<"matching muons with pT = "<<imuon->pt()<<" "<<iRA2muon->pt()<<std::endl;
      }
      //      else 	std::cout<<"NOT matching muons with pT = "<<imuon->pt()<<" "<<iRA2muon->pt()<<std::endl;
    }
    muonIsRA2[muTag].push_back(foundOverlap);

    //this is a pseudo-implementation of the greedy muon filter.
    //it does *not* actually faithfully reimplement it; that would have required RECO (I think)
    //instead of looking at all PFCandidates that are muons, I start with the muon list
    if (imuon->pfCandidateRef().isNonnull() ) {
      bool isIso = isIsolatedMuon(*imuon);
    
      if (isIso) {
	double totalCaloEnergy = imuon->pfCandidateRef()->rawEcalEnergy() +  imuon->pfCandidateRef()->rawHcalEnergy();
	float eOverP = totalCaloEnergy/imuon->pfCandidateRef()->p(); 
	muonEoverP[muTag].push_back( eOverP );
      }
      else   muonEoverP[muTag].push_back( -1 );
    }
    else   muonEoverP[muTag].push_back( -2 );

    //now storing ALL muons
    const  bool isGlobal = imuon->muonID("AllGlobalMuons");
    muonIsGlobalMuon[muTag].push_back(isGlobal);

    muonIsGlobalMuonPromptTight[muTag].push_back(imuon->muonID("GlobalMuonPromptTight"));

    //record pT and eta of all that pass 
    muonPt[muTag].push_back( imuon->pt() );
    muonEta[muTag].push_back( imuon->eta());
    muonPhi[muTag].push_back( imuon->phi());

//     if (muonDebug) {
//       std::cout<<"muon pT: "<<imuon->pt()<<" "<<imuon->innerTrack()->pt()<<std::endl;
//       std::cout<<"muon isolation pairs: "<< //these are confirmed match in 384 running over a data PATtuple
// 	imuon->ecalIso()<<" "<<imuon->isolationR03().emEt<<"\t"<<
// 	imuon->hcalIso()<<" "<<imuon->isolationR03().hadEt<<"\t"<<
// 	imuon->trackIso()<<" "<<imuon->isolationR03().sumPt<<std::endl;
//     }
    muonTrackIso[muTag].push_back( imuon->trackIso() );
    muonEcalIso[muTag].push_back( imuon->ecalIso() );
    muonHcalIso[muTag].push_back( imuon->hcalIso() );

    //do this stuff only for global muons
    if (isGlobal) {
      if (!imuon->combinedMuon().isNull()) {
	muonChi2[muTag].push_back(imuon->combinedMuon()->chi2());
	muonNdof[muTag].push_back(imuon->combinedMuon()->ndof());
      }
      else {
	muonChi2[muTag].push_back(0);
	muonNdof[muTag].push_back(0);
      }
      
      muonTrackd0[muTag].push_back(imuon->track()->d0() );
      muonTrackPhi[muTag].push_back(imuon->track()->phi() );
      
      muonNhits[muTag].push_back(imuon->numberOfValidHits());
      float ptratio = imuon->globalTrack()->pt() == 0 ? -99 : imuon->innerTrack()->pt()/imuon->globalTrack()->pt();
      muonPtRatio[muTag].push_back(ptratio);
    }
    else { //non-global
      muonChi2[muTag].push_back(0);
      muonNdof[muTag].push_back(0);
      muonTrackd0[muTag].push_back(-99 );
      muonTrackPhi[muTag].push_back(-99 );
      muonNhits[muTag].push_back(-99);
      muonPtRatio[muTag].push_back(-99);
    }

    // Muon veto
    bool passMuon = isGlobal ? muonId_(*imuon,iEvent) : false; //this includes isolation cut (and other things)
    muonPassID[muTag].push_back(passMuon);
    
    //store vtx z
    muonVtx_z[muTag].push_back(  imuon->vertex().z());
  } //end of loop over muons

  // it is really stupid that I am using a different style of logic for muons and electrons
  //But I can't decide which I like better

  if (eleDebug)   std::cout<<"--elec--"<<std::endl;

  edm::Handle<edm::View<pat::Electron> > electronHandle;
  iEvent.getByLabel(eTag,electronHandle);
  if (!electronHandle.isValid()) {std::cout<<"electrons are not valid!"<<std::endl;}
  const edm::View<pat::Electron> & electrons = *electronHandle;

  //get the RA2 electron list
  edm::Handle<edm::View<pat::Electron> > eleHandleRA2;
  iEvent.getByLabel(eleSelected_[il],eleHandleRA2);
  const edm::View<pat::Electron> & electronsRA2 = *eleHandleRA2;

  if (eleDebug) {
    std::cout<<" using RA2 electron list: "<<eleSelected_[il]<<std::endl
	     <<"N electrons, N RA2 electrons = "<<electrons.size()<<", "<<electronsRA2.size()<<std::endl;
  }


  for (edm::View<pat::Electron>::const_iterator ielectron = electrons.begin(); ielectron!=electrons.end(); ++ielectron) {

    //check if this muon is also on the RA2 list
    bool foundOverlap=false;
    for (edm::View<pat::Electron>::const_iterator iRA2ele = electronsRA2.begin(); iRA2ele!=electronsRA2.end(); ++iRA2ele) {

      //this method almost worked, but failed for the PF electrons
//       if ( ((reco::RecoCandidate*) &(*ielectron))-> overlap(  *iRA2ele ) ) {
// 	foundOverlap=true;
// 	std::cout<<"matching electrons with pT = "<<ielectron->pt()<<" "<<iRA2ele->pt()<<std::endl;
// 	//	std::cout<<"matching electrons with eT = "<<ielectron->et()<<" "<<iRA2ele->et()<<std::endl;
//       }
//       else {
// 	std::cout<<"NOT matching electrons with pT = "<<ielectron->pt()<<" "<<iRA2ele->pt()<<std::endl;
// 	//	std::cout<<"matching electrons with eT = "<<ielectron->et()<<" "<<iRA2ele->et()<<std::endl;
// 	}
//       }

      if ( (fabs(ielectron->pt() - iRA2ele->pt()) < tolerance_) 
	   && (fabs(ielectron->eta() - iRA2ele->eta()) < tolerance_) 
	   && (fabs(ielectron->phi() - iRA2ele->phi()) < tolerance_)) {
	foundOverlap=true;
	//	std::cout<<"matching electrons with pT = "<<ielectron->pt()<<" "<<iRA2ele->pt()<<std::endl;
      }
      //      else 	std::cout<<"NOT matching electrons with pT = "<<ielectron->pt()<<" "<<iRA2ele->pt()<<std::endl;
    }

    eleIsRA2[eTag].push_back(foundOverlap);

    eleIDLoose[eTag].push_back(ielectron->electronID( "eidLoose" ) );
    eleIDRobustTight[eTag].push_back(ielectron->electronID( "eidRobustTight" ) );

    eleEt[eTag].push_back( ielectron->et() );
    //    elePt[eTag].push_back( ielectron->pt() );
    eleEta[eTag].push_back( ielectron->eta());
    elePhi[eTag].push_back( ielectron->phi());
    eleTrackIso[eTag].push_back( ielectron->dr03TkSumPt() );
    eleEcalIso[eTag].push_back( ielectron->dr03EcalRecHitSumEt() );
    eleHcalIso[eTag].push_back( ielectron->dr03HcalTowerSumEt() );

    eledB[eTag].push_back( ielectron->dB());

    eleVtx_z[eTag].push_back( ielectron->vertex().z() );

    elePassID[eTag].push_back(electronId_(*ielectron) );
    
  } //end of loop over electrons
  if (eleDebug)   std::cout<<"--leptons done--"<<std::endl;

  //  leptonInfoFilled_=true;
}


//this function is not actually used! stick with the JetIDSelectionFunctor
//...and it would probably give the wrong output, since it operates on the corrected jets
bool
BasicTreeMaker::passJetId(const pat::Jet & jet) {

  //  if (jet.isCaloJet() ) std::cout<<"This is a calo jet"<<std::endl;
  //  if (jet.isJPTJet() ) std::cout<<"This is a JPT jet"<<std::endl;
  if ( !jet.isPFJet() ) {//std::cout<<"This is a PF jet"<<std::endl;
    
    reco::JetID myjetid = jet.jetID();
        
    std::cout<<"isJPT="<<jet.isJPTJet()<<"\t"<<jet.jetID().n90Hits<<"\t"<<jet.jetID().hitsInN90<<std::endl;
    float jet_n90Hits = jet.jetID().n90Hits;
    float jet_fHPD = jet.jetID().fHPD;
    float jet_emf = jet.emEnergyFraction();
    
    //std::cout<<jet_n90Hits<<"\t"<<jet_fHPD<<"\t"<<jet_emf<<std::endl;
    
    bool pass_n90Hits = jet_n90Hits > 1;
    bool pass_fHPD = jet_fHPD < 0.98;
    bool pass_emf = true;
    if ( fabs(jet.eta()) <2.6 ) {
      if (jet_emf <= 0.01) pass_emf=false;
    }
    else {                // HF
      if( jet_emf <= -0.9 ) pass_emf = false;
      if( jet.pt() > 80 && jet_emf >= 1 ) pass_emf = false;
    }
        
    return (pass_n90Hits && pass_fHPD && pass_emf);
  }

  return true;

}

void
BasicTreeMaker::fillJetInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup, unsigned int jetIndex)
{
  const bool jetDebug=false;
  if (jetDebug)  std::cout<<" == fillJetInfo "<<jetAlgorithmNames_[jetIndex]<<" =="<<std::endl; //debug
 
  //  if (jetInfoFilled_) return;

  JetCorrectionUncertainty *jecUncPF =0; //to be used inside the jet loop

  edm::Handle<edm::View<pat::Jet> > jetHandle;
  iEvent.getByLabel(jetAlgorithmNames_[jetIndex],jetHandle);
  const edm::View<pat::Jet> & jets = *jetHandle;

  //freya sort using et. don sorts with pt. stick with pt
  //further note on this -- RA2 selection sorts by pT
  std::vector<float> jetPtVector;
  for(edm::View<pat::Jet>::const_iterator ijet = jets.begin(); ijet != jets.end(); ++ijet) {
    jetPtVector.push_back( ijet->pt() );
  }

  std::vector<int> sortedJetIndices = IndexSorter< std::vector<float> >(jetPtVector,true)();

  //variables for use in the jet loop
  pat::strbitset ret1 = jetIdLoose_.getBitTemplate();
  pat::strbitset retpf = PFjetIdLoose_.getBitTemplate();
  pat::strbitset ret1t = jetIdTight_.getBitTemplate();
  pat::strbitset retpft = PFjetIdTight_.getBitTemplate();
  //now loop over sorted jets
  for (size_t ii = 0 ; ii< sortedJetIndices.size(); ++ii) {
    //get the index of the ii'th jet
    int jj = sortedJetIndices[ii];
    //now pull out that jet
    const pat::Jet & jet = jets[jj];
    
    //note that this code will fail if we ever run over 2 different PF jet collections in parallel!
    if (jet.isPFJet()) pfht_+= jet.pt(); //for tracking failure filter, compute PF HT with no cuts

    /*
      - canned code for getting the JEC uncertainties from the DB
      - note that it is hard-coded for PF -- this is kludge that I should eventually fix
      - it is unclear to me if this can be done less often than once per event. so stick to this for now
    */
    if ( jecUncPF == 0 && jet.isPFJet() ) { // do this only once per event and only for PF jets
      // handle the jet corrector parameters collection
      edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
      // get the jet corrector parameters collection from the global tag
      iSetup.get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl);
      // get the uncertainty parameters from the collection
      JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
      // instantiate the jec uncertainty object
      jecUncPF = new JetCorrectionUncertainty(JetCorPar);
    }

    //std::cout<<"about to get jet id"<<std::endl;
    
    bool passJetID = false, passTightJetID=false;
    if ( jet.isPFJet() ) {
      passJetID = PFjetIdLoose_(jet,retpf);
      //      try { 
      passTightJetID = PFjetIdTight_(jet,retpft); //} catch (...) {std::cout<<"Problem in PFjetIdTight!"<<std::endl;}
    }
    else if (jet.isJPTJet() || jet.isCaloJet() ) {
      passJetID = jetIdLoose_(jet, ret1);
      //      try { 
      passTightJetID = jetIdTight_(jet, ret1t); //} catch (...) {std::cout<<"Problem in jetIdTight!"<<std::endl;}
    }
    else {std::cout<<"Unknown jet type!"<<std::endl;}

    // === make cuts ===
    bool passLooseCuts = (jet.pt() > loosejetPtMin_) && (fabs(jet.eta())<loosejetEtaMax_);
    bool passTightCuts = (jet.pt() > jetPtMin_) && (fabs(jet.eta())<jetEtaMax_) && passJetID;

    if (jetDebug && jet.isPFJet()) {
      std::cout<<"jet "<<ii<< " JEC="<<jet.currentJECLevel()
	       <<" "<<jet.pt()<<" "<<jet.phi()<<" "<<jet.eta()<<" "<<passJetID<<std::endl;
      pat::Jet ujet = jet.correctedJet("Uncorrected");
      //       std::cout<<"jet "<<ii<< " JEC="<<ujet.currentJECLevel()
      // 	       <<" "<<ujet.pt()<<" "<<ujet.phi()<<" "<<ujet.eta()<<std::endl;
      std::cout<<"Uncorrected: "<< (ujet.neutralHadronEnergy() + ujet.HFHadronEnergy() ) / ujet.energy() 
	       << " "<<ujet.neutralEmEnergyFraction() 
	       << " "<<ujet.numberOfDaughters() 
	       << " "<<ujet.chargedHadronEnergyFraction()
	       << " "<<ujet.chargedMultiplicity() 
	       << " "<<ujet.chargedEmEnergyFraction()<<std::endl;
      std::cout<<"  corrected: "<< (jet.neutralHadronEnergy() + jet.HFHadronEnergy() ) / jet.energy() 
	       << " "<<jet.neutralEmEnergyFraction() 
	       << " "<<jet.numberOfDaughters() 
	       << " "<<jet.chargedHadronEnergyFraction()
	       << " "<<jet.chargedMultiplicity() 
	       << " "<<jet.chargedEmEnergyFraction()<<std::endl;

    }
    
    //let's enforce that the tight cuts are always a subset of the loose cuts
    if (passTightCuts && !passLooseCuts) assert(0);
    
    //now fill ntuple!
    /* some code used just printing some info about jet corrections during development
       std::vector<std::string> jecs = jet.availableJECLevels();
       for (unsigned int ijecs = 0; ijecs< jecs.size(); ijecs++)     std::cout<<jecs.at(ijecs)<<" ";
       std::cout<<std::endl;
       
       std::cout<<"using jeclevel = "<<jet.currentJECLevel()<<std::endl;
    */

    if (passLooseCuts) {
      //      looseJetIndex[jetAlgorithmTags_[jetIndex]].push_back( veryloosejetPtUncorr.size() - 1);
      loosejetPt[jetAlgorithmTags_[jetIndex]].push_back( jet.pt() );
      loosejetPz[jetAlgorithmTags_[jetIndex]].push_back( jet.pz() );
      loosejetPtUncorr[jetAlgorithmTags_[jetIndex]].push_back( jet.correctedJet("Uncorrected").pt() );
      loosejetEt[jetAlgorithmTags_[jetIndex]].push_back( jet.et() );
      loosejetE[jetAlgorithmTags_[jetIndex]].push_back( jet.energy() );
      loosejetEta[jetAlgorithmTags_[jetIndex]].push_back(jet.eta());
      loosejetPhi[jetAlgorithmTags_[jetIndex]].push_back(jet.phi());
      loosejetFlavor[jetAlgorithmTags_[jetIndex]].push_back( jet.partonFlavour() );

      //JEC uncertainties
      //for now we are hard-coding to only store them for PF jets
      if (jet.isPFJet() ) {
	if (fabs(jet.eta())<5) {
	  jecUncPF->setJetEta(jet.eta());
	  jecUncPF->setJetPt(jet.pt()); //corrected pT
	  const 	float uncplus = jecUncPF->getUncertainty(true);
	  loosejetJECUncPlus[jetAlgorithmTags_[jetIndex]].push_back( uncplus );
	  //due to weird behavior of this JetUncertainty class, need to reset the eta and pt
	  jecUncPF->setJetEta(jet.eta());
	  jecUncPF->setJetPt(jet.pt()); //corrected pT
	  const 	float uncminus = jecUncPF->getUncertainty(false);
	  loosejetJECUncMinus[jetAlgorithmTags_[jetIndex]].push_back( uncminus );
	  if (jetDebug)	std::cout<<"JEC uncertainty = "<<uncplus<<" "<<uncminus<<std::endl;
	}
	//code threw an exception for eta>6 jet
	else {
	  loosejetJECUncPlus[jetAlgorithmTags_[jetIndex]].push_back( 0 );
	  loosejetJECUncMinus[jetAlgorithmTags_[jetIndex]].push_back( 0 );
	}
      }

      //jet id
      loosejetPassLooseID[jetAlgorithmTags_[jetIndex]].push_back( passJetID);
      loosejetPassTightID[jetAlgorithmTags_[jetIndex]].push_back( passTightJetID);
      float hfrac=0;
      if ( jet.isCaloJet() || jet.isJPTJet()) hfrac= jet.energyFractionHadronic();
      else if (jet.isPFJet()) {
	hfrac= jet.chargedHadronEnergyFraction()+jet.neutralHadronEnergyFraction();
      }
      else {std::cout<<"Unknown jet type!"<<std::endl;}
      loosejetEnergyFracHadronic[jetAlgorithmTags_[jetIndex]].push_back(hfrac );

      if ( jet.genJet() != 0 ) {
	loosejetGenPt[jetAlgorithmTags_[jetIndex]].push_back( jet.genJet()->pt() );
	loosejetGenPhi[jetAlgorithmTags_[jetIndex]].push_back( jet.genJet()->phi() );
	loosejetGenEta[jetAlgorithmTags_[jetIndex]].push_back( jet.genJet()->eta() );

	loosejetInvisibleEnergy[jetAlgorithmTags_[jetIndex]].push_back( jet.genJet()->invisibleEnergy());
	
	if ( jet.genParticle() != 0) {
	  loosejetGenParticlePDGId[jetAlgorithmTags_[jetIndex]].push_back( jet.genParticle()->pdgId() );
	  
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

	//note! in this case it is important to fill all variables with dummy values!
	else { //no genparticle match
	  loosejetGenParticlePDGId[jetAlgorithmTags_[jetIndex]].push_back( -99 );
	}
      }
      else { //no genjet match
	loosejetInvisibleEnergy[jetAlgorithmTags_[jetIndex]].push_back( -99 );
	loosejetGenPt[jetAlgorithmTags_[jetIndex]].push_back( -99 );
	loosejetGenPhi[jetAlgorithmTags_[jetIndex]].push_back( -99 );
	loosejetGenEta[jetAlgorithmTags_[jetIndex]].push_back( -99 );
	loosejetGenParticlePDGId[jetAlgorithmTags_[jetIndex]].push_back( -99 );
      }
      
      for (unsigned int ib=0; ib<btagAlgorithmNames_.size(); ib++) {
	loosejetBTagDisc[jetAlgorithmTags_[jetIndex]][btagAlgorithmNames_[ib]].push_back( jet.bDiscriminator(btagAlgorithmNames_[ib]) );
      }


      //Beginnning Secondary Vertex Stuff
      if (jet.tagInfoSecondaryVertex() != 0) {
      int nSV = jet.tagInfoSecondaryVertex()->nVertices();
      loosejetNTracks[jetAlgorithmTags_[jetIndex]].push_back(jet.associatedTracks().size());
      loosejetNSV[jetAlgorithmTags_[jetIndex]].push_back(nSV);
      if(nSV>0)
	{
	  reco::TrackKinematics secondaryVertex(jet.tagInfoSecondaryVertex()->secondaryVertex(0) );
	  double vertexMass = secondaryVertex.vectorSum().M();
	  double weightedVertexMass = secondaryVertex.weightedVectorSum().M();
	  loosejetSVUnWeightedMass[jetAlgorithmTags_[jetIndex]].push_back(vertexMass);
	  loosejetSVWeightedMass[jetAlgorithmTags_[jetIndex]].push_back(weightedVertexMass);
	  double vertexMomentum = secondaryVertex.vectorSum().P();
	  double weightedVertexMomentum = secondaryVertex.weightedVectorSum().P();
	  double vertexDistance = sqrt((jet.tagInfoSecondaryVertex()->secondaryVertex(0).position() - jet.vertex()).mag2());
	  loosejetSVUnWeightedLifetime[jetAlgorithmTags_[jetIndex]].push_back(vertexMass*vertexDistance/vertexMomentum);
	  loosejetSVWeightedLifetime[jetAlgorithmTags_[jetIndex]].push_back(weightedVertexMass*vertexDistance/weightedVertexMomentum);
	  //Now I want the angle of the vertex with the boost direction of the jet.
	  //To boost into the jet's rest frame (where the jet's p4 is (E,px,py,pz) ) we can use \beta=(-px/E,-py/E,-pz/E)
	  //This transforms the jet's p4 into (m1,0,0,0)
	  //applying this transform to the vertex's p4 (E2,px2,py2,pz2) gives an awful mess
	  //However, what we want is the angle of the vertex location w.r.t. the boost axis, which is much nicer.  If we take the dot product of
	  //the momentum 3 vector of (E2,px2,py2,pz2) boosted into the rest frame of the jet with the boost vector \beta and then normalize so that 
	  //we are left with the cosine of the angle between them, the formula is (after asking Mathematica for some simplification):
	  // 
	  //         E*(\vec(p1).\vec(p2))-E2(\vec(p1)^2)
	  //       -----------------------------------------
	  //       sqrt(\vec(p1)^2*sqrt((p1.p2)^2-m1^2*m2^2)
      
	  double cosTheta = ( jet.p4().Vect().mag2() !=0 && pow(secondaryVertex.vectorSum().Dot(jet.p4()),2) != jet.p4().mass2()*secondaryVertex.vectorSum().mass2() ) ? 
	    (jet.p4().energy()*(secondaryVertex.vectorSum().Vect().Dot(jet.p4().Vect()))-secondaryVertex.vectorSum().energy()*(jet.p4().Vect().mag2()))/sqrt(jet.p4().Vect().mag2()*(pow(secondaryVertex.vectorSum().Dot(jet.p4()),2)-jet.p4().mass2()*secondaryVertex.vectorSum().mass2())) : 0 ;
	  
	  double weightedCosTheta = ( jet.p4().Vect().mag2() !=0 && pow(secondaryVertex.weightedVectorSum().Dot(jet.p4()),2) != jet.p4().mass2()*secondaryVertex.weightedVectorSum().mass2() ) ? 
	    (jet.p4().energy()*(secondaryVertex.weightedVectorSum().Vect().Dot(jet.p4().Vect()))-secondaryVertex.weightedVectorSum().energy()*(jet.p4().Vect().mag2()))/sqrt(jet.p4().Vect().mag2()*(pow(secondaryVertex.weightedVectorSum().Dot(jet.p4()),2)-jet.p4().mass2()*secondaryVertex.weightedVectorSum().mass2())) : 0 ;
	     
	  loosejetSVUnWeightedCosTheta[jetAlgorithmTags_[jetIndex]].push_back(cosTheta);
	  loosejetSVWeightedCosTheta[jetAlgorithmTags_[jetIndex]].push_back(weightedCosTheta);
	}
      else
	{
	  loosejetSVUnWeightedMass[jetAlgorithmTags_[jetIndex]].push_back(-100.);
	  loosejetSVWeightedMass[jetAlgorithmTags_[jetIndex]].push_back(-100.);
	  loosejetSVUnWeightedLifetime[jetAlgorithmTags_[jetIndex]].push_back(-100.);
	  loosejetSVWeightedLifetime[jetAlgorithmTags_[jetIndex]].push_back(-100.);
	  loosejetSVUnWeightedCosTheta[jetAlgorithmTags_[jetIndex]].push_back(-100.);
	  loosejetSVWeightedCosTheta[jetAlgorithmTags_[jetIndex]].push_back(-100.);
	}
      }
//       else {
// 	std::cout<<"Found no secondary vertex info!"<<std::endl;
//       }
      //Ending Secondary Vertex Stuff

    } //end of block for loose jets
    if (passTightCuts) {
      //allow us to reference things only stored for loose jets
      tightJetIndex[jetAlgorithmTags_[jetIndex]].push_back( loosejetPt[jetAlgorithmTags_[jetIndex]].size() - 1);
    } //end of tight cuts block
  } //end of loop over jets


  //this object is not always created, but it should be safe to delete a 0 pointer
  delete  jecUncPF;

}
		     

void
BasicTreeMaker::fillMetInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup) {


  for (unsigned int imettype = 0; imettype < metAlgorithmNames_.size(); imettype++) {
    edm::Handle<edm::View<pat::MET> > metHandle; 
    iEvent.getByLabel( metAlgorithmNames_[imettype], metHandle);
    MET[ metAlgorithmTags_[imettype] ] = metHandle->front().et() ;
    METphi[ metAlgorithmTags_[imettype] ] = metHandle->front().phi() ;
    METsig [ metAlgorithmTags_[imettype]]= metHandle->front().mEtSig();

    //      std::cout<<metAlgorithmNames_[imettype] << " "<<MET[ metAlgorithmTags_[imettype] ]<<std::endl;

    if (isMC_ && metHandle->front().genMET()) {
      GenMET[ metAlgorithmTags_[imettype] ] = metHandle->front().genMET()->et() ;
      GenMETphi[ metAlgorithmTags_[imettype] ] = metHandle->front().genMET()->phi() ;
    }
  }
    

}

//this is more or less working.
//next step is to make it fancier.
//first need to decide if we actually care about this info!
void
BasicTreeMaker::fillBTagDBInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //  if (!jetInfoFilled_) assert(0);

  //stealing generously from Freya et al.
  //for now just hard-code stuff

  std::vector<std::string> btaggingparamnames_;
  btaggingparamnames_.push_back("MCCaloSSVHEMb");
  btaggingparamnames_.push_back("MCCaloSSVHEMb");//yes, the duplication is intentional

  std::vector<std::string> btaggingparaminputtypes_;
  btaggingparaminputtypes_.push_back("BTAGBEFF");
  btaggingparaminputtypes_.push_back("BTAGBERR");

  std::map<std::string,float> bidParamsDiscCut_;

  //this is stuff that we would move out of this function when everything is finalized
  std::map<std::string,PerformanceResult::ResultType> btaggingparamtype_;
  btaggingparamtype_["BTAGBEFF"]=PerformanceResult::BTAGBEFF;
  btaggingparamtype_["BTAGBERR"]=PerformanceResult::BTAGBERR;
  btaggingparamtype_["BTAGCEFF"]=PerformanceResult::BTAGCEFF;
  btaggingparamtype_["BTAGCERR"]=PerformanceResult::BTAGCERR;
  btaggingparamtype_["BTAGLEFF"]=PerformanceResult::BTAGLEFF;
  btaggingparamtype_["BTAGLERR"]=PerformanceResult::BTAGLERR;
  btaggingparamtype_["BTAGNBEFF"]=PerformanceResult::BTAGNBEFF;
  btaggingparamtype_["BTAGNBERR"]=PerformanceResult::BTAGNBERR;
  btaggingparamtype_["BTAGBEFFCORR"]=PerformanceResult::BTAGBEFFCORR;
  btaggingparamtype_["BTAGBERRCORR"]=PerformanceResult::BTAGBERRCORR;
  btaggingparamtype_["BTAGCEFFCORR"]=PerformanceResult::BTAGCEFFCORR;
  btaggingparamtype_["BTAGCERRCORR"]=PerformanceResult::BTAGCERRCORR;
  btaggingparamtype_["BTAGLEFFCORR"]=PerformanceResult::BTAGLEFFCORR;
  btaggingparamtype_["BTAGLERRCORR"]=PerformanceResult::BTAGLERRCORR;
  btaggingparamtype_["BTAGNBEFFCORR"]=PerformanceResult::BTAGNBEFFCORR;
  btaggingparamtype_["BTAGNBERRCORR"]=PerformanceResult::BTAGNBERRCORR;
  btaggingparamtype_["BTAGNBERRCORR"]=PerformanceResult::BTAGNBERRCORR;
  btaggingparamtype_["MUEFF"]=PerformanceResult::MUEFF;
  btaggingparamtype_["MUERR"]=PerformanceResult::MUERR;
  btaggingparamtype_["MUFAKE"]=PerformanceResult::MUFAKE; 
  btaggingparamtype_["MUEFAKE"]=PerformanceResult::MUEFAKE;

  /*
FIXME for now i'm going to hardcode the "calo" here
  */

  edm::ESHandle<BtagPerformance> perfH;
  BinningPointByMap p;
  for (size_t ii=0; ii<btaggingparamnames_.size(); ii++){
     std::string beffstr = btaggingparamnames_[ii];
     //std::cout <<" Studying Beff with label "<<beffstr <<std::endl; //debug
     iSetup.get<BTagPerformanceRecord>().get(beffstr,perfH);  
     const BtagPerformance & pbeff = *(perfH.product());
     bidParamsDiscCut_[btaggingparamnames_[ii]]=pbeff.workingPoint().cut();
     //std::cout << " Beff cut value is : " << bidParamsDiscCut_[btaggingparamnames_[ii]] << std::endl; //debug
     // now loop over the jets:
     std::string lookupname=btaggingparamnames_[ii]+"-"+btaggingparaminputtypes_[ii];
     for (size_t ijet=0; ijet< loosejetEt["calo"].size(); ijet++) {
       p.reset();
       p.insert(BinningVariables::JetAbsEta,std::abs( loosejetEta["calo"][ijet] ));
       p.insert(BinningVariables::JetEt,loosejetEt["calo"][ijet]);
       float eff=0;
       //	 if(pbeff.isResultOk(btaggingparamtype_[jj].second,p))
       eff=pbeff.getResult(btaggingparamtype_[btaggingparaminputtypes_[ii]],p);
       std::cout << "value " << lookupname << " (PerformanceResult:" << btaggingparamtype_[btaggingparaminputtypes_[ii]] <<  ") for jet " << ijet << "(et,eta):("<<loosejetEt["calo"][ijet]<< "," << loosejetEta["calo"][ijet]<< ")  is " << eff << std::endl;
       //       jetSortedBIDParams_[lookupname + ID][ijet]=eff;

     }
  }

}

void
BasicTreeMaker::resetTreeVariables() {
  pfht_=0;

  //reset tree variables
  passTrigger.clear();
  hltPrescale.clear();

  bsx=bsy=bsz=0;

  passesBadPFMuonFilter=passesInconsistentMuonPFCandidateFilter=true;

  for (std::vector<std::string>::const_iterator ij=jetAlgorithmTags_.begin() ; ij!=jetAlgorithmTags_.end() ; ++ij) {
    tightJetIndex[*ij].clear();
    looseJetIndex[*ij].clear();
    //   jetPt.clear();
    //   jetEta.clear();
    //   jetPhi.clear();
    //   jetFlavor.clear();
    loosejetPt[*ij].clear();
    loosejetPz[*ij].clear();
    loosejetPtUncorr[*ij].clear();
    loosejetEt[*ij].clear();
    loosejetE[*ij].clear();
    loosejetEta[*ij].clear();
    loosejetPhi[*ij].clear();

    loosejetJECUncPlus[*ij].clear();
    loosejetJECUncMinus[*ij].clear();

    loosejetGenPt[*ij].clear();
    loosejetGenEta[*ij].clear();
    loosejetGenPhi[*ij].clear();

    loosejetPassLooseID[*ij].clear();
    loosejetPassTightID[*ij].clear();
    loosejetEnergyFracHadronic[*ij].clear();
    loosejetFlavor[*ij].clear();
    loosejetGenParticlePDGId[*ij].clear();
    loosejetInvisibleEnergy[*ij].clear();
    loosejetNTracks[*ij].clear();
    loosejetNSV[*ij].clear();
    loosejetSVWeightedMass[*ij].clear();
    loosejetSVUnWeightedMass[*ij].clear();
    loosejetSVWeightedLifetime[*ij].clear();
    loosejetSVUnWeightedLifetime[*ij].clear();
    loosejetSVWeightedCosTheta[*ij].clear();
    loosejetSVUnWeightedCosTheta[*ij].clear();

    for (unsigned int ib=0; ib<btagAlgorithmNames_.size(); ib++) {
      loosejetBTagDisc[*ij][btagAlgorithmNames_[ib]].clear();
    }
  }

  muonPtDiff_PF.clear();
  for (unsigned int il=0; il<eleAlgorithmNames_.size(); il++) {
    muonIsRA2[muonAlgorithmNames_[il]].clear();
    muonIsGlobalMuon[muonAlgorithmNames_[il]].clear();
    muonIsGlobalMuonPromptTight[muonAlgorithmNames_[il]].clear();
    muonPt[muonAlgorithmNames_[il]].clear();
    muonEta[muonAlgorithmNames_[il]].clear();
    muonPhi[muonAlgorithmNames_[il]].clear();
    muonTrackIso[muonAlgorithmNames_[il]].clear();
    muonEcalIso[muonAlgorithmNames_[il]].clear();
    muonHcalIso[muonAlgorithmNames_[il]].clear();

    muonEoverP[muonAlgorithmNames_[il]].clear();
    muonPtRatio[muonAlgorithmNames_[il]].clear();

    muonChi2[muonAlgorithmNames_[il]].clear();
    muonNdof[muonAlgorithmNames_[il]].clear();
    muonNhits[muonAlgorithmNames_[il]].clear();

    muonTrackd0[muonAlgorithmNames_[il]].clear();
    muonTrackPhi[muonAlgorithmNames_[il]].clear();

    muonPassID[muonAlgorithmNames_[il]].clear();

    muonVtx_z[muonAlgorithmNames_[il]].clear();

    eleEt[eleAlgorithmNames_[il]].clear();
    //    elePt[eleAlgorithmNames_[il]].clear();
    eleEta[eleAlgorithmNames_[il]].clear();
    elePhi[eleAlgorithmNames_[il]].clear();
    eleTrackIso[eleAlgorithmNames_[il]].clear();
    eleEcalIso[eleAlgorithmNames_[il]].clear();
    eleHcalIso[eleAlgorithmNames_[il]].clear();

    eledB[eleAlgorithmNames_[il]].clear();

    eleVtx_z[eleAlgorithmNames_[il]].clear();

    eleIDLoose[eleAlgorithmNames_[il]].clear();
    eleIDRobustTight[eleAlgorithmNames_[il]].clear();
    elePassID[eleAlgorithmNames_[il]].clear();
    eleIsRA2[eleAlgorithmNames_[il]].clear();
  }
  for (unsigned int il=0; il<tauAlgorithmNames_.size(); il++) {
    tauPt[tauAlgorithmNames_[il]].clear();
    tauEta[tauAlgorithmNames_[il]].clear();
    tauPhi[tauAlgorithmNames_[il]].clear();
    tauTaNC[tauAlgorithmNames_[il]].clear();

    for (unsigned int ialg=0; ialg<tauidAlgorithmNames_.size(); ialg++) {
      tauID[tauAlgorithmNames_[il]][tauidAlgorithmNames_[ialg]].clear();
    }
  }

  pv_pass=false;
  pv_isFake.clear();
  pv_z.clear();
  pv_ndof.clear();
  pv_chi2.clear();
  pv_rho.clear();


  for (unsigned int im=0; im<metAlgorithmTags_.size() ; im++) {
    MET[metAlgorithmTags_[im]]=-99;
    METphi[metAlgorithmTags_[im]]=-99;
    METsig[metAlgorithmTags_[im]]=-99;

    GenMET[metAlgorithmTags_[im]]=-99;
    GenMETphi[metAlgorithmTags_[im]]=-99;
  }
  trackPt.clear();
  trackEta.clear();
  trackPhi.clear();
  SUSY_nb=0;
  qScale=0;
  mcWeight=0;
  ptHat=0;
  SumPtOverHT=0;
  pdfWeights.clear(); 
  if (isMC_) pdfWeights.reserve(45); //this number is hard-coded twice in this class FIXME
  topDecayCode.clear();
  ZDecayMode=-99;
  flavorHistory=-99;

  susy_tanBeta=-99;
  susy_A0=-99;
  susy_crossSection=-99;
  susy_m0=-99;
  susy_m12=-99;
  susy_mu=-99;
  susy_run=-99;

  runNumber=0;
  eventNumber=0;
  lumiSection=0;

}

// ------------ method called to for each event  ------------
void
BasicTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  
  resetTreeVariables();

  //  cout<<" ==== "<<endl;
  Heventcount_->Fill(0.5);
  

  //start analyzin'
  //  cutResults.push_back(true); //"Inclusive" step

  //event info
  runNumber   = iEvent.run();
  eventNumber = iEvent.eventAuxiliary().event() ;
  lumiSection = iEvent.getLuminosityBlock().luminosityBlock();

  //  std::cout<<"~~~~~ "<<runNumber<<" "<<lumiSection<<" "<<eventNumber<<std::endl; //debug!

  // todo -- should put this in its own method, and check out sal's code -- much fancier!
  //Get the beam spot
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
  reco::BeamSpot beamSpot = *beamSpotHandle;
  bsx = beamSpot.x0();
  bsy = beamSpot.y0();
  bsz = beamSpot.z0();

  //trigger
  fillTriggerInfo(iEvent,iSetup);

  //primary vertex

  pv_pass= pvSelector_( iEvent ) ;
  fillPVInfo(iEvent,iSetup);

  //fills jets cuts, as well as MET, and the values of DeltaPhi
  for (unsigned int i=0; i<jetAlgorithmNames_.size(); i++) {
    fillJetInfo( iEvent,iSetup, i);
  }
  fillMetInfo( iEvent,iSetup);

  /* for now we don't bother to do this */
  //fill b efficiency info (must be done after jets are filled)
  //fillBTagDBInfo(iEvent,iSetup);

  //fill muon, electron
  for (unsigned int i=0; i<eleAlgorithmNames_.size(); i++) {
    fillLeptonInfo(iEvent,iSetup, i);
  }
  //fill tau
  for (unsigned int i=0; i<tauAlgorithmNames_.size(); i++) {
    fillTauInfo(iEvent,iSetup, i);
  }


  fillTrackInfo(iEvent,iSetup);
  if (isMC_)    fillMCInfo(iEvent,iSetup);

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

// ------------ method called once per run  ------------
void
BasicTreeMaker::beginRun(const edm::Run & iRun, const edm::EventSetup & iSetup)
{
  if (doPrescale_) {
  //compiles but does not run (at least on a test MC PAT tuple)
  bool hltChanged=false;
  if (  hltConfig_.init(iRun, iSetup, "HLT", hltChanged) ) {
    if (hltChanged) {
      std::cout<<"beginRun: The HLT config has changed!"<<std::endl;
    }
  }
  else {
    std::cout<<"ERROR -- something went wrong in hltConfig.init()!"<<std::endl;
  }
  }
}

void BasicTreeMaker::findHLTProcessName(const edm::Event& iEvent) {
  using namespace edm;

  //run this only once
  if (processName_!="") return;

  //some code from Josh Bendavid to automagically get the trigger tag (in MC)
  //https://hypernews.cern.ch/HyperNews/CMS/get/physTools/1791/1/1/1/1/1/2.html

  Handle<trigger::TriggerEvent> triggerEventHLT;
  iEvent.getByLabel("hltTriggerSummaryAOD", triggerEventHLT);

  processName_ = triggerEventHLT.provenance()->processName();
  std::cout<<"[BasicTreeMaker::findHLTProcessName] Will use trigger tag = "<<processName_<<std::endl;

}

// ------------ method called once each job just before starting event loop  ------------
void 
BasicTreeMaker::beginJob()
{
  using namespace std;
  //nfail_=0;

  //== create tree branches ==

  //tree that is filled only once per job
  infotree_->Branch("btagAlgorithms",&btagAlgorithmNames_);
  infotree_->Branch("tauidAlgorithms",&tauidAlgorithmNames_);
  infotree_->Branch("triggerList",&triggersOfInterest_);

  //tree that is filled for each event
  tree_->Branch("runNumber",&runNumber,"runNumber/l");
  tree_->Branch("lumiSection",&lumiSection,"lumiSection/l");
  tree_->Branch("eventNumber",&eventNumber,"eventNumber/l");

  tree_->Branch("bsx",&bsx,"bsx/F");
  tree_->Branch("bsy",&bsy,"bsy/F");
  tree_->Branch("bsz",&bsz,"bsz/F");

  //  tree_->Branch("cutResults",&cutResults);
  tree_->Branch("passTrigger",&passTrigger); //stores results for triggers listed in triggerList
  tree_->Branch("hltPrescale",&hltPrescale);
  //primary vertex branches
  tree_->Branch("pv_pass",&pv_pass,"pv_pass/O"); //output of PVSelector
  tree_->Branch("pv_isFake",&pv_isFake);
  tree_->Branch("pv_z",&pv_z);
  tree_->Branch("pv_rho",&pv_rho);
  tree_->Branch("pv_chi2",&pv_chi2);
  tree_->Branch("pv_ndof",&pv_ndof);

  //muon cleaning
  tree_->Branch("muonPtDiff_PF",&muonPtDiff_PF);
  tree_->Branch("passesBadPFMuonFilter",&passesBadPFMuonFilter,"passesBadPFMuonFilter/O");
  tree_->Branch("passesInconsistentMuonPFCandidateFilter",&passesInconsistentMuonPFCandidateFilter,"passesInconsistentMuonPFCandidateFilter/O");

  //jet branches
  for (vector<string>::const_iterator ij=jetAlgorithmTags_.begin() ; ij!=jetAlgorithmTags_.end() ; ++ij) {

    string tail="_";
    tail += (*ij);

    tree_->Branch( (string("tightJetIndex")+tail).c_str(),&tightJetIndex[*ij]);
    tree_->Branch( (string("looseJetIndex")+tail).c_str(),&looseJetIndex[*ij]);

    //loose jet branches
    tree_->Branch( (string("loosejetPt")+tail).c_str(),&loosejetPt[*ij]);
    tree_->Branch( (string("loosejetPz")+tail).c_str(),&loosejetPz[*ij]);
    tree_->Branch( (string("loosejetPtUncorr")+tail).c_str(),&loosejetPtUncorr[*ij]);
    tree_->Branch( (string("loosejetEt")+tail).c_str(),&loosejetEt[*ij]);
    tree_->Branch( (string("loosejetE")+tail).c_str(),&loosejetE[*ij]);
    tree_->Branch( (string("loosejetEta")+tail).c_str(),&loosejetEta[*ij]);
    tree_->Branch( (string("loosejetPhi")+tail).c_str(),&loosejetPhi[*ij]);
    tree_->Branch( (string("loosejetPassLooseID")+tail).c_str(),&loosejetPassLooseID[*ij]);
    tree_->Branch( (string("loosejetPassTightID")+tail).c_str(),&loosejetPassTightID[*ij]);
    tree_->Branch( (string("loosejetEnergyFracHadronic")+tail).c_str(),&loosejetEnergyFracHadronic[*ij]);

    tree_->Branch( (string("loosejetJECUncPlus")+tail).c_str(),&loosejetJECUncPlus[*ij]);
    tree_->Branch( (string("loosejetJECUncMinus")+tail).c_str(),&loosejetJECUncMinus[*ij]);

    tree_->Branch( (string("loosejetFlavor")+tail).c_str(),&loosejetFlavor[*ij]);
    
    tree_->Branch( (string("loosejetGenPt")+tail).c_str(),&loosejetGenPt[*ij]);
    tree_->Branch( (string("loosejetGenEta")+tail).c_str(),&loosejetGenEta[*ij]);
    tree_->Branch( (string("loosejetGenPhi")+tail).c_str(),&loosejetGenPhi[*ij]);

    tree_->Branch( (string("loosejetGenParticlePDGId")+tail).c_str(),&loosejetGenParticlePDGId[*ij]);
    tree_->Branch( (string("loosejetInvisibleEnergy")+tail).c_str(),&loosejetInvisibleEnergy[*ij]);

    tree_->Branch( (string("loosejetNTracks")+tail).c_str(),&loosejetNTracks[*ij]);
    tree_->Branch( (string("loosejetNSV")+tail).c_str(),&loosejetNSV[*ij]);
    tree_->Branch( (string("loosejetSVUnWeightedMass")+tail).c_str(),&loosejetSVUnWeightedMass[*ij]);
    tree_->Branch( (string("loosejetSVWeightedMass")+tail).c_str(),&loosejetSVWeightedMass[*ij]);
    tree_->Branch( (string("loosejetSVUnWeightedLifetime")+tail).c_str(),&loosejetSVUnWeightedLifetime[*ij]);
    tree_->Branch( (string("loosejetSVWeightedLifetime")+tail).c_str(),&loosejetSVWeightedLifetime[*ij]);
    tree_->Branch( (string("loosejetSVUnWeightedCosTheta")+tail).c_str(),&loosejetSVUnWeightedCosTheta[*ij]);
    tree_->Branch( (string("loosejetSVWeightedCosTheta")+tail).c_str(),&loosejetSVWeightedCosTheta[*ij]);

    //ROOT does not allow us to put a map of vectors into a tree
    //so we have to go with this approach
    for (unsigned int ib=0; ib<btagAlgorithmNames_.size(); ib++) {
      string bname = "loosejetBTagDisc_";
      bname += btagAlgorithmNames_[ib];
      bname += tail;
      tree_->Branch(bname.c_str(),&loosejetBTagDisc[*ij][btagAlgorithmNames_[ib]]);
    }
   
  }

  for (unsigned int im=0; im<metAlgorithmTags_.size(); im++) {
    string metname=metAlgorithmTags_[im];
    metname += "MET";
    string thirdarg=metname;
    thirdarg += "/F";
    tree_->Branch(metname.c_str(),&MET[metAlgorithmTags_[im]],thirdarg.c_str());

    metname=metAlgorithmTags_[im];
    metname += "METphi";
    thirdarg=metname;
    thirdarg += "/F";
    tree_->Branch(metname.c_str(),&METphi[metAlgorithmTags_[im]],thirdarg.c_str());

    metname=metAlgorithmTags_[im];
    metname += "METsig";
    thirdarg=metname;
    thirdarg += "/F";
    tree_->Branch(metname.c_str(),&METsig[metAlgorithmTags_[im]],thirdarg.c_str());

    metname=metAlgorithmTags_[im];
    metname += "GenMET";
    thirdarg=metname;
    thirdarg += "/F";
    tree_->Branch(metname.c_str(),&GenMET[metAlgorithmTags_[im]],thirdarg.c_str());

    metname=metAlgorithmTags_[im];
    metname += "GenMETphi";
    thirdarg=metname;
    thirdarg += "/F";
    tree_->Branch(metname.c_str(),&GenMETphi[metAlgorithmTags_[im]],thirdarg.c_str());
  }

  tree_->Branch("trackPt",&trackPt);
  tree_->Branch("trackEta",&trackEta);
  tree_->Branch("trackPhi",&trackPhi);
  tree_->Branch("SumPtOverHT",&SumPtOverHT,"SumPtOverHT/F");

  for (unsigned int il=0; il<tauAlgorithmNames_.size(); il++) {
    string tail="_";
    //string itail="_"; //only need this for integer fields
    if (tauAlgorithmNames_[il].find( "PF")!=string::npos) {
      tail+="PF";
      //      itail+="PF/I";
    }
    else {
      tail="";
      //      itail="/I";
    }

    tree_->Branch( (string("tauPt")+tail).c_str(),&tauPt[tauAlgorithmNames_[il]]);
    tree_->Branch( (string("tauEta")+tail).c_str(),&tauEta[tauAlgorithmNames_[il]]);
    tree_->Branch( (string("tauPhi")+tail).c_str(),&tauPhi[tauAlgorithmNames_[il]]);
    tree_->Branch( (string("tauTaNC")+tail).c_str(),&tauTaNC[tauAlgorithmNames_[il]]);

    for (unsigned int ialg=0; ialg<tauidAlgorithmNames_.size(); ialg++) {
      string bname = "tauID_";
      bname += tauidAlgorithmNames_[ialg];
      bname += tail;
      tree_->Branch(bname.c_str(),&tauID[tauAlgorithmNames_[il]][tauidAlgorithmNames_[ialg]]);
    }
  }

  for (unsigned int il=0; il<muonAlgorithmNames_.size(); il++) {
    string tail="_";
    string itail="_";
    if (muonAlgorithmNames_[il].find( "PF")!=string::npos) {
      tail+="PF";
      itail+="PF/I";
    }
    else {
      tail="";
      itail="/I";
    }

    tree_->Branch( (string("muonIsRA2")+tail).c_str(),&muonIsRA2[muonAlgorithmNames_[il]]);
    tree_->Branch( (string("muonIsGlobalMuon")+tail).c_str(),&muonIsGlobalMuon[muonAlgorithmNames_[il]]);
    tree_->Branch( (string("muonIsGlobalMuonPromptTight")+tail).c_str(),&muonIsGlobalMuonPromptTight[muonAlgorithmNames_[il]]);

    tree_->Branch( (string("muonEoverP")+tail).c_str(),&muonEoverP[muonAlgorithmNames_[il]]);
    tree_->Branch( (string("muonPtRatio")+tail).c_str(),&muonPtRatio[muonAlgorithmNames_[il]]);

    tree_->Branch( (string("muonPt")+tail).c_str(),&muonPt[muonAlgorithmNames_[il]]);
    tree_->Branch( (string("muonEta")+tail).c_str(),&muonEta[muonAlgorithmNames_[il]]);
    tree_->Branch( (string("muonPhi")+tail).c_str(),&muonPhi[muonAlgorithmNames_[il]]);
    tree_->Branch( (string("muonTrackIso")+tail).c_str(),&muonTrackIso[muonAlgorithmNames_[il]]);
    tree_->Branch( (string("muonEcalIso")+tail).c_str(),&muonEcalIso[muonAlgorithmNames_[il]]);
    tree_->Branch( (string("muonHcalIso")+tail).c_str(),&muonHcalIso[muonAlgorithmNames_[il]]);

    tree_->Branch( (string("muonChi2")+tail).c_str(),&muonChi2[muonAlgorithmNames_[il]]);
    tree_->Branch( (string("muonNdof")+tail).c_str(),&muonNdof[muonAlgorithmNames_[il]]);
    tree_->Branch( (string("muonNhits")+tail).c_str(),&muonNhits[muonAlgorithmNames_[il]]);
    tree_->Branch( (string("muonTrackd0")+tail).c_str(),&muonTrackd0[muonAlgorithmNames_[il]]);
    tree_->Branch( (string("muonTrackPhi")+tail).c_str(),&muonTrackPhi[muonAlgorithmNames_[il]]);
    tree_->Branch( (string("muonPassID")+tail).c_str(),&muonPassID[muonAlgorithmNames_[il]]);

    tree_->Branch( (string("muonVtx_z")+tail).c_str(),&muonVtx_z[muonAlgorithmNames_[il]]);


    tree_->Branch((string("eleIsRA2")+tail).c_str(),&eleIsRA2[eleAlgorithmNames_[il]]);
    tree_->Branch((string("eleEt")+tail).c_str(),&eleEt[eleAlgorithmNames_[il]]);
    //    tree_->Branch((string("elePt")+tail).c_str(),&elePt[eleAlgorithmNames_[il]]);
    tree_->Branch((string("eleEta")+tail).c_str(),&eleEta[eleAlgorithmNames_[il]]);
    tree_->Branch((string("elePhi")+tail).c_str(),&elePhi[eleAlgorithmNames_[il]]);
    tree_->Branch((string("eleTrackIso")+tail).c_str(),&eleTrackIso[eleAlgorithmNames_[il]]);
    tree_->Branch((string("eleEcalIso")+tail).c_str(),&eleEcalIso[eleAlgorithmNames_[il]]);
    tree_->Branch((string("eleHcalIso")+tail).c_str(),&eleHcalIso[eleAlgorithmNames_[il]]);

    tree_->Branch((string("eledB")+tail).c_str(),&eledB[eleAlgorithmNames_[il]]);

    tree_->Branch((string("eleVtx_z")+tail).c_str(),&eleVtx_z[eleAlgorithmNames_[il]]);

    tree_->Branch((string("eleIDLoose")+tail).c_str(),&eleIDLoose[eleAlgorithmNames_[il]]);
    tree_->Branch((string("eleIDRobustTight")+tail).c_str(),&eleIDRobustTight[eleAlgorithmNames_[il]]);
    tree_->Branch((string("elePassID")+tail).c_str(),&elePassID[eleAlgorithmNames_[il]]);
  }

  tree_->Branch("SUSY_nb",&SUSY_nb,"SUSY_nb/I");
  tree_->Branch("qScale",&qScale,"qScale/D");
  tree_->Branch("ptHat",&ptHat,"ptHat/D");
  tree_->Branch("mcWeight",&mcWeight,"mcWeight/D");
  tree_->Branch("pdfWeights",&pdfWeights);
  tree_->Branch("topDecayCode",&topDecayCode);
  tree_->Branch("ZDecayMode",&ZDecayMode,"ZDecayMode/I");
  tree_->Branch("flavorHistory",&flavorHistory,"flavorHistory/I");

  tree_->Branch("susy_tanBeta",&susy_tanBeta,"susy_tanBeta/D");
  tree_->Branch("susy_A0",&susy_A0,"susy_A0/D");
  tree_->Branch("susy_crossSection",&susy_crossSection,"susy_crossSection/D");
  tree_->Branch("susy_m0",&susy_m0,"susy_m0/D");
  tree_->Branch("susy_m12",&susy_m12,"susy_m12/D");
  tree_->Branch("susy_mu",&susy_mu,"susy_mu/D");
  tree_->Branch("susy_run",&susy_run,"susy_run/D");

  thetimer_->start();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BasicTreeMaker::endJob() {
  using namespace std;

  thetimer_->stop();

  infotree_->Fill();

  cout<<"Job took [CPU/Wall] seconds: "<<thetimer_->cpuTime()<<" / "<<thetimer_->realTime()<<endl;
  cout<<"events/sec = "<<Heventcount_->GetEntries()<<" / "<<thetimer_->realTime()<<" = "<<Heventcount_->GetEntries() /thetimer_->realTime()<<endl;  


}

//define this as a plug-in
DEFINE_FWK_MODULE(BasicTreeMaker);
