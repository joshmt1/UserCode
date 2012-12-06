// -*- C++ -*-
//
// Package:    MinDeltaRToJet
// Class:      MinDeltaRToJet
// 
/**\class MinDeltaRToJet MinDeltaRToJet.cc PatTools/MinDeltaRToJet/src/MinDeltaRToJet.cc

 - given a list of jets and a PFCandidate, find the minDeltaR between the PF candidate and the jets in the list
[excluding the jet that contains that PF candidate]

 - produce a single float with that value

 -- OR --
   - given a list of jets and a genParticle, find the minDeltaR between the genParticle and the jets in the list.
   trick is to remove the reco'd version of the genParticle from the jet list
The only idea I have is to impose a lower limit on minDeltaR. Hard-coded to 0.3 for now

*/
//
// Original Author:  Joshua Thompson,6 R-029,+41227678914,
// $Id: MinDeltaRToJet.cc,v 1.3 2012/11/13 22:54:08 joshmt Exp $
//
//


// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaR.h"
//#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

//
// class declaration
//

class MinDeltaRToJet : public edm::EDProducer {
   public:
      explicit MinDeltaRToJet(const edm::ParameterSet&);
      ~MinDeltaRToJet();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  
      // ----------member data ---------------------------
  edm::InputTag candSrc_,genCandSrc_,jetSrc_;

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
MinDeltaRToJet::MinDeltaRToJet(const edm::ParameterSet& iConfig) :
  candSrc_(iConfig.getParameter<edm::InputTag>("pfCandSource")),
  genCandSrc_(iConfig.getParameter<edm::InputTag>("GenCandSource")),
  jetSrc_(iConfig.getParameter<edm::InputTag>("JetSource"))
{
  //register your products
  // Examples
  produces<float >("minDeltaRToJet");//.setBranchAlias("pfcands_trkiso"); //last part needed?
  
}


MinDeltaRToJet::~MinDeltaRToJet()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MinDeltaRToJet::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  const double tolerance = 0.01;
   
  edm::Handle< edm::View<reco::PFCandidate> > pfcand;
  iEvent.getByLabel(candSrc_,pfcand);

  edm::Handle< edm::View<reco::GenParticle> > gencand;
  iEvent.getByLabel(genCandSrc_,gencand);

  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(jetSrc_,jets);

  std::auto_ptr< float > minDeltaRToJet(new float);

  const bool usePf = !(pfcand.failedToGet());
  const size_t loopmax = usePf ? pfcand->size() : gencand->size() ;
  float mindr = 999;
  for (size_t k=0; k<loopmax; k++) {
    //if multiple elements in input list, we'll find the min over all of them
    //    if (k>0) {std::cout<<" [MinDeltaRToJet] cand list is bigger than one element"<<std::endl; break;}

    double muoneta,muonphi;
    if (usePf) {
      muoneta = (*pfcand)[k].eta();
      muonphi = (*pfcand)[k].phi();
    }
    else { //we call it muon but really it's a generic particle
      muoneta = (*gencand)[k].eta();
      muonphi = (*gencand)[k].phi();
    }
    //    std::cout<<"[muon] "<<(*pfcand)[k].pt()<<" "<<(*pfcand)[k].eta()<<" "<<(*pfcand)[k].phi()<<std::endl;

    for (edm::View<pat::Jet>::const_iterator jet = jets->begin(); jet != jets->end(); ++jet) {
      float dr = reco::deltaR( jet->eta(), jet->phi(), muoneta,muonphi);

      if ( dr < mindr && usePf) { //
	//loop over jet constituents to check if this jet is actually one and the same as the muon
	//we do this with a *very* tight DR and pT match
	//	std::cout<<"[jet] "<<jet->pt()<<" "<<jet->eta()<<" "<<jet->phi()<<std::endl;
	bool usethisjet=true;
	for (size_t ijc = 0; ijc<jet->getPFConstituents().size(); ijc++ ) {
	  reco::PFCandidate themuon = (*pfcand)[k];
	  //check DeltaR and pT. they should be darn near identical
	  if ( reco::deltaR(muoneta,muonphi, jet->getPFConstituents().at(ijc)->eta(),jet->getPFConstituents().at(ijc)->phi())<tolerance
	       && std::abs( jet->getPFConstituents().at(ijc)->pt() - themuon.pt())<tolerance) {
	    usethisjet=false;
	    break; //break out of the loop over pfconstituents
	  }
	}
	if (usethisjet)	mindr=dr;
      }
      else if (!usePf) {
	//here the logic is simple -- don't accept a dr less than cutoff of 0.3
	if ( dr < mindr && dr>= 0.3) mindr=dr;
      }
    }
  }

  *minDeltaRToJet = mindr;

  //   std::cout<<" mindr to jet = "<<mindr<<std::endl;

  iEvent.put(minDeltaRToJet,"minDeltaRToJet");
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
MinDeltaRToJet::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MinDeltaRToJet::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
MinDeltaRToJet::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
MinDeltaRToJet::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MinDeltaRToJet::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MinDeltaRToJet::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MinDeltaRToJet::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}



//define this as a plug-in
DEFINE_FWK_MODULE(MinDeltaRToJet);
