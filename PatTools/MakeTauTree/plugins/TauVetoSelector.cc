#include "SandBox/Skims/src/IndirectTauVeto.cc"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
//#include "DataFormats/PatCandidates/interface/Tau.h"

#include <iostream>

class TauVetoSelector : public edm::EDFilter {
  
public:
  
  explicit TauVetoSelector(const edm::ParameterSet & iConfig);
  ~TauVetoSelector();
  
private:
  
  virtual bool filter(edm::Event & iEvent, const edm::EventSetup & iSetup);
  
  edm::InputTag jetSrc_, metSrc_;
  bool doTauVeto_;

};


TauVetoSelector::TauVetoSelector(const edm::ParameterSet & iConfig) {
  jetSrc_         = iConfig.getParameter<edm::InputTag>("JetSource");
  metSrc_         = iConfig.getParameter<edm::InputTag>("MetSource");
  doTauVeto_      = iConfig.getParameter<bool>("DoTauVeto");
  produces<std::vector<pat::Jet> >("");
}


TauVetoSelector::~TauVetoSelector() {
}


bool TauVetoSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // read in the objects
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(jetSrc_, jets);
  edm::Handle<edm::View<pat::MET> > mets;
  iEvent.getByLabel(metSrc_, mets);

  //output product
  std::auto_ptr<std::vector<pat::Jet> > prod(new std::vector<pat::Jet>());
  //loop over jets
  for (edm::View<pat::Jet>::const_iterator jj = jets->begin(); jj != jets->end(); ++jj) {

    //hard-coded for CSVT...not elegant at all
    bool notbjet = !(jj->bDiscriminator("combinedSecondaryVertexBJetTags") > 0.898);

    bool ptcut = jj->pt()>15;

    bool etacut = fabs( jj->eta() ) < 2.4;

    //does it NOT have a b tag, and does it have enough pt?
    if (notbjet && ptcut && etacut) {
      //if not, then see whether it is a tau
      std::vector<pat::Jet*> thisjet;
      pat::Jet * jet = new pat::Jet(*jj);
      thisjet.push_back( jet ); //correct?
      bool istau = tauVeto(  thisjet,mets->front().et() ,mets->front().phi());
      delete jet;
      
      if (istau) { //add this jet to the product
	prod->push_back(pat::Jet(*jj));
      }
    }
  }

  //  std::cout<<" Tau Veto -- "<<prod->size()<<std::endl;

  //special test!
//   std::cout<<" --- "<<std::endl;
//   edm::Handle<edm::View<pat::Tau> > taus;
//   iEvent.getByLabel("patTausPFchs", taus);
//   for (edm::View<pat::Jet>::const_iterator jet = jets->begin(); jet != jets->end(); ++jet) {
//     if (jet->pt()>15) std::cout<<"jet "<<jet->pt()<<"\t"<<jet->eta()<<"\t"<<jet->phi()<<std::endl;
//   }
//   for (edm::View<pat::Tau>::const_iterator tau = taus->begin(); tau != taus->end(); ++tau) {
//     if (tau->pt()>15) std::cout<<"tau "<<tau->pt()<<"\t"<<tau->eta()<<"\t"<<tau->phi()<<std::endl;
//   }

//another test
//   edm::Handle<edm::View<reco::PFMET> > pfmets;
//   iEvent.getByLabel("pfMet", pfmets);
//   std::cout<<" pat::MET = "<<mets->front().et()<<" ; "<<" reco::PFMET = "<<pfmets->front().et()<<std::endl;

  bool result = (doTauVeto_ ? (prod->size() == 0) : true);
  // store in the event
  iEvent.put(prod);
  return result;

}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(TauVetoSelector);
