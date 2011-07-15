// -*- C++ -*-
//
// Package:    MPT
// Class:      MPT
// 
/**\class MPT MPT.cc Ntuples/MPT/src/MPT.cc

 Description: simple calculation of the track-based missing momentum

 Implementation:
     results of the calculation are stored in a class SimpleEt that simply holds two floats (pt and phi)
*/
//
// Original Author:  Joshua Thompson,6 R-029,+41227678914,
//         Created:  Thu Jul 14 15:14:59 CEST 2011
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "MyDataFormats/RA2bDataFormats/interface/SimpleEt.h"
//
// class declaration
//

class MPT : public edm::EDProducer {
   public:
      explicit MPT(const edm::ParameterSet&);
      ~MPT();

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

  double minpt_;

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
MPT::MPT(const edm::ParameterSet& iConfig) :
  minpt_ (iConfig.getParameter<double>("minpt"))
{
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/

  produces< SimpleEt >();

   //now do what ever other initialization is needed
  
}


MPT::~MPT()
{
 
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MPT::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<reco::TrackCollection> generalTracks;
  iEvent.getByLabel("generalTracks", generalTracks);

  float mptx=0,mpty=0;

  for (reco::TrackCollection::const_iterator trit=generalTracks->begin(); trit!=generalTracks->end() ; ++trit) {

    if (trit->pt() >= minpt_  && fabs(trit->eta()) < 5) {
      mptx -= trit->pt() * cos(trit->phi());
      mpty -= trit->pt() * sin(trit->phi());
    }
  }

  //std::cout<<"MPT["<<int(minpt_)<<"] = "<<sqrt(mptx*mptx + mpty*mpty)  <<std::endl;

  std::auto_ptr<SimpleEt> thempt(new SimpleEt(-99));
  thempt->pt = sqrt(mptx*mptx + mpty*mpty) ;
  thempt->phi = atan2(mpty,mptx);

  iEvent.put(thempt);

}

// ------------ method called once each job just before starting event loop  ------------
void 
MPT::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MPT::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
MPT::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
MPT::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MPT::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MPT::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MPT::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MPT);
