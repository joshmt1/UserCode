// -*- C++ -*-
//
// Package:    FindEventNumber
// Class:      FindEventNumber
// 
/**\class FindEventNumber FindEventNumber.cc CUSUSY/FindEventNumber/src/FindEventNumber.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Joshua Thompson,6 R-029,+41227678914,
//         Created:  Sat Oct 30 21:09:21 CEST 2010
// $Id: FindEventNumber.cc,v 1.1 2010/11/26 09:06:08 joshmt Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class FindEventNumber : public edm::EDFilter {
   public:
      explicit FindEventNumber(const edm::ParameterSet&);
      ~FindEventNumber();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------

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
FindEventNumber::FindEventNumber(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

}


FindEventNumber::~FindEventNumber()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
FindEventNumber::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   unsigned int   runNumber   = iEvent.run();
   //unsigned int   eventNumber = iEvent.eventAuxiliary().event() ;
   //   unsigned int   lumiSection = iEvent.getLuminosityBlock().luminosityBlock();
   
   //hard code a list of event numbers!
   //   if ( eventNumber == 817 ) return true;
   //   if ( eventNumber == 25463 ) return true;

   if (runNumber == 148032) return true;

   return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
FindEventNumber::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FindEventNumber::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(FindEventNumber);
