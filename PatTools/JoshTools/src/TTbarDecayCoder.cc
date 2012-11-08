// -*- C++ -*-
//
// Package:    TTbarDecayCoder
// Class:      TTbarDecayCoder
// 
/**\class TTbarDecayCoder TTbarDecayCoder.cc PatTools/TTbarDecayCoder/src/TTbarDecayCoder.cc

 Description:  a rather ugly tool for analyzing the genparticles and categorizing the decays of any tops
Should work for ttbar and T2tt. Not implemented for single top.


*/
//
// Original Author:  Joshua Thompson,6 R-029,+41227678914,
//         Created:  Wed Nov  7 14:21:41 CET 2012
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

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

//
// class declaration
//

class TTbarDecayCoder : public edm::EDProducer {
   public:
      explicit TTbarDecayCoder(const edm::ParameterSet&);
      ~TTbarDecayCoder();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  
  int getTopDecayCategory(int code1, int code2);
  int findTopDecayMode( const reco::Candidate & cand);// float& taupt, float& tauvisiblept, float &taueta, float &tauphi, int &nChargedTauDaughters) {

      // ----------member data ---------------------------
  edm::InputTag genParticleSrc_;

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
TTbarDecayCoder::TTbarDecayCoder(const edm::ParameterSet& iConfig) :
  genParticleSrc_(iConfig.getParameter<edm::InputTag>("GenParticleSource"))
{
  //register your products
  // Examples
  produces<int >("ttbarDecayCode");//.setBranchAlias("pfcands_trkiso"); //last part needed?
  
}


TTbarDecayCoder::~TTbarDecayCoder()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
TTbarDecayCoder::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(genParticleSrc_,genParticles);

  int ntops=0;
  int topDecayCode[2]={0,0};

  for (size_t k = 0 ; k<genParticles->size(); k++) {
    const reco::Candidate & TCand = (*genParticles)[ k ];
    if (abs(TCand.pdgId())== 6 && TCand.status()==3) { //find t quark
      int topcode = findTopDecayMode(TCand);//,taupt,tauptvis,taueta,tauphi,nprong);
      if (ntops<=1) {
	topDecayCode[ntops] = topcode;
	ntops++;
      }
      else { //should not happen
	std::cout<<"ntops = "<<ntops<<std::endl;
      }
    }
  }

  std::auto_ptr< int >  ttbarDecayCode(new int);

  *ttbarDecayCode =   getTopDecayCategory(topDecayCode[0], topDecayCode[1]);

  iEvent.put(ttbarDecayCode,"ttbarDecayCode");

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
TTbarDecayCoder::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TTbarDecayCoder::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
TTbarDecayCoder::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TTbarDecayCoder::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TTbarDecayCoder::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TTbarDecayCoder::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TTbarDecayCoder::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


int TTbarDecayCoder::getTopDecayCategory(int code1, int code2) {
  // copied from some old code:
  //http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/joshmt/CUSUSY/BasicTreeTools/basicLoop.h?revision=1.83&view=markup
  /*
    reclassifying one set of arbitrary integers with another set of arbitrary integers.
    ugly, but i don't have a better idea.
    i will use enums to make things a little more clear
  */

  //here are the codes defined in code1, code2
  enum TopDecayCodes {kUnknown=-1, kNoB = 0, kHadronic = 1, kElectron = 2, kMuon=3,kTauHadronic=4,kTauElectron=5,kTauMuon=6,kTauMisc=7 };
  enum TopDecayCategory {kTTbarUnknown=0,kAllLeptons=1,kAllHadronic=2,kOneElectron=3,kOneMuon=4,kOneTauE=5,kOneTauMu=6,kOneTauHadronic=7,kAllTau=8,kTauPlusLepton=9, nTopCategories=10};
  //note that this code is completely "fragile" against changes in the TopDecayCodes scheme.

  int code=-1;
  unsigned int ntop=0;
  if (code1>0) ++ntop;
  if (code2>0) ++ntop;

  if (ntop==2) { //this is what we expect for ttbar

    //i'm pretty sure that the 'misc' category is only filled by hadronic decays
    if (code1==kTauMisc) code1=kHadronic;
    if (code2==kTauMisc) code2=kHadronic;
    
    if ( (code1==kElectron ||code1==kMuon) && (code2==kElectron ||code2==kMuon) ) code=kAllLeptons;
    else if ( code1==kHadronic && code2==kHadronic ) code=kAllHadronic;
    else if ( (code1==kElectron && code2==kHadronic) || (code2==kElectron && code1==kHadronic) ) code=kOneElectron;
    else if ( (code1==kMuon && code2==kHadronic) || (code2==kMuon && code1==kHadronic) ) code=kOneMuon;
    else if ( (code1==kTauElectron && code2==kHadronic) || (code2==kTauElectron && code1==kHadronic) ) code=kOneTauE;
    else if ( (code1==kTauMuon && code2==kHadronic) || (code2==kTauMuon && code1==kHadronic) ) code=kOneTauMu;
    else if ( (code1==kTauHadronic && code2==kHadronic) || (code2==kTauHadronic && code1==kHadronic) ) code=kOneTauHadronic;
    //this logic depends on all of the highest categories being tau
    else if ( code1>=kTauHadronic && code2>=kTauHadronic ) code=kAllTau;
    else if ( code1>=kTauHadronic && (code2==kElectron || code2==kMuon) ) code=kTauPlusLepton;
    else if ( code2>=kTauHadronic && (code1==kElectron || code1==kMuon) ) code=kTauPlusLepton;
    else {
      code=kTTbarUnknown;
      std::cout<<"TopCode = "<<code1<<"\t"<<code2<<std::endl;
    }
  }
  //  else if (ntop==1) { //maybe we expect this for single top? //never cared enough to look
  //
  //  }
  //  else {
  //    
  //  }

  return code;
}


int TTbarDecayCoder::findTopDecayMode( const reco::Candidate & cand) {//, float& taupt, float& tauvisiblept, float &taueta, float &tauphi, int &nChargedTauDaughters) {

  //the tau stuff was a specific hack for a specific set of code.
  //i'm going to leave it in,
  //the results won't go anywhere for now but they might be useful later
  float taupt,taueta,tauphi,tauvisiblept;
  int nChargedTauDaughters;

  taupt = -1;
  tauvisiblept = 0;

  taueta=-99;
  tauphi=-99;
  nChargedTauDaughters = 0;

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
	    taupt = mytau->pt();
	    taueta = mytau->eta();
	    tauphi = mytau->phi();
	    //   	    std::cout<<" -- found a tau"<<std::endl;
	    for (size_t k=0; k<mytau->numberOfDaughters(); k++) {
	      //	      std::cout<<mytau->daughter(k)->pdgId()<<"\t"<<mytau->daughter(k)->status()<<std::endl;
	      const reco::Candidate *mytau2 = mytau->daughter(k);
	      if (abs(mytau2->pdgId()) == 15) { 
		
		for (size_t l=0; l<mytau2->numberOfDaughters(); l++) {
		  //	  		  std::cout<<"\t"<<mytau2->daughter(l)->pdgId()<<"\t"<<mytau2->daughter(l)->status()<<std::endl;

		  int taudau = abs(mytau2->daughter(l)->pdgId());
		  if (taudau == 11) { found_tau_elec=true; nChargedTauDaughters++;}
		  else if (taudau == 13) {found_tau_muon=true; nChargedTauDaughters++;}
		  else if (taudau ==213 || taudau==111 || taudau==113 ||taudau==211||taudau==130
			   || taudau==223 ||taudau==321 ||taudau ==323 ||taudau ==310||taudau==221 ||taudau==20213) {
		    found_tau_had=true;
		    unsigned int ntaud = mytau2->daughter(l)->numberOfDaughters();
		    if (ntaud == 0) {
		      //this catches cases of tau -> pi nu
		      if ( mytau2->daughter(l)->charge() != 0) nChargedTauDaughters++;
		    }
		    else { //this catches cases where the tau does something fancier
		      for (size_t m=0;	m< ntaud ; m++) {
			//			std::cout<<"\t\t"<<mytau2->daughter(l)->daughter(m)->pdgId()<<"\t"<<mytau2->daughter(l)->daughter(m)->status()<<std::endl;
			if (mytau2->daughter(l)->daughter(m)->charge() != 0) nChargedTauDaughters++;
		      }
		    }
		  }
		  else if (taudau == 14 || taudau==16 || taudau==12) { } //do nothing for neutrinos
		  else if (taudau == 22  ) { } //do nothing for photons 
		  else if (taudau == 24) { //weird case ... seems to always correspond to a tau->h decay; anyway put in the full machinery, more or less
		    unsigned int ntaud = mytau2->daughter(l)->numberOfDaughters();
		    for (size_t m=0;	m< ntaud ; m++) {
		      //		      std::cout<<"\t\t"<<mytau2->daughter(l)->daughter(m)->pdgId()<<"\t"<<mytau2->daughter(l)->daughter(m)->status()<<std::endl;
		      int Wdau = abs(mytau2->daughter(l)->daughter(m)->pdgId());
		      if (Wdau==11) found_tau_elec=true;
		      else if (Wdau==13) found_tau_muon=true;
		      else if (Wdau ==213 || Wdau==111 || Wdau==113 ||Wdau==211||Wdau==130
			   || Wdau==223 ||Wdau==321 ||Wdau ==323 ||Wdau ==310||Wdau==221 ||Wdau==20213) found_tau_had=true;

		      if (mytau2->daughter(l)->daughter(m)->charge() != 0) nChargedTauDaughters++;
		    }
		    
		  }
		  else if (taudau==15) { //did not see this in TTbar, but did see it in T2tt (madgraph v pythia?)
		    // std::cout<<"found tau as daughter of tau"<<std::endl; //debug info
		    //seems to be tau -> tau+photon
		    //again, copy and paste the code from above for classifying tau daughters
		    const reco::Candidate *newtau = mytau2->daughter(l);

		    for (size_t m=0; m<newtau->numberOfDaughters(); m++) {
		      //  std::cout<<"\t"<<newtau->daughter(m)->pdgId()<<"\t"<<newtau->daughter(m)->status()<<std::endl;
		      int newtaudau = abs(newtau->daughter(m)->pdgId());

		      //copy and paste of the code above...would be better to do this more elegantly!
		      if (newtaudau == 11) { found_tau_elec=true; nChargedTauDaughters++;}
		      else if (newtaudau == 13) {found_tau_muon=true; nChargedTauDaughters++;}
		      else if (newtaudau ==213 || newtaudau==111 || newtaudau==113 ||newtaudau==211||newtaudau==130
			       || newtaudau==223 ||newtaudau==321 ||newtaudau ==323 ||newtaudau ==310||newtaudau==221 ||newtaudau==20213) {
			found_tau_had=true;
			//did not put the nChargedTauDaughters code in here....hence that value will be wrong in this case
		      }
		      else if (newtaudau == 14 || newtaudau==16 || newtaudau==12) { } //do nothing for neutrinos
		      else if (newtaudau == 22  ) { } //do nothing for photons 
		      else std::cout<<"WARNING -- Unknown [new] tau daughter found = "<<newtaudau<<std::endl; //debug info

		    }
		    // std::cout<<"\t~~done~~"<<std::endl; //debug info

		  }
		  else {
		    std::cout<<"WARNING -- Unknown tau daughter found = "<<taudau<<std::endl; //debug info
		  }
		  

		  //now we want to sum the 'visible' pT of the tau.
		  //for now that means anything other than neutrinos.
		  if ( !(taudau == 14 || taudau==16 || taudau==12) ) { // NOT neutrinos
		    tauvisiblept += mytau2->daughter(l)->pt();
		  }

		}
		//		std::cout<<" nChargedTauDaughters = "<<nChargedTauDaughters<<std::endl;
	      }
	      else if (abs(mytau2->pdgId()) ==22) { //photon
		//do nothing except add to tau vis pt (should i do that?)
		tauvisiblept += mytau2->pt();
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


//define this as a plug-in
DEFINE_FWK_MODULE(TTbarDecayCoder);
