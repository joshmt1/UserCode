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
// $Id: TTbarDecayCoder.cc,v 1.2 2012/11/11 18:08:50 joshmt Exp $
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
  int findTopDecayMode( const reco::Candidate & cand,float& taupt, float& tauvisiblept, float &taueta, float &tauphi, int &nChargedTauDaughters) ;
  void analyzeTauDecays(const reco::Candidate *mytau,bool &found_tau_elec,bool & found_tau_muon,bool & found_tau_had, int & nChargedTauDaughters, float &tauvisiblept);

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
  produces< std::vector<int> >("topDecayCode");
  produces< std::vector<float> >("tauGenPt");
  produces< std::vector<float> >("tauGenVisPt");
  produces< std::vector<int> >("tauGenNProng");
  produces< std::vector<float> >("tauGenEta");
  produces< std::vector<float> >("tauGenPhi");
  
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
  //  int topDecayCode[2]={0,0};

  std::auto_ptr< std::vector<int> >  topDecayCode(new std::vector<int>);
  std::auto_ptr< std::vector<float> > tauGenPt(new std::vector<float>);
  std::auto_ptr< std::vector<float> > tauGenVisPt(new std::vector<float>);
  std::auto_ptr< std::vector<int> > tauGenNProng(new std::vector<int>);
  std::auto_ptr< std::vector<float> > tauGenEta(new std::vector<float>);
  std::auto_ptr< std::vector<float> > tauGenPhi(new std::vector<float>);

  for (size_t k = 0 ; k<genParticles->size(); k++) {
    const reco::Candidate & TCand = (*genParticles)[ k ];
    if (abs(TCand.pdgId())== 6 && TCand.status()==3) { //find t quark

      float taupt=-1;
      float tauptvis=-1;
      float taueta=-99;
      float tauphi=-99;
      int nprong=-1;

      int topcode = findTopDecayMode(TCand,taupt,tauptvis,taueta,tauphi,nprong);
      //      std::cout<<"[main] "<<topcode<<std::endl;
      if (ntops<=1) {
	topDecayCode->push_back(topcode);
	tauGenPt->push_back(taupt);
	tauGenNProng->push_back(nprong);
	tauGenVisPt->push_back(tauptvis);
	tauGenEta->push_back(taueta);
	tauGenPhi->push_back(tauphi);
	ntops++;
      }
      else { //should not happen
	std::cout<<"ntops = "<<ntops<<std::endl;
      }
    }
  }

  std::auto_ptr< int >  ttbarDecayCode(new int);

  *ttbarDecayCode =   getTopDecayCategory(topDecayCode->at(0), topDecayCode->at(1));

  iEvent.put(ttbarDecayCode,"ttbarDecayCode");
  iEvent.put(topDecayCode,"topDecayCode");
  iEvent.put(tauGenPt,"tauGenPt");
  iEvent.put(tauGenVisPt,"tauGenVisPt");
  iEvent.put(tauGenNProng,"tauGenNProng");
  iEvent.put(tauGenEta,"tauGenEta");
  iEvent.put(tauGenPhi,"tauGenPhi");

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

void TTbarDecayCoder::analyzeTauDecays(const reco::Candidate *mytau,bool &found_tau_elec,bool & found_tau_muon,bool & found_tau_had, int & nChargedTauDaughters, float &tauvisiblept) {

  //    std::cout<<"[analyzeTauDecays] "<<mytau->pdgId()<<std::endl; //debug info

  int nrecursivecalls=0;

  for (size_t m=0; m<mytau->numberOfDaughters(); m++) {
    //         std::cout<<"\t"<<mytau->daughter(m)->pdgId()<<"\t"<<mytau->daughter(m)->status()<<std::endl;
    int newtaudau = abs(mytau->daughter(m)->pdgId());
    
    if (newtaudau == 11) { found_tau_elec=true; nChargedTauDaughters++; tauvisiblept += mytau->daughter(m)->pt();}
    else if (newtaudau == 13) {found_tau_muon=true; nChargedTauDaughters++; tauvisiblept += mytau->daughter(m)->pt();}
    else if (newtaudau ==213 || newtaudau==111 || newtaudau==113 ||newtaudau==211||newtaudau==130
	     || newtaudau==223 ||newtaudau==321 ||newtaudau ==323 ||newtaudau ==310||newtaudau==221 ||newtaudau==20213) {
      found_tau_had=true;
      // nChargedTauDaughters code
      unsigned int ntaud = mytau->daughter(m)->numberOfDaughters();
      if (ntaud == 0) {
	//this catches cases of tau -> pi nu
	if ( mytau->daughter(m)->charge() != 0)  nChargedTauDaughters++;
	if ( abs(mytau->daughter(m)->pdgId()) != 14 && abs(mytau->daughter(m)->pdgId()) != 12 && abs(mytau->daughter(m)->pdgId()) != 16)  tauvisiblept += mytau->daughter(m)->pt();
      }
      else { //this catches cases where the tau does something fancier
	for (size_t mm=0;	mm< ntaud ; mm++) {
	  //			std::cout<<"\t\t"<<mytau2->daughter(l)->daughter(m)->pdgId()<<"\t"<<mytau2->daughter(l)->daughter(m)->status()<<std::endl;
	  if (mytau->daughter(m)->daughter(mm)->charge() != 0) nChargedTauDaughters++;

	  unsigned int apdgid = abs(mytau->daughter(m)->daughter(mm)->pdgId());
	  if ( apdgid != 14 &&  apdgid!= 12 &&  apdgid!= 16)  tauvisiblept += mytau->daughter(m)->daughter(mm)->pt();

	}
      }

    }
    else if (newtaudau == 14 || newtaudau==16 || newtaudau==12) { } //do nothing for neutrinos
    else if (newtaudau == 15  ) {  //the tau's daughter is a tau, so analyze the decays of *that* tau
      //            std::cout<<"[analyzeTauDecays] tau's daughter is tau"<<std::endl; //debug info
      nrecursivecalls++;
      analyzeTauDecays( mytau->daughter(m),found_tau_elec,found_tau_muon,found_tau_had, nChargedTauDaughters,tauvisiblept);
    }
    else if (newtaudau == 22  ) {  tauvisiblept += mytau->daughter(m)->pt();} //do nothing for photons, except sum up the visible pt
    else if (newtaudau == 24  ) { //weird case where the tau decays into a W
      //analyze the decays of the W
      //      std::cout<<"[analyzeTauDecays] tau's daughter is W"<<std::endl; //debug info
      nrecursivecalls++;
     analyzeTauDecays( mytau->daughter(m),found_tau_elec,found_tau_muon,found_tau_had, nChargedTauDaughters,tauvisiblept);
    }
    else std::cout<<"WARNING -- Unknown tau daughter found = "<<newtaudau<<std::endl; //debug info
  }

  if (nrecursivecalls>1) std::cout<<"WARNING -- nrecursions = "<<nrecursivecalls<<std::endl;

}

int TTbarDecayCoder::findTopDecayMode( const reco::Candidate & cand, float& taupt, float& tauvisiblept, float &taueta, float &tauphi, int &nChargedTauDaughters) {

  //the tau stuff was a specific hack for a specific set of code.
  //i'm going to leave it in,
  //the results won't go anywhere for now but they might be useful later
//   float taupt,taueta,tauphi,tauvisiblept;
//   int nChargedTauDaughters;

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
	    analyzeTauDecays( mytau,found_tau_elec, found_tau_muon, found_tau_had, nChargedTauDaughters,tauvisiblept);
	  }
	  else if (wdau == 14 || wdau==16 || wdau==12) { } //do nothing for neutrinos
	  else std::cout<<"WARNING -- W daughter is not as expected! "<<wdau<<std::endl;
	  
	}
	//skip neutrinos completely....
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
      std::cout<<"WARNING -- something weird happened in tau decay detection!"<<std::endl; //i have seen this before in the case of a 'decay' to pi+ e+e-
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
