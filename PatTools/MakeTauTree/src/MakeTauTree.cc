// -*- C++ -*-
//
// Package:    MakeTauTree
// Class:      MakeTauTree
// 
/**\class MakeTauTree MakeTauTree.cc PatTools/MakeTauTree/src/MakeTauTree.cc

 Description: read PAT output and make a simple root tree
(like my "reducedTrees")

Goal is to study the decays of ttbar, in particular the tau rejection possibilities

Update -- I will now expand this code to try to put all sorts of useful quantities into the tree...

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Joshua Thompson,6 R-029,+41227678914,
//         Created:  Wed Jul 25 15:22:44 CEST 2012
// $Id: MakeTauTree.cc,v 1.3 2012/08/27 14:06:50 joshmt Exp $
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

   //added by me
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/METReco/interface/MET.h"
   //#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "DataFormats/Math/interface/deltaR.h"

//tau veto is precomputed but we want access to some of the pieces
//#include "SandBox/Skims/src/IndirectTauVeto.cc"
#include "StevenCode/LeptonSelection/src/IndirectTau.cc" //new location within steven's code

//ROOT includes
#include "TTree.h"

//
// class declaration
//

class MakeTauTree : public edm::EDAnalyzer {
   public:
      explicit MakeTauTree(const edm::ParameterSet&);
      ~MakeTauTree();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  void fillMCInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void doTauJetMatching(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  //copied and pasted (more or less) from code i wrote years ago....
  int findTopDecayMode( const reco::Candidate & cand, float& taupt, float& tauvisiblept, float &taueta, float &tauphi, int & nChargedTauDaughters);
  int getTopDecayCategory(int code1, int code2);

  void isolatedTrackVeto(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void isolatedTrackDetails(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void POGtaus(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  void clearVariables() ;

  // ----------member data ---------------------------

  // -- config info
  edm::InputTag pogtauSrc_,jetSrc_, metSrc_,electronSrc_,muonSrc_,tauSrc_,genParticleSrc_,vtxSrc_;
  bool storeAllJetInfo_;

  double matchradius_;

  // ntuple stuff

  TTree* tree_;

  //jet counters
  int njets70;
  int njets50;
  int njets30;
  int nCSVMbjets30;
  int nCSVTbjets30;

  //store info on all reco jets
  int njets10;
  float jetPt[100];
  float jetEta[100];
  float jetPhi[100];
  float jetCSV[100];
  float jetMT[100];
  int jetChm[100];
  float jetLrm[100];
  bool jetPassTauId[100];
  bool jetMatchGenTau[100];
  float jetTauDr[100]; //to study matching

  //event variables
  float MET;
  float METphi;
  int nVetoMuons;
  int nVetoElectrons;
  int nVetoTaus;
  int nIsolatedTracks;

  float maxPtIsoTrack_Iso;
  float maxPtIsoTrack_Pt;
  float minIsoIsoTrack_Iso;
  float minIsoIsoTrack_Pt;



  int nLooseTaus15;
  int nVLooseTaus15;
  int nMediumTaus15;
  int nLooseTaus20;
  int nVLooseTaus20;
  int nMediumTaus20;
  float DeltaPhiMetJet1;
  float DeltaPhiMetJet2;
  float DeltaPhiMetJet3;
  float minDeltaPhiMetBM;
  float minDeltaPhiMetBT;
  int nGoodPV;

  int lowChmJetChm;
  float lowChmJetLrm;
  float lowChmJetCSV;
  float lowChmJetMT;
  float lowChmJetPt;

  //some variables for DeltaPhi (b, MET)
  //info for top decays and tau decays
  int topDecayCode[2];
  int ttbarDecayCode;

  float tauGenPt[2];
  float tauGenVisPt[2];
  float tauGenEta[2];
  float tauGenPhi[2];
  int tauGenNProng[2];

  //info for reco jets that are matched to gen taus
  float taujetPt[2];
  float taujetEta[2];
  float taujetPhi[2];
  float taujetCSV[2];
  float taujetMT[2];
  int taujetChm[2];
  float taujetLrm[2];



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
MakeTauTree::MakeTauTree(const edm::ParameterSet& iConfig) :
  pogtauSrc_(iConfig.getParameter<edm::InputTag>("TauSource")), //these are the pat::Taus
  jetSrc_(iConfig.getParameter<edm::InputTag>("JetSource")),
  metSrc_(iConfig.getParameter<edm::InputTag>("MetSource")),
  electronSrc_(iConfig.getParameter<edm::InputTag>("ElectronVetoSource")),
  muonSrc_(iConfig.getParameter<edm::InputTag>("MuonVetoSource")),
  tauSrc_(iConfig.getParameter<edm::InputTag>("TauVetoSource")), //these are jets to be used in the indirect tau veto
  genParticleSrc_(iConfig.getParameter<edm::InputTag>("GenParticleSource")),
  vtxSrc_(iConfig.getParameter<edm::InputTag>("VertexSource")),
  storeAllJetInfo_(iConfig.getParameter<bool>("StoreAllJetInfo")),
  matchradius_(0.3),
  tree_(0)
{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  //branch creation for tree done in beginJob instead of here
  tree_=  fs->make<TTree>("reducedTree","reducedTree");

}


MakeTauTree::~MakeTauTree()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

void MakeTauTree::clearVariables() {

   taujetPt[0]=-99;
   taujetEta[0]=-99;
   taujetPhi[0]=-99;
   taujetCSV[0]=-99;
   taujetMT[0]=-99;
   taujetChm[0]=-99;
   taujetLrm[0]=-99;

   taujetPt[1]=-99;
   taujetEta[1]=-99;
   taujetPhi[1]=-99;
   taujetCSV[1]=-99;
   taujetMT[1]=-99;
   taujetChm[1]=-99;
   taujetLrm[1]=-99;

   njets10=-1;
   for (int i=0;i<100;i++) {
      jetPt[i]=-99;
      jetEta[i]=-99;
      jetPhi[i]=-99;
      jetCSV[i]=-99;
      jetMT[i]=-99;
      jetChm[i]=-99;
      jetLrm[i]=-99;
      jetTauDr[i]=-99;
      jetPassTauId[i]=false;
      jetMatchGenTau[i]=false;
   }
   
   DeltaPhiMetJet1=99;
   DeltaPhiMetJet2=99;
   DeltaPhiMetJet3=99;
   minDeltaPhiMetBM=99; //set higher than physical range
   minDeltaPhiMetBT=99; //set higher than physical range

}

// ------------ method called for each event  ------------
void
MakeTauTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //  std::cout<<" == event =="<<std::endl;

  clearVariables(); //got to be careful when we make a tree using global variables

  //order of the pieces here is fragile.
  //MET must come first, because other calculations depend on it
  edm::Handle<edm::View<reco::MET> > mets;
  iEvent.getByLabel(metSrc_, mets);
  MET = mets->front().et();
  METphi = mets->front().phi();

  //MC info is needed by some other stuff (gen tau info)...so it should come after MET
  fillMCInfo(iEvent,iSetup);  //info about decays of top, etc

  //do gen tau <--> reco jet matching //must come after fillMCInfo()
  doTauJetMatching(iEvent,iSetup);

  njets70=0;
  njets50=0;
  njets30=0;
  nCSVMbjets30=0;
  nCSVTbjets30=0;
  njets10=0;

  int njetsAllEta=0;

  lowChmJetChm =1000; //initialize big because we're looking for the min value
  lowChmJetLrm = -1;
  lowChmJetCSV = -1;
  lowChmJetMT = -1;
  lowChmJetPt =-1;

  edm::Handle<edm::View<pat::Jet> > jets; //jet collection had better be sorted by pT!
  iEvent.getByLabel(jetSrc_,jets);
  for (edm::View<pat::Jet>::const_iterator jet = jets->begin(); jet != jets->end(); ++jet) {

    //    std::cout<<"jet -- "<<jet->pt()<<"\t"<<jet->eta()<<std::endl;

    //loose eta cut for these
    if (jet->pt() >=10 && fabs(jet->eta())<5) {
      jetPt[njets10] = jet->pt();
      jetEta[njets10] = jet->eta();
      jetPhi[njets10]=jet->phi();
      jetCSV[njets10]=jet->bDiscriminator("combinedSecondaryVertexBJetTags");
      jetMT[njets10] =  sqrt(2*(jet->pt()*MET)*(1-cos( fabs(reco::deltaPhi(jet->phi(),METphi)) )));

      //      pat::Jet * thisjet = new pat::Jet(*jet); //just a trick to get the right data format
      jetChm[njets10] = getCHMdr(*jet,0.3,1);
      jetLrm[njets10] = getLRM(*jet);
      //      std::vector<pat::Jet*> jetv; jetv.push_back(thisjet);
      std::vector<pat::Jet> jetv; jetv.push_back(*jet);
      jetPassTauId[njets10] = tauVeto(jetv,MET,METphi) && (jet->pt()>15) && (jetCSV[njets10]<=0.898) && fabs(jet->eta())<2.4;
      //      delete thisjet;

      //finally, does the jet match a gen tau?
      double dr1 = reco::deltaR( jet->eta(), jet->phi(), tauGenEta[0],tauGenPhi[0]);
      double dr2 = reco::deltaR( jet->eta(), jet->phi(), tauGenEta[1],tauGenPhi[1]);
      bool match1 = dr1 < matchradius_;
      bool match2 = dr2 < matchradius_;
      jetMatchGenTau[njets10] = match1 || match2;

      jetTauDr[njets10] = (dr1<dr2) ? dr1 : dr2;
      
      //now look for the 'lowest chm' jet, with pT>15 and |eta|<2.4
      //this is just copying jet-by-jet info into an event-by-event quantity for convenience
      if (jet->pt() >=15 && fabs(jet->eta())<2.4) {
	if ( jetChm[njets10] < lowChmJetChm) {
	  lowChmJetChm = jetChm[njets10];
	  lowChmJetLrm = jetLrm[njets10];
	  lowChmJetCSV = jetCSV[njets10];
	  lowChmJetMT =  jetMT[njets10];
	  lowChmJetPt =  jetPt[njets10];
	}
      }
      
      
      njets10++;
    }


    //use all jets within |eta|<5. 2011 note says 4.7; that might have been an artifact of JEC problems
    //there is a nasty detail here about whether or not to apply jetID to these jets
    //logic says no, but jet collection might say yes!
    //moot if we veto any event that has a jet failing jet id (which is recommended in 2012)
    if (jet->pt() >=30 && fabs(jet->eta())<5) {
      ++njetsAllEta;
      if (njetsAllEta==1) DeltaPhiMetJet1 = reco::deltaPhi(jet->phi(),METphi);
      else if (njetsAllEta==2) DeltaPhiMetJet2 = reco::deltaPhi(jet->phi(),METphi);
      else if (njetsAllEta==3) DeltaPhiMetJet3 = reco::deltaPhi(jet->phi(),METphi);
    }

    if ( fabs(jet->eta()) < 2.4 ) {
      float pt = jet->pt();
      if ( pt >= 30) {
	++njets30;
	if (jet->bDiscriminator("combinedSecondaryVertexBJetTags") > 0.679) {
	  ++nCSVMbjets30;
	  if ( reco::deltaPhi(jet->phi(),METphi) < minDeltaPhiMetBM) minDeltaPhiMetBM=reco::deltaPhi(jet->phi(),METphi);
	}
	if (jet->bDiscriminator("combinedSecondaryVertexBJetTags") > 0.898) {
	  ++nCSVTbjets30;
	  if ( reco::deltaPhi(jet->phi(),METphi) < minDeltaPhiMetBT) minDeltaPhiMetBT=reco::deltaPhi(jet->phi(),METphi);
	}
	if (pt>=50) 	  ++njets50;
	if (pt>=70) 	  ++njets70;
      }
    }
  }

  //the first of these is the Lowette bare-bones implementation.
  //the 2nd is the direct hooberman version that gives a bit more info.
  isolatedTrackVeto(iEvent,iSetup);
  isolatedTrackDetails(iEvent,iSetup);
  POGtaus(iEvent,iSetup);

  edm::Handle<edm::View<reco::GsfElectron> > electrons;
  iEvent.getByLabel(electronSrc_, electrons);
  nVetoElectrons = electrons->size();

  edm::Handle<edm::View<reco::Muon> > muons;
  iEvent.getByLabel(muonSrc_, muons);
  nVetoMuons = muons->size();

  edm::Handle<edm::View<pat::Jet> > taus;
  iEvent.getByLabel(tauSrc_,taus);
  nVetoTaus = taus->size();

  edm::Handle<edm::View<reco::Vertex> > vertices;
  iEvent.getByLabel(vtxSrc_, vertices); //this one those be customized
  nGoodPV = vertices->size(); // now we know this is wrong. no quality cuts applied!

  //fill tree
  tree_->Fill();
}

int MakeTauTree::getTopDecayCategory(int code1, int code2) {
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


int MakeTauTree::findTopDecayMode( const reco::Candidate & cand, float& taupt, float& tauvisiblept, float &taueta, float &tauphi, int &nChargedTauDaughters) {

  //not sure if this is best way to do it, but let's see...
  //could generalize to other leptons
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
	    //	    std::cout<<" -- found a tau"<<std::endl;
	    for (size_t k=0; k<mytau->numberOfDaughters(); k++) {
	      //	      std::cout<<mytau->daughter(k)->pdgId()<<"\t"<<mytau->daughter(k)->status()<<std::endl;
	      const reco::Candidate *mytau2 = mytau->daughter(k);
	      if (abs(mytau2->pdgId()) == 15) { 
		
		for (size_t l=0; l<mytau2->numberOfDaughters(); l++) {
		  //		  std::cout<<"\t"<<mytau2->daughter(l)->pdgId()<<"\t"<<mytau2->daughter(l)->status()<<std::endl;

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

//good candidate for somehow being broken into its own file,
// to avoid all of this copying and pasting in the future
void MakeTauTree::fillMCInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //std::cout<<" == fillMCInfo =="<<std::endl;

  //snipping out genEventInfoProduct stuff, etc
  //focus only on top decay for now

  // === get some info about the decay structure in the event ===

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(genParticleSrc_,genParticles);

  int ntops=0;
  topDecayCode[0] = 0;
  topDecayCode[1] = 0;

  tauGenPt[0]=0;
  tauGenPt[1]=0;
  tauGenVisPt[0]=0;
  tauGenVisPt[1]=0;

  tauGenNProng[0]=-1;
  tauGenNProng[1]=-1;

  tauGenEta[0]= -99;
  tauGenEta[1]= -99;
  tauGenPhi[0]= -99;
  tauGenPhi[1]= -99;

  for (size_t k = 0 ; k<genParticles->size(); k++) {
    const reco::Candidate & TCand = (*genParticles)[ k ];
    if (abs(TCand.pdgId())== 6 && TCand.status()==3) { //find t quark
      float taupt=-1;
      float tauptvis=-1;
      float taueta=-99;
      float tauphi=-99;
      int nprong=-1;
      int topcode = findTopDecayMode(TCand,taupt,tauptvis,taueta,tauphi,nprong);
      if (ntops<=1) {
	tauGenNProng[ntops]=nprong;
	topDecayCode[ntops] = topcode;
	tauGenPt[ntops]=taupt;
	tauGenVisPt[ntops]=tauptvis;
	tauGenEta[ntops]=taueta;
	tauGenPhi[ntops]=tauphi;
	ntops++;
      }
      else { //should not happen
	std::cout<<"ntops = "<<ntops<<std::endl;
      }
      //std::cout<<topcode<<"\t"<<taupt<<std::endl;

    }
    //also snip Z part for now
  }

  
  ttbarDecayCode =  getTopDecayCategory(topDecayCode[0], topDecayCode[1]);

}


void MakeTauTree::doTauJetMatching(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //run after fillMCInfo has already filled the (global) info on the gen-level taus
  //now we look for reco'd jets that match those taus

  bool foundmatch[2]={false,false};

  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(jetSrc_,jets);
  for (edm::View<pat::Jet>::const_iterator jet = jets->begin(); jet != jets->end(); ++jet) {

    if (jet->pt() < 5) continue; //have no idea what to use here, but jet collection probably already has a cut applied

    double dr1 = reco::deltaR( jet->eta(), jet->phi(), tauGenEta[0],tauGenPhi[0]);
    double dr2 = reco::deltaR( jet->eta(), jet->phi(), tauGenEta[1],tauGenPhi[1]);
    bool match1 = dr1 < matchradius_;
    bool match2 = dr2 < matchradius_;
    //now we know if we have a reco <--> gen match
    int ii;
    if (match1 && match2) { //don't know how often this will happen. take the better match
      ii = (dr1<dr2) ? 0 : 1;
    }
    else if (match1) {
      ii=0;
    }
    else if (match2) {
      ii=1;
    }

    if (match1||match2) {
      //this if will cause us to take only the first match found...(jets should be sorted by pt, so highest pt match)
      if (!foundmatch[ii]) {
	taujetPt[ii] = jet->pt();
	taujetEta[ii] = jet->eta();
	taujetPhi[ii] = jet->phi();
	taujetCSV[ii] = jet->bDiscriminator("combinedSecondaryVertexBJetTags");
	//fabs should not be needed inside of cos() but leave it there anyway
	//MET and METphi should already be filled for this event (ntuple variables), but this is somewhat dangerous I suppose
	taujetMT[ii] = sqrt(2*(jet->pt()*MET)*(1-cos( fabs(reco::deltaPhi(jet->phi(),METphi)) ))); //calculate transverse mass

	//	pat::Jet * thisjet = new pat::Jet(*jet); //just a trick to get the right data format
	taujetChm[ii] = getCHMdr(*jet, 0.3, 1.); //calculate CHM
	taujetLrm[ii] = getLRM(*jet); //LRM
	//delete thisjet; //clean up!

	foundmatch[ii]=true;
      }
    }

  }

}

//this calls steven's redo of the ben hooberman code
void MakeTauTree::isolatedTrackVeto(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<reco::PFCandidateCollection > isolatedcands;
  iEvent.getByLabel("trackIsolationSelector",isolatedcands);

  nIsolatedTracks =  isolatedcands->size();

}

//requires running Ben Hooberman's version of the code in the path
void MakeTauTree::isolatedTrackDetails(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  maxPtIsoTrack_Iso=-1;
  maxPtIsoTrack_Pt=-1;
  minIsoIsoTrack_Iso=1e9;
  minIsoIsoTrack_Pt=-1;

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
    //    std::cout<<ii<<"\t"<<cand_dzpv->at(ii)<<" "<<cand_pt->at(ii)<<" "<<cand_trkiso->at(ii)<<" "<<cand_chg->at(ii)<<std::endl;
    if ( cand_chg->at(ii) != 0  && (fabs(cand_dzpv->at(ii)) < 0.05)) { //charged and compatible with PV
      //      if ( (cand_pt->at(ii) > 10) && (cand_trkiso->at(ii)/cand_pt->at(ii) <0.1) ) nIsolatedTracks++;

      if (cand_pt->at(ii) > 10) { //look for track above pT threshold that is most isolated
	if (cand_trkiso->at(ii)/cand_pt->at(ii) < minIsoIsoTrack_Iso) {
	  minIsoIsoTrack_Iso = cand_trkiso->at(ii)/cand_pt->at(ii);
	  minIsoIsoTrack_Pt = cand_pt->at(ii);
	}
      }

      //now look for track under iso threshold that has highest pT
      if (cand_trkiso->at(ii)/cand_pt->at(ii) <0.10) { 
	if (cand_pt->at(ii) >maxPtIsoTrack_Pt) {
	  maxPtIsoTrack_Iso = cand_trkiso->at(ii)/cand_pt->at(ii);
	  maxPtIsoTrack_Pt = cand_pt->at(ii);
	}
      }

    }
  }

}

void MakeTauTree::POGtaus(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  nLooseTaus15=0;
  nVLooseTaus15=0;
  nMediumTaus15=0;
  nLooseTaus20=0;
  nVLooseTaus20=0;
  nMediumTaus20=0;

  edm::Handle<edm::View<pat::Tau> > taus;
  iEvent.getByLabel(pogtauSrc_,taus);
  //  std::cout<<" == pog taus =="<<std::endl;
  for (edm::View<pat::Tau>::const_iterator tau = taus->begin(); tau != taus->end(); ++tau) {
    //    std::cout<<"\t"<<tau->pt()<<"\t"<<tau->tauID("byLooseCombinedIsolationDeltaBetaCorr")<<" "<<tau->tauID("byVLooseCombinedIsolationDeltaBetaCorr") <<std::endl;
    if ( tau->pt() > 15 ) {
      if ( tau->tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0) nLooseTaus15++;
      if ( tau->tauID("byVLooseCombinedIsolationDeltaBetaCorr") > 0) nVLooseTaus15++;
      if ( tau->tauID("byMediumCombinedIsolationDeltaBetaCorr") > 0) nMediumTaus15++;
      
      if ( tau->pt() > 20 ) {
	if ( tau->tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0) nLooseTaus20++;
	if ( tau->tauID("byVLooseCombinedIsolationDeltaBetaCorr") > 0) nVLooseTaus20++;
	if ( tau->tauID("byMediumCombinedIsolationDeltaBetaCorr") > 0) nMediumTaus20++;
	
      }
    }
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
MakeTauTree::beginJob()
{
  tree_->Branch("nGoodPV",&nGoodPV,"nGoodPV/I");

  tree_->Branch("njets70",&njets70,"njets70/I");
  tree_->Branch("njets50",&njets50,"njets50/I");
  tree_->Branch("njets30",&njets30,"njets30/I");
  tree_->Branch("nCSVMbjets30",&nCSVMbjets30,"nCSVMbjets30/I");
  tree_->Branch("nCSVTbjets30",&nCSVTbjets30,"nCSVTbjets30/I");

  tree_->Branch("MET",&MET,"MET/F");
  tree_->Branch("METphi",&METphi,"METphi/F");

  tree_->Branch("DeltaPhiMetJet1",&DeltaPhiMetJet1,"DeltaPhiMetJet1/F");
  tree_->Branch("DeltaPhiMetJet2",&DeltaPhiMetJet2,"DeltaPhiMetJet2/F");
  tree_->Branch("DeltaPhiMetJet3",&DeltaPhiMetJet3,"DeltaPhiMetJet3/F");

  tree_->Branch("minDeltaPhiMetBM",&minDeltaPhiMetBM,"minDeltaPhiMetBM/F");
  tree_->Branch("minDeltaPhiMetBT",&minDeltaPhiMetBT,"minDeltaPhiMetBT/F");

  tree_->Branch("lowChmJetChm",&lowChmJetChm,"lowChmJetChm/I");
  tree_->Branch("lowChmJetLrm",&lowChmJetLrm,"lowChmJetLrm/F");
  tree_->Branch("lowChmJetCSV",&lowChmJetCSV,"lowChmJetCSV/F");
  tree_->Branch("lowChmJetMT",&lowChmJetMT,"lowChmJetMT/F");
  tree_->Branch("lowChmJetPt",&lowChmJetPt,"lowChmJetPt/F");

  tree_->Branch("nVetoElectrons",&nVetoElectrons,"nVetoElectrons/I");
  tree_->Branch("nVetoMuons",&nVetoMuons,"nVetoMuons/I");
  tree_->Branch("nVetoTaus",&nVetoTaus,"nVetoTaus/I");

  tree_->Branch("nLooseTaus15",&nLooseTaus15,"nLooseTaus15/I");
  tree_->Branch("nVLooseTaus15",&nVLooseTaus15,"nVLooseTaus15/I");
  tree_->Branch("nMediumTaus15",&nMediumTaus15,"nMediumTaus15/I");

  tree_->Branch("nLooseTaus20",&nLooseTaus20,"nLooseTaus20/I");
  tree_->Branch("nVLooseTaus20",&nVLooseTaus20,"nVLooseTaus20/I");
  tree_->Branch("nMediumTaus20",&nMediumTaus20,"nMediumTaus20/I");

  tree_->Branch("nIsolatedTracks",&nIsolatedTracks,"nIsolatedTracks/I");

  tree_->Branch("maxPtIsoTrack_Iso",&maxPtIsoTrack_Iso,"maxPtIsoTrack_Iso/F");
  tree_->Branch("maxPtIsoTrack_Pt",&maxPtIsoTrack_Pt,"maxPtIsoTrack_Pt/F");
  tree_->Branch("minIsoIsoTrack_Iso",&minIsoIsoTrack_Iso,"minIsoIsoTrack_Iso/F");
  tree_->Branch("minIsoIsoTrack_Pt",&minIsoIsoTrack_Pt,"minIsoIsoTrack_Pt/F");


  tree_->Branch("topDecayCode",&topDecayCode,"topDecayCode[2]/I");
  tree_->Branch("ttbarDecayCode",&ttbarDecayCode,"ttbarDecayCode/I");

  tree_->Branch("tauGenPt",&tauGenPt,"tauGenPt[2]/F");
  tree_->Branch("tauGenVisPt",&tauGenVisPt,"tauGenVisPt[2]/F");
  tree_->Branch("tauGenEta",&tauGenEta,"tauGenEta[2]/F");
  tree_->Branch("tauGenPhi",&tauGenPhi,"tauGenPhi[2]/F");
  tree_->Branch("tauGenNProng",&tauGenNProng,"tauGenNProng[2]/I");

  tree_->Branch("taujetPt",&taujetPt,"taujetPt[2]/F");
  tree_->Branch("taujetEta",&taujetEta,"taujetEta[2]/F");
  tree_->Branch("taujetPhi",&taujetPhi,"taujetPhi[2]/F");
  tree_->Branch("taujetCSV",&taujetCSV,"taujetCSV[2]/F");
  tree_->Branch("taujetMT",&taujetMT,"taujetMT[2]/F");
  tree_->Branch("taujetChm",&taujetChm,"taujetChm[2]/I");
  tree_->Branch("taujetLrm",&taujetLrm,"taujetLrm[2]/F");

  //for reasons of code simplicity, this is the *only* place that we turn off part of the code with StoreAllJetInfo
  if (storeAllJetInfo_) {
    tree_->Branch("njets10",&njets10,"njets10/I");
    tree_->Branch("jetPt",jetPt,"jetPt[njets10]/F");
    tree_->Branch("jetEta",jetEta,"jetEta[njets10]/F");
    tree_->Branch("jetPhi",jetPhi,"jetPhi[njets10]/F");
    tree_->Branch("jetCSV",jetCSV,"jetCSV[njets10]/F");
    tree_->Branch("jetMT",jetMT,"jetMT[njets10]/F");
    tree_->Branch("jetChm",jetChm,"jetChm[njets10]/I");
    tree_->Branch("jetLrm",jetLrm,"jetLrm[njets10]/F");
    tree_->Branch("jetPassTauId",jetPassTauId,"jetPassTauId[njets10]/O");
    tree_->Branch("jetMatchGenTau",jetMatchGenTau,"jetMatchGenTau[njets10]/O");
    tree_->Branch("jetTauDr",jetTauDr,"jetTauDr[njets10]/F");//to study matching
  }

}

// ------------ method called once each job just after ending the event loop  ------------
void 
MakeTauTree::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MakeTauTree::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
MakeTauTree::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MakeTauTree::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MakeTauTree::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MakeTauTree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MakeTauTree);
