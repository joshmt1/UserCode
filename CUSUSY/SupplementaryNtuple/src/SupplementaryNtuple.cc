#include <iostream>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "TTree.h"

#include "CUSUSY/SupplementaryNtuple/interface/SupplementaryNtuple.h"

SupplementaryNtuple::SupplementaryNtuple(const edm::ParameterSet& iConfig):
  genParticles_(iConfig.getParameter<edm::InputTag>("genParticles")),
  susytree_(0)
{

}


SupplementaryNtuple::~SupplementaryNtuple()
{

}
// ------------ method called once each job just before starting event loop  ------------
void 
SupplementaryNtuple::beginJob()
{

  susytree_=new TTree("susytree","susytree");
  susytree_->Branch("SUSY_nb", &SUSY_nb, "SUSY_nb/I");

}


// ------------ method called once each job just after ending the event loop  ------------
void 
SupplementaryNtuple::endJob() { 
  
  std::cout << "+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=" << std::endl;
//   std::cout << "\n\nSummary:"<< std::endl;
//   //  std::cout << "number of events processed: " << histocontainer_["eventcount"]->GetEntries() << std::endl;
  std::cout << "number of events added to SUSY tree: " << susytree_->GetEntries() << std::endl;
  std::cout << "\n+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=" << std::endl;

}

// bool
// SupplementaryNtuple::checkMother( int mother ) {


// }


//return the number of SUSY mothers
int
SupplementaryNtuple::findSUSYMaternity( const reco::Candidate & cand ) {

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

// ------------ method called for each event  ------------
void
SupplementaryNtuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // ==== MC truth information ====

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(genParticles_,genParticles);

  //  std::cout<<" -- event MC truth tree -- "<<std::endl; //for debugging

  //variables in the tree
  SUSY_nb=0;

  for (size_t k = 0 ; k<genParticles->size(); k++) {

    const reco::Candidate & TCand = (*genParticles)[ k ];

//     //for debugging
//     int motherid=-1;
//     if ( TCand.numberOfMothers() >0) { motherid = TCand.mother()->pdgId();}
//     std::cout<< TCand.status()<<"\t"<<TCand.pdgId()<<"\t"<<motherid<<std::endl;;

    if (abs(TCand.pdgId())== 5 && TCand.status()==3) { //find b quark
      int nSUSY=findSUSYMaternity(TCand);
      //if we have one b, and it has nSUSY greater than 1, we still want to just count the b once
      if (  nSUSY>0 ) SUSY_nb++;
    }
  }
  susytree_->Fill();
  
}


//define this as a plug-in
DEFINE_FWK_MODULE(SupplementaryNtuple);
