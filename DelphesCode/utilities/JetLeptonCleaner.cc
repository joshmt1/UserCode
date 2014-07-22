#include "JetLeptonCleaner.hh"

#include <iostream>
#include <cassert>

#include "classes/DelphesClasses.h"
#include "external/mt2analysis/Utilities.hh"

using namespace std;

#include <vector>

JetLeptonCleaner::JetLeptonCleaner(TClonesArray* jets,TClonesArray* els,TClonesArray* mus,TClonesArray* photons) :
  isolationCut_(0.2),//was 0.2
  drmin_(0.3) //or 0.5?
{


  //  cout<<" ----- "<<endl;
  //  cout<<jets->GetEntries()<<" " <<els->GetEntries()<<" "<<mus->GetEntries()<<" "      <<photons->GetEntries()<<" "<<endl;

  applyIsolation<Electron>(els);
  applyIsolation<Muon>(mus);
  if (photons!=0) applyIsolation<Photon>(photons);

  //  cout<<jets->GetEntries()<<" " <<els->GetEntries()<<" "<<mus->GetEntries()<<" "      <<photons->GetEntries()<<" "<<endl;

  cleanJets<Electron>(jets,els);
  cleanJets<Muon>(jets,mus);
  if (photons!=0) cleanJets<Photon>(jets,photons);

  //  cout<<jets->GetEntries()<<" " <<els->GetEntries()<<" "<<mus->GetEntries()<<" "      <<photons->GetEntries()<<" "<<endl;
  //  cout<<" ----- "<<endl;

}

JetLeptonCleaner::~JetLeptonCleaner() {

}

template <class T> void JetLeptonCleaner::cleanJets(TClonesArray* jets,TClonesArray* leps) {

  vector<Jet*> good;
  
  const  int nj= jets->GetEntries();
  
  for (int ijet = 0 ; ijet < nj ; ++ijet) {
    Jet * jet = (Jet*) jets->At(ijet);
    //see if jet is near a lepton/photon
    bool isgood=true;
    for (int ilep = 0 ; ilep < leps->GetEntries() ; ++ilep) {
      T * lep = (T*) leps->At(ilep);
      if ( Util::GetDeltaR(jet->Eta,lep->Eta,jet->Phi,lep->Phi) < drmin_) {
	isgood=false;
	//	cout<<" jet removed! dr = "<< Util::GetDeltaR(jet->Eta,lep->Eta,jet->Phi,lep->Phi)<<endl;
      }
    }
    if (isgood) good.push_back(jet);
  }
	 
  if (nj != (int)good.size()) { //if nothing has changed, no need to clear and refill
    jets->Delete();//Clear() leads to huge memory leaks!
    for (size_t k=0;k<good.size();k++) {
      new( (*jets)[k]) Jet( *good.at(k)); //see http://root.cern.ch/root/html/TClonesArray.html
    }
  }

}

template <class T> void JetLeptonCleaner::applyIsolation(TClonesArray* collection) {

  vector<T*> good;

  int n=collection->GetEntries();
  for (int i = 0 ; i < n ; ++i) {
    T * e = (T*) collection->At(i);
    if (e->IsolationVar <= isolationCut_)       good.push_back(e);
  }

  if (n != (int)good.size()) { //if nothing has changed, no need to clear and refill
    collection->Delete(); //Clear() leads to huge memory leaks!
    for (size_t k=0;k<good.size();k++) {
      new( (*collection)[k]) T( *good.at(k)); //see http://root.cern.ch/root/html/TClonesArray.html
    }
  }

}

