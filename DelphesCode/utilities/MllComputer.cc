#include "MllComputer.hh"

#include <iostream>
#include <cassert>
#include "classes/DelphesClasses.h"
#include "external/mt2analysis/Utilities.hh"

using namespace std;

MllComputer::MllComputer(TClonesArray* el, TClonesArray* mu) :
  el_(el),
  mu_(mu),
  foundGood_(false),
  maxEta_(-99)
{
}

MllComputer::~MllComputer() {
}

float MllComputer::GetMll() {

  //get highest pT *good* lepton, where good means in eta acceptance
  //(in the real analysis, good also implies some quality cuts)
  if (!foundGood_) findGoodLeptons();
  if (lep_pt_.size() <2) return -2;

  //
  set< pair<float, pair<LeptonFlavor,int> > >::reverse_iterator rit;
  int ii=0;
  int    lead_lepton_charge=0;
  TLorentzVector l1;
  for (rit=lep_pt_.rbegin(); rit != lep_pt_.rend(); ++rit) {

    if (ii==0) { //find charge of lead lepton
      lead_lepton_charge =  GetCharge(rit);
      l1 = Get4Vector(rit);
    }
    else if (ii>0) { //for non-leading leptons
      //find 2nd highest pT lepton with opposite charge and DR>0.3 compared to first lepton
      int ch2 = GetCharge(rit); 
      //      cout<<"\tcharges "<<lead_lepton_charge<<" "<<ch2<<endl;
      if (lead_lepton_charge != ch2) { //OS lepton pair
	//compute DR
	TLorentzVector l2 = Get4Vector(rit);
	float dr=Util::GetDeltaR(l2.Eta(),l1.Eta(),l2.Phi(),l1.Phi());
	//	cout<<"\t\tDR="<<dr<<endl;
	if (dr>0.3) { //passes dr cut
	  TLorentzVector ll = l1+l2;
	  //	  cout<<"Mll = "<<ll.M()<<endl;
	  if ( std::abs(l2.Eta()) > std::abs(l1.Eta()) ) maxEta_ = std::abs(l2.Eta());
	  else maxEta_ = std::abs(l1.Eta());
	  return ll.M();
	}
      }
    }

    ++ii;
  }
 


  return -1;
}

void MllComputer::findGoodLeptons() {

  //find good leptons and also make a joint list of leptons (e/mu combined) sorted by pt

  //also reject pairs with either lepton in eta of 1.4-1.6 

  for (int i = 0 ; i < el_->GetEntries() ; i++) {
    Electron *el = (Electron*) el_->At(i);
    if (el->PT < 20 ) continue;
    float abseta = std::abs(el->Eta);
    if ( abseta > 2.4) continue; 
    if ( abseta>1.4 && abseta<1.6) continue;//remove eta range of 1.4-1.6
    //TODO remove barrel/endcap transition
    lep_pt_.insert( make_pair( el->PT, make_pair(kElectron,i))); //set sorts by pt, low to high
  }

  for (int i = 0 ; i < mu_->GetEntries() ; i++) {
    Muon *mu = (Muon*) mu_->At(i);
    if (mu->PT < 20 ) continue;
    float abseta = std::abs(mu->Eta);
    if ( abseta > 2.4) continue; 
     if ( abseta>1.4 && abseta<1.6) continue;//remove eta range of 1.4-1.6
    lep_pt_.insert( make_pair( mu->PT, make_pair(kMuon,i))); //set sorts by pt, low to high
  }
 

  foundGood_=true;
}

int MllComputer::GetCharge( set< pair<float, pair<LeptonFlavor,int> > >::reverse_iterator & rit) {
  int ch=0;

  int lepton_index = rit->second.second;
  LeptonFlavor lepton_flavor = rit->second.first;
  if ( lepton_flavor==kElectron)
    ch=   ((Electron*) el_->At(lepton_index))->Charge;
  else if ( lepton_flavor==kMuon)
    ch=   ((Muon*) mu_->At(lepton_index))->Charge;
  else assert(0);

  return ch;
}

TLorentzVector MllComputer::Get4Vector( set< pair<float, pair<LeptonFlavor,int> > >::reverse_iterator & rit) {
  TLorentzVector v;

  int lepton_index = rit->second.second;
  LeptonFlavor lepton_flavor = rit->second.first;
  if ( lepton_flavor==kElectron)
    v= ((Electron*) el_->At(lepton_index))->P4();
  else if ( lepton_flavor==kMuon)
    v=  ((Muon*) mu_->At(lepton_index))->P4();
  else assert(0);

  return v;

}

