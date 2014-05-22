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
  maxEta_(-99),
  isSF_(false),
  p4_l1_(0),
  p4_l2_(0),
  l1_flavor_(-99),
  l2_flavor_(-99)
{
}

MllComputer::~MllComputer() {
 delete p4_l1_;
 delete p4_l2_;
}

float MllComputer::GetMee_Test() {
  //not for serious use
  //just trying to match a plot made by somebody else

  /*
  int charge1=0;
  TLorentzVector e1;
  for (int i = 0 ; i < el_->GetEntries() ; i++) {
    Electron *el = (Electron*) el_->At(i);
    if (charge1==0) {
      charge1 = el->Charge;
      e1=el->P4();
    }
    else { //we already found the first electron
      if ( el->Charge == -charge1) {  //found OS leptons
	TLorentzVector e2 = el->P4();
	TLorentzVector ll = e1+e2;
	return ll.M();
      }
      //if same sign, ignore and move on
    }
  }
  */

  //insert batool's code directly
  Electron *elec1, *elec2;
  // If event contains at least 2 electrons
  if(el_->GetEntries() > 1)
    {
      // Take first two electrons
      elec1 = (Electron *) el_->At(0);
      elec2 = (Electron *) el_->At(1);

      //  their invariant mass
      if ((elec1->Charge+elec2->Charge) !=0) return -10;
      return ((elec1->P4()) + (elec2->P4())).M();
    }

  return 0;
}

float MllComputer::GetMll() {

  //get highest pT *good* lepton, where good means in eta acceptance
  //(in the real analysis, good also implies some quality cuts)
  if (!foundGood_) findGoodLeptons();
  if (lep_pt_.size() <2) return -2;

  if (p4_l1_!=0) delete p4_l1_;
  if (p4_l2_!=0) delete p4_l2_;

  //
  set< pair<float, pair<LeptonFlavor,int> > >::reverse_iterator rit;
  int ii=0;
  int    lead_lepton_charge=0;
  LeptonFlavor l1fl=kMuon;
  //  TLorentzVector l1;
  for (rit=lep_pt_.rbegin(); rit != lep_pt_.rend(); ++rit) {

    if (ii==0) { // lead lepton info
      lead_lepton_charge =  GetCharge(rit);
      p4_l1_ = new TLorentzVector( Get4Vector(rit) );
      l1fl = rit->second.first;
      l1_flavor_ = (l1fl==kElectron) ? 11 : 13;
      l1_flavor_ *= lead_lepton_charge;
    }
    else if (ii>0) { //for non-leading leptons
      //find 2nd highest pT lepton with opposite charge and DR>0.3 compared to first lepton
      int ch2 = GetCharge(rit); 
      //      cout<<"\tcharges "<<lead_lepton_charge<<" "<<ch2<<endl;
      if (lead_lepton_charge != ch2) { //OS lepton pair
	//compute DR
	TLorentzVector l2 = Get4Vector(rit);
	float dr=Util::GetDeltaR(l2.Eta(),p4_l1_->Eta(),l2.Phi(),p4_l1_->Phi());
	//	cout<<"\t\tDR="<<dr<<endl;
	if (dr>0.3) { //passes dr cut
	  TLorentzVector ll = (*p4_l1_)+l2;
	  //	  cout<<"Mll = "<<ll.M()<<endl;
	  if ( std::abs(l2.Eta()) > std::abs(p4_l1_->Eta()) ) maxEta_ = std::abs(l2.Eta());
	  else maxEta_ = std::abs(p4_l1_->Eta());
	  LeptonFlavor l2fl = rit->second.first;
	  isSF_ = (l1fl == l2fl);
	  p4_l2_ = new TLorentzVector(l2);
	  l2_flavor_ = (l2fl==kElectron) ? 11 : 13;
	  l2_flavor_ *= ch2;
	  return ll.M();
	}
      }
    }

    ++ii;
  }
 


  return -1;
}

TLorentzVector* MllComputer::GetLeptonP4(int index) {
  if (index==1) return p4_l1_;
  if (index==2) return p4_l2_;
  return 0;
}

int MllComputer::GetLeptonFlavor(int index) {
  if (index==1) return l1_flavor_;
  if (index==2) return l2_flavor_;
  return 0;
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

