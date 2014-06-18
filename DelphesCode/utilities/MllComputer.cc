#include "MllComputer.hh"

#include <iostream>
#include <cassert>
#include "classes/DelphesClasses.h"
#include "external/mt2analysis/Utilities.hh"


using namespace std;

MllComputer::MllComputer(TClonesArray* el, TClonesArray* mu) :
  minpt_(20),
  maxetacut_(2.4),
  removegap_(true),
  maxreliso_(0.15),
  randomizeLeptons_(false),
  seed_(1234),
  el_(el),
  mu_(mu),
  rand_(0),
  foundGood_(false),
  maxEta_(-99),
  isSF_(false),
  p4_l1_(0),
  p4_l2_(0),
  l1_iso_(-99),
  l2_iso_(-99),
  l1_flavor_(-99),
  l2_flavor_(-99)
{
}

MllComputer::~MllComputer() {
 delete p4_l1_;
 delete p4_l2_;
 delete rand_;
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

float MllComputer::GetMll_random() {

  if (rand_==0) rand_ =  new TRandom3(seed_); //do *NOT* always use same seed

  //get Mll *not* selecting leptons by pt sorting
  //too bad that i need to duplicate some code, but i think it is much cleaner this way
  if (!foundGood_) findGoodLeptons();
  if (lep_pt_.size() <2) return -2;

  if (p4_l1_!=0) delete p4_l1_;
  if (p4_l2_!=0) delete p4_l2_;

  //convert set to vector
  vector< std::pair<float, std::pair<LeptonFlavor,int> > > lep_pt;

  set< pair<float, pair<LeptonFlavor,int> > >::iterator it;
  for (it=lep_pt_.begin(); it != lep_pt_.end(); ++it) {
    lep_pt.push_back( *it);
  }

  //pick out the two random leptons to use
  set<int> used;
  unsigned int n = lep_pt.size();
  int i1 = rand_->Integer(n);
  used.insert(i1);
  if ( lep_pt[i1].second.first ==kElectron) {
    l1_flavor_ = 11 * ((Electron*)el_->At(lep_pt[i1].second.second))->Charge;
    l1_iso_ =  ((Electron*)el_->At(lep_pt[i1].second.second))->IsolationVar;
    p4_l1_ = new TLorentzVector( ((Electron*)el_->At(lep_pt[i1].second.second))->P4());
  }
  else if (lep_pt[i1].second.first==kMuon) {
    l1_flavor_ = 13 * ((Muon*)mu_->At(lep_pt[i1].second.second))->Charge;
    l1_iso_ =  ((Muon*)mu_->At(lep_pt[i1].second.second))->IsolationVar;
    p4_l1_ = new TLorentzVector( ((Muon*)mu_->At(lep_pt[i1].second.second))->P4());
  }
  int i2=-1;
  while (used.size()<n) {

    i2 = rand_->Integer(n);
    if (i1==i2) {i2=-1; continue;}

    if ( lep_pt[i2].second.first ==kElectron) {
      l2_flavor_ = 11 * ((Electron*)el_->At(lep_pt[i2].second.second))->Charge;
      l2_iso_ = ((Electron*)el_->At(lep_pt[i2].second.second))->IsolationVar;
      p4_l2_ = new TLorentzVector( ((Electron*)el_->At(lep_pt[i2].second.second))->P4()); 
    }
    else if (lep_pt[i2].second.first==kMuon) {
      l2_flavor_ = 13 * ((Muon*)mu_->At(lep_pt[i2].second.second))->Charge;
      l2_iso_ = ((Muon*)mu_->At(lep_pt[i2].second.second))->IsolationVar;
      p4_l2_ = new TLorentzVector( ((Muon*)mu_->At(lep_pt[i2].second.second))->P4());
    }
    //ensure opposite charge!
    if ( l1_flavor_/std::abs(l1_flavor_) == l2_flavor_/std::abs(l2_flavor_) ) {
      used.insert(i2);
      i2=-1;
      continue;
    }

    //ok, i2 is an OS lepton
    break;
  }

  if (i2<0 )  {
    //clear other stuff?
    return -1;
  }


  isSF_ = (std::abs(l1_flavor_)==std::abs(l2_flavor_));
  maxEta_ = (std::abs(p4_l1_->Eta() )>std::abs(p4_l2_->Eta())) ? std::abs(p4_l1_->Eta() ):std::abs(p4_l2_->Eta());

  //boy, this code was really hard-coded like mad with set in mind
  //need:
  /*
  float maxEta_;//largest value found in event
  bool isSF_;
  TLorentzVector* p4_l1_;
  TLorentzVector* p4_l2_;
  int l1_flavor_;
  int l2_flavor_;
  */

  TLorentzVector mll=(*p4_l1_)+(*p4_l2_);
  return mll.M();

}

float MllComputer::GetMll() {

  if (randomizeLeptons_) return GetMll_random();

  //get highest pT *good* lepton, where good means in eta acceptance
  //(in the real analysis, good also implies some quality cuts)
  if (!foundGood_) findGoodLeptons();
  if (lep_pt_.size() <2) return -2;

  if (p4_l1_!=0) delete p4_l1_;
  if (p4_l2_!=0) delete p4_l2_;

  float mll = -1;
  extra_leptons_.clear();
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
      l1_iso_ = GetIsolation(rit);
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
	if (dr>0.3 && mll<0) { //passes dr cut
	  TLorentzVector ll = (*p4_l1_)+l2;
	  //	  cout<<"Mll = "<<ll.M()<<endl;
	  if ( std::abs(l2.Eta()) > std::abs(p4_l1_->Eta()) ) maxEta_ = std::abs(l2.Eta());
	  else maxEta_ = std::abs(p4_l1_->Eta());
	  LeptonFlavor l2fl = rit->second.first;
	  isSF_ = (l1fl == l2fl);
	  p4_l2_ = new TLorentzVector(l2);
	  l2_flavor_ = (l2fl==kElectron) ? 11 : 13;
	  l2_flavor_ *= ch2;
	  l2_iso_ = GetIsolation(rit);
	  mll = ll.M();
	}
	else {	//lepton that won't go into mLL: add to extra leptons
	  extra_leptons_.insert(*rit);
	}
      }
      else {	//lepton that won't go into mLL: add to extra leptons
	extra_leptons_.insert(*rit);
      }
    }

    ++ii;
  }
 


  return mll;
}

TLorentzVector* MllComputer::GetLeptonP4(int index) {
  if (index==1) return p4_l1_;
  if (index==2) return p4_l2_;
  return 0;
}


TLorentzVector MllComputer::GetExtraLeptonP4(int index) {
  set< pair<float, pair<LeptonFlavor,int> > >::reverse_iterator rit;
  int i=0;
  for (rit=extra_leptons_.rbegin(); rit != extra_leptons_.rend(); ++rit) {
    if (i==index) {
      return Get4Vector(rit);
    }
    ++i;
  }

  TLorentzVector d;
  return d;
}

int MllComputer::GetLeptonFlavor(int index) {
  if (index==1) return l1_flavor_;
  if (index==2) return l2_flavor_;
  return 0;
}

float MllComputer::GetLeptonIsolation(int index) {
  if (index==1) return l1_iso_;
  if (index==2) return l2_iso_;
  return 0;
}

void MllComputer::findGoodLeptons() {

  //find good leptons and also make a joint list of leptons (e/mu combined) sorted by pt

  //also reject pairs with either lepton in eta of 1.4-1.6 

  for (int i = 0 ; i < el_->GetEntries() ; i++) {
    Electron *el = (Electron*) el_->At(i);
    if (el->PT < minpt_ ) continue;
    if (el->IsolationVar >maxreliso_) continue;
    float abseta = std::abs(el->Eta);
    if ( abseta > maxetacut_) continue; 
    if ( removegap_ && abseta>1.4 && abseta<1.6) continue;//remove eta range of 1.4-1.6
    //TODO remove barrel/endcap transition
    lep_pt_.insert( make_pair( el->PT, make_pair(kElectron,i))); //set sorts by pt, low to high
  }

  for (int i = 0 ; i < mu_->GetEntries() ; i++) {
    Muon *mu = (Muon*) mu_->At(i);
    if (mu->PT < minpt_ ) continue;
    if (mu->IsolationVar >maxreliso_) continue;
    float abseta = std::abs(mu->Eta);
    if ( abseta > maxetacut_) continue; 
     if ( removegap_ && abseta>1.4 && abseta<1.6) continue;//remove eta range of 1.4-1.6
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

float MllComputer::GetIsolation( set< pair<float, pair<LeptonFlavor,int> > >::reverse_iterator & rit) {
  float val=0;

  int lepton_index = rit->second.second;
  LeptonFlavor lepton_flavor = rit->second.first;
  if ( lepton_flavor==kElectron)
    val=   ((Electron*) el_->At(lepton_index))->IsolationVar;
  else if ( lepton_flavor==kMuon)
    val=   ((Muon*) mu_->At(lepton_index))->IsolationVar;
  else assert(0);

  return val;
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

