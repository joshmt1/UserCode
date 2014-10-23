#include "FindHemispheres.hh"

#include <iostream>
#include <cassert>
#include "classes/DelphesClasses.h"

using namespace std;


FindHemispheres::FindHemispheres( TClonesArray * jets) : 
  maxIter_(100),
  jets_(jets),
  done_(false)
{

}

FindHemispheres::~FindHemispheres() { }


void FindHemispheres::find() {
  if (done_) return;
  const bool debug=false;

  /*
    
  find initial axes: directions of two jets with largest invariant mass, assuming jets are massless
  
  <set>s: jets in hemi 1, jets in hemi 2
  
  while(true) {
  for each of the other jets
  find Lund distance for each hemisphere
  associate jet with hemisphere with smaller Lund distance (do NOT touch the axes)
  
  recompute axes using sum of momenta of jets in the hemisphere
  
  see if jets in sets changed compared to previous iteration; if no, break
  if yes, update the <set>s
  }
  
  
  */
  
  //find initial axes
  float maxDijetInvMass = 0;
  int seed1=-1,seed2=-1;
  for (int ijet = 0 ; ijet<jets_->GetEntries() ; ++ijet) {
    if ( !isGoodJet(ijet) ) continue;
    for (int jjet = 0 ; jjet<jets_->GetEntries() ; ++jjet) {
      if (ijet==jjet) continue;
      if ( !isGoodJet(jjet) ) continue;

      float mjj=getDijetInvMass( ijet,jjet);
      if ( mjj>maxDijetInvMass) {
	maxDijetInvMass = mjj;
	seed1 = ijet;
	seed2 = jjet;
      }
    }
  }


  if (debug)  cout<<"Got initial axes jets = "<<seed1<<" "<<seed2<<endl;
  if (seed1==-1 || seed2==-1) return;

  set<int> hemi1,hemi2;
  hemi1.insert(seed1);
  hemi2.insert(seed2);
  int counter = 0;
  while (true) {
    if (debug) cout<<"it "<<counter<<endl;

    TLorentzVector axis1 = getAxis( hemi1);
    TLorentzVector axis2 = getAxis( hemi2);

    set<int> newhemi1,newhemi2;
    //    for (int ijet = 0 ; ijet<jets_->GetEntries() ; ++ijet) {
    for (int ijet = jets_->GetEntries()-1 ; ijet>=0 ; ijet--) {
      //on first iteration only, skip the seed jets
      if ( ijet==seed1) {newhemi1.insert(ijet); continue;}
      if ( ijet==seed2) {newhemi2.insert(ijet); continue;}

      if ( !isGoodJet(ijet) ) continue;
      
      double lund1 = getLundDistance( axis1, ijet);
      double lund2 = getLundDistance( axis2, ijet);
      if (debug)  cout<<"lund distances = "<<lund1<<" "<<lund2<<endl;

      //problem to be understood or solved -- sometimes one of the hemispheres is empty.
      //need to reject this as a possibility!
      // so look for the case where we're about to end the loop with an empty set and prevent it
      //      if ( (ijet+1 == jets_->GetEntries()) && (newhemi1.size()==0 || newhemi2.size()==0)) {
      if ( (ijet == 0) && (newhemi1.size()==0 || newhemi2.size()==0)) {
	if ( newhemi1.size() == 0) newhemi1.insert(ijet);
	else	if ( newhemi2.size() == 0) newhemi2.insert(ijet);
      }
      else { //otherwise use the default method
	if (lund1<lund2) newhemi1.insert(ijet);
	else newhemi2.insert(ijet);
      }
    }


    if (debug) {
      for ( set<int>::iterator ij=newhemi1.begin() ; ij!=newhemi1.end() ; ++ij) cout<<*ij<<" ";
      cout<<endl;
      for ( set<int>::iterator ij=newhemi2.begin() ; ij!=newhemi2.end() ; ++ij) cout<<*ij<<" ";
      cout<<endl;
    }

    if ( (newhemi1==hemi1 && newhemi2==hemi2) || (newhemi1==hemi2 && newhemi2==hemi1)) break; //done!
    if (++counter > maxIter_) {/*cout<<"Hit iteration max!"<<endl;*/ break;}
    hemi1 = newhemi1;
    hemi2 = newhemi2;
    seed1=-1; seed2=-1; // disable the 'first iteration' test at the beginning of the loop
  }

  hemi1_ = hemi1;
  hemi2_ = hemi2;
  //redundant to the calculation inside the loop, but whatever
  pjet1_ = getAxis(hemi1_);
  pjet2_ = getAxis(hemi2_);

  done_=true;
}

TLorentzVector FindHemispheres::getAxis( const set<int> & jets ) {
  //sum the momenta in jets to get an axis

  TLorentzVector ax;
  for ( set<int>::iterator jet = jets.begin(); jet!=jets.end(); ++jet) {
    TLorentzVector jet4v = ((Jet*) jets_->At( *jet))->P4();
    ax+= jet4v;
  }
  return ax;

}

double FindHemispheres::getLundDistance( const TLorentzVector & ax, int ijet) {
  //calculate lund distance between axis and ijet

  TLorentzVector jet = ((Jet*) jets_->At(ijet))->P4();
  return  ((ax.E() - ax.P()*cos( ax.Angle(jet.Vect()) ))*ax.E()) / pow( ax.E() + jet.E(),2);

}

double FindHemispheres::getDijetInvMass( int i, int j) {
  //AN says "two (massless) jets which have the largest invariant mass"
  //not sure how literally I should take that massless part. Ignore it for now.

  TLorentzVector j1 = ((Jet*) jets_->At(i))->P4();
  TLorentzVector j2 = ((Jet*) jets_->At(j))->P4();

  TLorentzVector jj = j1+j2;
  return jj.M();

}

double FindHemispheres::Px(int hemi) {
  assert(hemi == 1 || hemi==2);
  if (!done_) find();
  if (hemi1_.size()==0 || hemi2_.size()==0) return -1;
  if (hemi==1) return pjet1_.Px();
  return pjet2_.Px();
}
double FindHemispheres::Py(int hemi) {
  assert(hemi == 1 || hemi==2);
  if (!done_) find();
  if (hemi1_.size()==0 || hemi2_.size()==0) return -1;
  if (hemi==1) return pjet1_.Py();
  return pjet2_.Py();
}

bool FindHemispheres::isGoodJet(int jj) {

/*
jet quality cuts: 20 GeV, not sure about eta, etc
*/

  return  ( (((Jet*) jets_->At(jj))->PT)>20);
}

vector<int> FindHemispheres::getGrouping() {
  if (!done_) find();

  vector<int> g;

  //  if (hemi1_.size()==0 || hemi2_.size()==0) return -1;
  for (int ijet = 0 ; ijet<jets_->GetEntries() ; ++ijet) {
   
    if ( hemi1_.count( ijet)) g.push_back( 1);
    else    if ( hemi2_.count( ijet)) g.push_back( 2);
    else g.push_back(0); //should not happen, i think 

  }
  return g;
}
