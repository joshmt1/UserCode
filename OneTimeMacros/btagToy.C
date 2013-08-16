#include <vector>
#include <iostream>
#include <utility>

#include "TRandom3.h"

// toy for b-tag probability formula testing
enum tagStatus {k0=0, kL=1, kM=2,kT=3,kMax=4};
enum jetFlavor {klf = 1, kc = 4, kb=5}; //change these?

TRandom3 * rnd_=0;

using namespace std;

//a 'jet' would be a (flavor, tag status) pair
//an 'event' would be a vector of jets

//a 'dataset' is a vector of 'events'

tagStatus nextTagStatus(tagStatus ts) {

  if       (ts == k0 ) return kL;
  else  if (ts == kL ) return kM;
  else  if (ts == kM ) return kT;
  else  if (ts == kT ) return kMax;
  else assert(0);

  return k0;

}

double getEffMC( const jetFlavor flav, const tagStatus tag) {

  if (tag==kMax) return 0;

  double eff=1; //this is for the k0 (untagged) case

  if (flav==klf) {
    if      (tag==kL) eff = 0.1;
    else if (tag==kM) eff = 0.01;
    else if (tag==kT) eff = 0.001;
  }
  else if (flav==kc) {
    //i have no idea if these are right
    if      (tag==kL) eff = 0.3;
    else if (tag==kM) eff = 0.2;
    else if (tag==kT) eff = 0.1;
  }
  else if (flav==kb) {
    //again, just a guess
    if      (tag==kL) eff = 0.8;
    else if (tag==kM) eff = 0.7;
    else if (tag==kT) eff = 0.6;
  }
  return eff;

}

bool closureTestMode_=false;
double getSF(const jetFlavor flav, const tagStatus tag) {

  //for special test of Pawandeep's code
  if (closureTestMode_)  return 1;

  //normal code

  if (tag==k0) return 1;
  if (tag==kMax) return 0;
  double sf = 1;


  if (flav==klf) {
    if      (tag==kL) sf = 1.1;
    else if (tag==kM) sf = 1.3;
    else if (tag==kT) sf = 1.5;
  }
  else if (flav==kc) {
    if      (tag==kL) sf = 0.95;
    else if (tag==kM) sf = 0.9;
    else if (tag==kT) sf = 0.8;
  }
  else if (flav==kb) {
    if      (tag==kL) sf = 0.95;
    else if (tag==kM) sf = 0.9;
    else if (tag==kT) sf = 0.8;
  }
  return sf;

}

tagStatus getTagStatus(const bool isMC, const jetFlavor flav ) {
  //here we figure out if the jet is tagged or not

  double eff_L=getEffMC(flav,kL);
  double eff_M=getEffMC(flav,kM);
  double eff_T=getEffMC(flav,kT);

  if (!isMC) {
    eff_L *= getSF(flav,kL);
    eff_M *= getSF(flav,kM);
    eff_T *= getSF(flav,kT);
  }

  assert(eff_M>eff_T);
  assert(eff_L>eff_M);

  double r= rnd_->Uniform() ;

  tagStatus tag=k0;

  if (r <= eff_L) tag = kL;
  if (r <= eff_M) tag = kM;
  if (r <= eff_T) tag = kT;

  return tag;
}

float chanceOfJet5_ = 0.50;
vector<pair<jetFlavor,tagStatus> > generateTtbar(const bool isMC) {

  vector<pair<jetFlavor,tagStatus> > ttbarEvent;

  //generate 2 b, and 50% c+lf/lf+lf
  jetFlavor b1 = kb;
  jetFlavor b2 = kb;

  jetFlavor Wj1,Wj2;
  if (rnd_->Uniform() >0.5) { //50% chance of charm
    Wj1 = kc;
  }
  else {
    Wj1 = klf;
  }
  Wj2 = klf; //2nd W jet always light

  //now, need to get the tag status for these jets
  ttbarEvent.push_back(make_pair( b1, getTagStatus(isMC, b1)));
  ttbarEvent.push_back(make_pair( b2, getTagStatus(isMC, b2)));

  ttbarEvent.push_back(make_pair( Wj1, getTagStatus(isMC, Wj1)));
  ttbarEvent.push_back(make_pair( Wj2, getTagStatus(isMC, Wj2)));

  //see if there's a 5th jet
  if (rnd_->Uniform() <chanceOfJet5_) {
    jetFlavor extraJet;
    double r = rnd_->Uniform();
    //these numbers are not very important; they determine the true flavor composition of the extra jet
    if (r < 0.4) extraJet = klf;
    else if (r < 0.7) extraJet = kc;
    else extraJet = kb;
    ttbarEvent.push_back(make_pair( extraJet, getTagStatus(isMC, extraJet)));
  }

  return ttbarEvent;

}

double getEventWeight_LMT(const vector<pair<jetFlavor,tagStatus> > & event) {

  double prodMC=1;
  double prodD =1;
  //loop over jets
  for (size_t ijet = 0; ijet<event.size(); ++ijet) {

    pair<jetFlavor,tagStatus> thisjet = event.at(ijet);

    double eff_J      = getEffMC( thisjet.first,thisjet.second);
    double eff_Jplus1 = getEffMC( thisjet.first, nextTagStatus(thisjet.second));
    double SF_J      = getSF( thisjet.first,thisjet.second);
    double SF_Jplus1 = getSF(thisjet.first,nextTagStatus(thisjet.second));
    double StimesEff_J = SF_J*eff_J;
    double StimesEff_Jp1 = SF_Jplus1*eff_Jplus1;
    prodMC *= (eff_J - eff_Jplus1 );
    prodD*= (StimesEff_J - StimesEff_Jp1);

  }

  return prodD/prodMC;

}


double getEventWeight_2b(const vector<pair<jetFlavor,tagStatus> > & event) {

  double prodMC=1;
  double prodD =1;
  //loop over jets
  for (size_t ijet = 0; ijet<event.size(); ++ijet) {

    pair<jetFlavor,tagStatus> thisjet = event.at(ijet);

    if (thisjet.second == kT) {
      prodMC *= getEffMC( thisjet.first,thisjet.second);
      prodD  *= (getSF(thisjet.first,thisjet.second)*getEffMC( thisjet.first,thisjet.second));
    }
    else if (thisjet.second == kL || thisjet.second == k0) { //not tagged as medium
      prodMC *= ( 1.0- getEffMC( thisjet.first,kM));
      prodD *= ( 1.0- getSF(thisjet.first,kM)*getEffMC( thisjet.first,kM));
    }
    else assert(0); //nothing else should be possible!

  }

  return prodD/prodMC;

}


void countInBins( const vector< vector<pair<jetFlavor,tagStatus> >  > & dataset,double & n2b,double & n3b,double & n4b,const bool reweight=false) {

  //loop over events, classifying them as 2b, 3b, 4b and print output to screen
  //also return number counted in each category (via variables passed by reference)

  n2b=0;
  n3b=0;
  n4b=0;

  //loop over events
  for (size_t iev = 0; iev<dataset.size(); ++iev) {
    //for the event, need to count T,M,L
    int nL=0,nM=0,nT=0;
    vector<pair<jetFlavor,tagStatus> > event = dataset.at(iev);
    for (size_t ijet = 0; ijet<event.size(); ++ijet) {
      pair<jetFlavor,tagStatus> jet = event.at(ijet);
      if ( jet.second == kL) nL++;
      else if ( jet.second == kM) nM++;
      else if ( jet.second == kT) nT++;
    }

    double w = 1;
    if      ( nT == 2 && nM == 0 ) {
      if (reweight) w = getEventWeight_2b(event); //TODO
      n2b += w;
    }
    else if ( nT == 2 && nM == 1 && nL==0) {
      if (reweight) w = getEventWeight_LMT(event);
      n3b += w;
    }
    //this is a bit tricky! i think it's correct
    else if ( nT >= 2 && nM >=1 && (nT+nM+nL>=4)) {
      if (reweight) w = getEventWeight_LMT(event);
      n4b += w;
    }
  }

  cout<<"2b 3b 4b  "<<n2b<<" "<<n3b<<" "<<n4b<<endl;

}

double nSigmaPoisson(double a, double b) {

  double diff = a-b;
  //add in quadrature: sqrt( sqrt(a)^2 + sqrt(b)^2 )
  double err = sqrt(a +b);

  return diff/err;
}

double PercentDiff(double a, double b) {

  double diff = a-b;
  return 100*diff / a;

}

void go(const ULong64_t ngen = 1000000) {
  rnd_ = new TRandom3(123444321); //some seed

  vector< vector<pair<jetFlavor,tagStatus> >  > MC;
  vector< vector<pair<jetFlavor,tagStatus> >  > data;

  for (ULong64_t iev = 0; iev <ngen; ++iev) {
    vector<pair<jetFlavor,tagStatus> > dataevent  = generateTtbar(false);
    data.push_back(dataevent);

    vector<pair<jetFlavor,tagStatus> > mcevent  = generateTtbar(true);
    MC.push_back(mcevent);
  }

  double n2b_mc,n3b_mc,n4b_mc;
  double n2b_data,n3b_data,n4b_data;
  cout<<" == MC   =="<<endl;
  countInBins(MC,n2b_mc,n3b_mc,n4b_mc);
  cout<<" == MC reweight  =="<<endl;
  countInBins(MC,n2b_mc,n3b_mc,n4b_mc,true);
  cout<<" == data =="<<endl;
  countInBins(data,n2b_data,n3b_data,n4b_data);

  cout<<nSigmaPoisson(n2b_mc,n2b_data)<<" "<<nSigmaPoisson(n3b_mc,n3b_data)<<" "<<nSigmaPoisson(n4b_mc,n4b_data)<<endl;
  cout<<PercentDiff(n2b_mc,n2b_data)<<" "<<PercentDiff(n3b_mc,n3b_data)<<" "<<PercentDiff(n4b_mc,n4b_data)<<endl;

}
