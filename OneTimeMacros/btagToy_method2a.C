#include <vector>
#include <iostream>
#include <utility>

#include "TRandom3.h"

/*
test out the Owen/Pawandeep implementation of "Method 2a"

note that this code treats tagStatus as mutually exclusive.
i.e. if a jet is tagged as kT, then it is not kM
*/

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

tagStatus promote(tagStatus ts) {
  return nextTagStatus(ts);
}
tagStatus demote(tagStatus ts) {

  if       (ts == k0 ) assert(0);
  else  if (ts == kL ) return k0;
  else  if (ts == kM ) return kL;
  else  if (ts == kT ) return kM;
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


pair<jetFlavor,tagStatus> retagjet( const pair<jetFlavor,tagStatus> & jetin) {

  //need MC eff and SF for the jet
  double sft = getSF(jetin.first, kT);
  double sfm = getSF(jetin.first, kM);
  double sfl = getSF(jetin.first, kL);

  double efft = getEffMC(jetin.first, kT);
  double effm = getEffMC(jetin.first, kM);
  double effl = getEffMC(jetin.first, kL);

  //verify owen's assumptions

  // 0 <= eff_T <= eff_M <= eff_L <= 1
  //this is hard-coded to be true
  //now check sf*eff
  if (sfl*effl > 1) sfl=1.0/effl;
  if (sfm*effm > sfl*effl) sfm = sfl*effl / effm; 
  if (sft*efft > sfm*effm) sft = sfm*effm / efft;
  if (sft*efft < 0) sft=0;

  //now adjust the jet
  //the jet starts with the tag it has
  tagStatus newtag = jetin.second;
  //now take pawandeep's code

  //the -- and ++ operators used by pawandeep should work but let's be safer
           	if ( (sft<1)&&(sfm<1)&&(sfl<1) ) {
               		if ((newtag==kT) && (rnd_->Rndm() < (1-sft)) ) newtag = demote(newtag);
               		if ((newtag==kM) && (rnd_->Rndm() < ((1-sfm)/(1-sft*(efft/effm)))) ) newtag = demote(newtag);
               		if ((newtag==kL) && (rnd_->Rndm() < ((1-sfl)/(1-sfm*(effm/effl)))) ) newtag = demote(newtag);
           	}//1 
           	else if ( (sft<1)&&(sfm<1)&&(sfl>1) ) {
               		if ((newtag==kT) && (rnd_->Rndm() < (1-sft)) ) newtag = demote(newtag);
               		if ((newtag==kM) && (rnd_->Rndm() < ((1-sfm)/(1-sft*(efft/effm)))) ) newtag = demote(newtag);
               		if ((newtag==k0) && (rnd_->Rndm() < ((sfl-1)/((1./effl)-1))) ) newtag = promote(newtag);
           	}//2
           	else if ( (sft>1)&&(sfm<1)&&(sfl<1) ) {
               		if ((newtag==kM) && (rnd_->Rndm() < ((sft-1)/((effm/efft)-1))) ) newtag = promote(newtag);
               		if ((newtag==kM) && (rnd_->Rndm() < ((1-sfm)/(1-sft*(efft/effm)))) ) newtag = demote(newtag);
               		if ((newtag==kL) && (rnd_->Rndm() < ((1-sfl)/(1-sfm*(effm/effl)))) ) newtag = demote(newtag);
           	}//3
           	else if ( (sft>1)&&(sfm<1)&&(sfl>1) ) {
               		if ((newtag==kM) && (rnd_->Rndm() < ((sft-1)/((effm/efft)-1))) ) newtag = promote(newtag);
               		if ((newtag==kM) && (rnd_->Rndm() < ((1-sfm)/(1-sft*(efft/effm)))) ) newtag = demote(newtag);
               		if ((newtag==k0) && (rnd_->Rndm() < ((sfl-1)/((1./effl)-1))) ) newtag = promote(newtag);
           	}//4
           	else if ( (sft<1)&&(sfm>1)&&(sfl<1) ) {
               		if ((newtag==kT) && (rnd_->Rndm() < (1-sft)) ) newtag = demote(newtag);
               		if ((newtag==kL) && (rnd_->Rndm() < ((sfm-1)/((effl/effm)-1))) ) newtag = promote(newtag);
               		if ((newtag==kL) && (rnd_->Rndm() < ((1-sfl)/(1-sfm*(effm/effl)))) ) newtag = demote(newtag);
           	}//5
           	else if ( (sft<1)&&(sfm>1)&&(sfl>1) ) {
               		if ((newtag==kT) && (rnd_->Rndm() < (1-sft)) ) newtag = demote(newtag);
               		if ((newtag==k0) && (rnd_->Rndm() < ((sfl-1)/((1./effl)-1))) ) newtag = promote(newtag);
               		if ((newtag==kL) && (rnd_->Rndm() < ((sfm-1)/(sfl*(effl/effm)-1))) ) newtag = promote(newtag);
           	}//6
           	else if ( (sft>1)&&(sfm>1)&&(sfl<1) ) {
               		if ((newtag==kL) && (rnd_->Rndm() < (sfm-1)/((effl/effm)-1)) ) newtag = promote(newtag);
               		if ((newtag==kM) && (rnd_->Rndm() < ((sft-1)/(sfm*(effm/efft)-1))) ) newtag = promote(newtag);
               		if ((newtag==kL) && (rnd_->Rndm() < ((1-sfl)/(1-sfm*(effm/effl)))) ) newtag = demote(newtag);
           	}//7
           	else if ( (sft>1)&&(sfm>1)&&(sfl>1) ) {
               		if ((newtag==k0) && (rnd_->Rndm() < (sfl-1)/((1./effl)-1)) ) newtag = promote(newtag);
               		if ((newtag==kL) && (rnd_->Rndm() < ((sfm-1)/(sfl*(effl/effm)-1))) ) newtag = promote(newtag);
               		if ((newtag==kM) && (rnd_->Rndm() < ((sft-1)/(sfm*(effm/efft)-1))) ) newtag = promote(newtag);
           	}//8
           
		return make_pair(jetin.first,newtag);
}

vector<pair<jetFlavor,tagStatus> > retagevent( const vector<pair<jetFlavor,tagStatus> > & eventin) {

  vector<pair<jetFlavor,tagStatus> > eventout;
  //loop over jets
  for (size_t ijet = 0; ijet<eventin.size() ; ijet++) {
    pair<jetFlavor,tagStatus> thisjet = eventin.at(ijet);

    pair<jetFlavor,tagStatus> adjustedjet = retagjet( thisjet);
    eventout.push_back(adjustedjet);
  }
  return eventout;

}

vector< vector<pair<jetFlavor,tagStatus> >  > method2a( const  vector< vector<pair<jetFlavor,tagStatus> >  > & input) {
  vector< vector<pair<jetFlavor,tagStatus> >  > outputevents;
  // loop over events
  for ( size_t ievent = 0; ievent < input.size(); ievent++) {
    vector<pair<jetFlavor,tagStatus> > thisevent = input.at(ievent);

    vector<pair<jetFlavor,tagStatus> > adjustedevent = retagevent(thisevent);

    outputevents.push_back(adjustedevent);
  }
  return outputevents;

}


//(newly corrected)
void countInBins( const vector< vector<pair<jetFlavor,tagStatus> >  > & dataset/*,double & n2b,double & n3b,double & n4b,const bool reweight=false*/) {

  //loop over events, classifying them as 2b, 3b, 4b and print output to screen
  //also return number counted in each category (via variables passed by reference)

  ULong64_t n2b=0;
  ULong64_t n3b=0;
  ULong64_t n4b=0;

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

    if      ( nT == 2 && nM == 0 ) {
      n2b++;
    }
    else if (nT >=2 && (nM+nT)==3 && nL==0 ) {
      n3b++;
    }
    else if ( nT >= 2 && (nM+nT >=3) && (nT+nM+nL>=4)) {
      n4b++;
    }
  }

  cout<<"2b 3b 4b  "<<n2b<<" "<<n3b<<" "<<n4b<<endl;

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

  //now adjust the MC
  vector< vector<pair<jetFlavor,tagStatus> >  > MCpost = method2a(MC);

  cout<<" raw MC "<<endl;
  countInBins(MC);
  cout<<" MC after recipe "<<endl;
  countInBins(MCpost);
  cout<<" 'data'"<<endl;
  countInBins(data);


}
