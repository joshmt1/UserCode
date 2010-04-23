/*
Goal of this class is to make a class to hold a number of histograms

The base data structure is a map using a string (the name of the histogram)
as a key, that references to a pointer that holds the histogram

We also have a secondary data structure to map by both string and
TFile*, to allow this class to hold a bunch of histograms with the
same name

Utilities should be provided that call Sumw2() and do other useful things
Creating a histogram should be as easy calling a function and passing the name, title, and
other standard arguments

*/

#include <iostream>

#include "HistHolder.h"

HistHolder::HistHolder() :
  select_("*"),
  reject_("")
{

}

HistHolder::~HistHolder()
{
  /* causes a crash
  for ( std::map< std::string, TH1F*>::const_iterator i = histHolder_.begin() ; i!= histHolder_.end() ; ++i )
    delete i->second;
  */
  reset();
}

void HistHolder::reset() {
  //  std::cout<<"HistHolder reset()"<<std::endl;
  histHolder_.clear();
  histHolder2_.clear();
  histHolderP_.clear();
  histHolderP2_.clear();
}

void
HistHolder::Print(TString opt) {

  opt.ToUpper();
  if (opt=="V") {
    std::cout<<"histHolder_ size = "<<histHolder_.size()<<std::endl;
    for ( std::map< std::string, TH1F*>::const_iterator i = histHolder_.begin() ; i!= histHolder_.end() ; ++i ) {
      std::cout<<i->first<<"\t"<<i->second<<std::endl;
    }
    std::cout<<"histHolder2_ size = "<<histHolder2_.size()<<std::endl;
    for ( std::map< std::string, TH2F*>::const_iterator i = histHolder2_.begin() ; i!= histHolder2_.end() ; ++i ) {
      std::cout<<i->first<<"\t"<<i->second<<std::endl;
    }
    std::cout<<"histHolderP_ size = "<<histHolderP_.size()<<std::endl;
    for ( std::map< std::pair<std::string, TFile*>, TH1F*>::const_iterator i = histHolderP_.begin() ; i!= histHolderP_.end() ; ++i ) {
      std::cout<<i->first.first<<" ; "<<i->first.second<<"\t"<<i->second<<std::endl;
    }
    std::cout<<"histHolderP2_ size = "<<histHolderP2_.size()<<std::endl;
    for ( std::map< std::pair<std::string, TFile*>, TH2F*>::const_iterator i = histHolderP2_.begin() ; i!= histHolderP2_.end() ; ++i ) {
      std::cout<<i->first.first<<" ; "<<i->first.second<<"\t"<<i->second<<std::endl;
    }
      
  }
  else {
    std::cout<<"histHolder_ size = "<<histHolder_.size()<<std::endl;
    std::cout<<"histHolder2_ size = "<<histHolder2_.size()<<std::endl;
    std::cout<<"histHolderP_ size = "<<histHolderP_.size()<<std::endl;
    std::cout<<"histHolderP2_ size = "<<histHolderP2_.size()<<std::endl;
  }
  
}

void HistHolder::Reset() {


  for ( std::map< std::string, TH1F*>::const_iterator i = histHolder_.begin() ; i!= histHolder_.end() ; ++i ) {
    i->second->Reset();
  }

  for ( std::map< std::string, TH2F*>::const_iterator i = histHolder2_.begin() ; i!= histHolder2_.end() ; ++i ) {
    i->second->Reset();
  }

  for ( std::map< std::pair<std::string, TFile*>, TH1F*>::const_iterator i = histHolderP_.begin() ; i!= histHolderP_.end() ; ++i ) {
    i->second->Reset();
  }

  for ( std::map< std::pair<std::string, TFile*>, TH2F*>::const_iterator i = histHolderP2_.begin() ; i!= histHolderP2_.end() ; ++i ) {
    i->second->Reset();
  }
      
}

void
HistHolder::make(std::string name, std::string title, Int_t nbins, Double_t min, Double_t max) {

  TH1F* hist = new TH1F(name.c_str(),title.c_str(),nbins,min,max);
  histHolder_[name] = hist;

}

void 
HistHolder::make2(std::string name, std::string title, Int_t nx, Double_t minx, Double_t maxx, 
		  Int_t ny, Double_t miny, Double_t maxy) {
  TH2F* hist = new TH2F(name.c_str(),title.c_str(),nx,minx,maxx,ny,miny,maxy);
  histHolder2_[name] = hist;
}

Int_t
HistHolder::load(TString name, TFile* file) {

  Int_t code=0;

  TObject* h = file->Get(name);
  if (h==0) {std::cout<<"Problem loading "<<name<<std::endl; return code;}

  if ( h->InheritsFrom("TH2F") ) {
    histHolderP2_[make_pair(std::string(name.Data()),file)] = (TH2F*) h;
    code=2;
  }
  else if ( h->InheritsFrom("TH1F") ) {
    histHolderP_[make_pair(std::string(name.Data()),file)] = (TH1F*) h;
    code=1;
  }
  else {
    std::cout<<"Could not determine type of histogram "<<name<<std::endl;
  }

  return code;
}

// void
// HistHolder::load2(TString name, const TFile* file) {

//   TH2F* hist = ((TH2F*)file->Get(name));
//   histHolder2_[std::string(name.Data())] = hist;

// }

void
HistHolder::make(std::string name, std::string title, Int_t nbins, Double_t min, Double_t max, std::string xtitle) {

  make(name,title,nbins,min,max);
  histHolder_[name]->SetXTitle(xtitle.c_str());
}

void
HistHolder::make(std::string name, std::string title, Int_t nbins, Double_t min, Double_t max,
		 std::string xtitle, std::string ytitle) {

  make(name,title,nbins,min,max,xtitle);
  histHolder_[name]->SetYTitle(ytitle.c_str());
}

void
HistHolder::Write() {
  //FIXME need to add histHolder2 and other histHolders
  for ( std::map< std::string, TH1F*>::const_iterator i = histHolder_.begin() ; i!= histHolder_.end() ; ++i ) {
    if ( passesFilter( i->first))
      i->second->Write();
  }
  
}

void HistHolder::Sumw2() {

  for ( std::map< std::string, TH1F*>::const_iterator i = histHolder_.begin() ; i!= histHolder_.end() ; ++i ) {
    if ( passesFilter( i->first))
      i->second->Sumw2();
  }

  for ( std::map< std::pair<std::string, TFile*>, TH1F*>::const_iterator i = histHolderP_.begin() ;
	i!= histHolderP_.end() ; ++i ) {
    
    if ( passesFilter(i->first))
      i->second->Sumw2();
  }

}

void HistHolder::SetMaximum( Double_t max) {
  for ( std::map< std::string, TH1F*>::const_iterator i = histHolder_.begin() ; i!= histHolder_.end() ; ++i ) {
    if ( passesFilter(i->first)) 
      i->second->SetMaximum(max);
  }

  for ( std::map< std::pair<std::string, TFile*>, TH1F*>::const_iterator i = histHolderP_.begin() ;
	i!= histHolderP_.end() ; ++i ) {
    if ( passesFilter( i->first))
      i->second->SetMaximum(max);
  }

}

void HistHolder::SetMinimum( Double_t min) {
  for ( std::map< std::string, TH1F*>::const_iterator i = histHolder_.begin() ; i!= histHolder_.end() ; ++i ) {
    if ( passesFilter(i->first)) 
      i->second->SetMinimum(min);
  }
  
  for ( std::map< std::pair<std::string, TFile*>, TH1F*>::const_iterator i = histHolderP_.begin() ;
	i!= histHolderP_.end() ; ++i ) {
    if ( passesFilter(i->first))
      i->second->SetMinimum(min);
  }
  
}

Float_t HistHolder::GetMaximum() {

  Float_t max = -99999999;

  for ( std::map< std::string, TH1F*>::const_iterator i = histHolder_.begin() ; i!= histHolder_.end() ; ++i ) {
    if ( passesFilter(i->first)) {
      if ( i->second->GetMaximum() > max ) max = i->second->GetMaximum();
    }
  }

  for ( std::map< std::pair<std::string, TFile*>, TH1F*>::const_iterator i = histHolderP_.begin() ;
	i!= histHolderP_.end() ; ++i ) {
    if ( passesFilter(i->first)) {
      if ( i->second->GetMaximum() > max ) max = i->second->GetMaximum();
    }
  }
  
  return max;
}

TString HistHolder::getMaximumName() {

  //most code copied from GetMaximum()
  Float_t max = -99999999;
  TString nameofmax="";

  for ( std::map< std::string, TH1F*>::const_iterator i = histHolder_.begin() ; i!= histHolder_.end() ; ++i ) {
    if ( passesFilter(i->first)) {
      if ( i->second->GetMaximum() > max ) {
	max = i->second->GetMaximum();
	nameofmax = TString(i->first);
      }
    }
  }

  for ( std::map< std::pair<std::string, TFile*>, TH1F*>::const_iterator i = histHolderP_.begin() ;
	i!= histHolderP_.end() ; ++i ) {
    if ( passesFilter(i->first)) {
      if ( i->second->GetMaximum() > max ) {
	max = i->second->GetMaximum();
	nameofmax = TString(i->first.first);
      }
    }
  }
  
  return nameofmax;
}

void HistHolder::normalize() {
  for ( std::map< std::string, TH1F*>::const_iterator i = histHolder_.begin() ; i!= histHolder_.end() ; ++i ) {
    if ( passesFilter(i->first)) {
      double factor = i->second->Integral();
      i->second->Scale(1.0/factor);
    }
  }

  for ( std::map< std::pair<std::string, TFile*>, TH1F*>::const_iterator i = histHolderP_.begin() ;
	i!= histHolderP_.end() ; ++i ) {
    if ( passesFilter( i->first)) {
      double factor = i->second->Integral();
      i->second->Scale(1.0/factor);
    }
  }
  
}

// utility functions

bool HistHolder::passesFilter(const std::string mystr, TFile* filep) {

  //if reject_ is empty, then reject nothing
  //if select_ is empty, then accept everything
  if (select_=="") select_="*";

  std::string mystr2="";
  if (filep!=0) mystr2 = std::string(filep->GetName());

  bool rejected=false;
  if (reject_ != "") {
    bool rejected1 = (mystr.find(reject_) != std::string::npos);
    bool rejected2=false;
    if (filep !=0) {
      rejected2 = (mystr2.find(reject_) != std::string::npos);
    }
    rejected = rejected1 || rejected2;
  }

  if (rejected) return false;

  if (select_ == "*") return true;

//   std::cout<<"input string = "<<mystr<<std::endl;
//   std::cout<<"compare str  = "<<select_<<std::endl;
//   std::cout<<"found it = "<<(  mystr.find(select_) != std::string::npos )<<std::endl;

  bool pass1 = mystr.find(select_) != std::string::npos ;
  bool pass2=false;
  if ( filep!=0) {
    pass2 = (mystr2.find(select_) != std::string::npos);
  }

  return pass1 || pass2;

}


FileHolder::FileHolder() {}

FileHolder::~FileHolder() {
  files_.clear();

}
