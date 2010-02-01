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
#include <map>
#include <string>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>

/*
to do:
all sort of protections are needed...

also, i should make sure that all interfaces accept all 3 possibilites:
const char* ("hello world")
std::string
TString
*/

class HistHolder {
public :

  HistHolder();
  virtual ~HistHolder();

  //convention is to use lowercase for methods that are custom to this class
  //and uppercase for those that share a name with a normal ROOT method

  void make(std::string name, std::string title, Int_t nbins, Double_t min, Double_t max);
  void make(std::string name, std::string title, Int_t nbins, Double_t min, Double_t max,std::string xtitle);
  void make(std::string name, std::string title, Int_t nbins, Double_t min, Double_t max,
	    std::string xtitle, std::string ytitle);

  void make2(std::string name, std::string title, Int_t nx, Double_t minx, Double_t maxx, 
	     Int_t ny, Double_t miny, Double_t maxy);

  //load histograms from files
  //note that if the file is closed, the histo pointer may go away
  void load(TString name, TFile* file);
  //  void load2(TString name, TFile* file);

 //select a group of histograms for the next operation
  void select(const std::string str) {select_=str;}
  void reject(const std::string str) {reject_=str;}

  //these methods operate on the selected histograms
  //use * to work with all (default)

  void normalize(); //scale to unit area
  void Write();
  void Sumw2();
  void SetMinimum(Double_t min);
  void SetMaximum(Double_t max);
  Float_t GetMaximum();

  TH1F* operator[](std::string name) {return histHolder_[name]; }

  TH1F* find(std::string name) {return histHolder_[name];}
  TH1F* find(std::string name, TFile* file) {return histHolderP_[make_pair(name,file)];}
  TH2F* find2(std::string name) {return histHolder2_[name];}
 
  TH1F* find(TString name) {return find(std::string(name.Data()));}
  TH1F* find(TString name, TFile* file) {return find( std::string(name.Data()),file); }
  TH2F* find2(TString name) {return find2( std::string(name.Data()));} 
 
  void print();

private :
  bool passesFilter(const std::string mystr, TFile* filep=0);
  bool passesFilter( std::pair<std::string , TFile*> mypair) {return passesFilter( mypair.first,mypair.second);}

  std::map< std::string, TH1F*> histHolder_;
  std::map< std::string, TH2F*> histHolder2_;

  std::map< std::pair< std::string, TFile*>, TH1F* > histHolderP_;

  std::string select_;
  std::string reject_;
};

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
}

void
HistHolder::print() {

  std::cout<<"histHolder_ size = "<<histHolder_.size()<<std::endl;
  std::cout<<"histHolder2_ size = "<<histHolder2_.size()<<std::endl;
  std::cout<<"histHolderP_ size = "<<histHolderP_.size()<<std::endl;

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

void
HistHolder::load(TString name, TFile* file) {

  TH1F* hist = ((TH1F*)file->Get(name));
  histHolderP_[make_pair(std::string(name.Data()),file)] = hist;

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
  //FIXME need to add histHolder2
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

// an auxilliary class
// problem = in CINT I cannot use a std::vector<TFile*>
// this is supposed to be an easy solution, basically just a minimal interface to std::vector

class FileHolder {
public :
  FileHolder();
  virtual ~FileHolder();

  void add(TFile* filep) {files_.push_back(filep);}
  unsigned int size() {return files_.size();}

  TFile* at(unsigned int n) { return files_.at(n); }

private :

  std::vector<TFile*> files_;

};

FileHolder::FileHolder() {}

FileHolder::~FileHolder() {}
