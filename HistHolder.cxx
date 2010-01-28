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
#include <map>
#include <string>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>

#define DUMMYSTR "ignoreme123"

/*
to do:
all sort of protections are needed...
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
  //void load(TString name, const TFile* file);
  //  void load2(TString name, const TFile* file);

 //select a group of histograms for the next operation
  void select(std::string name) {select_=name;}

  //these methods operate on the selected histograms
  //use * to work with all 
  void Write();
  void Sumw2();
  void SetMaximum(Double_t max);

  TH1F* operator[](std::string name) {return histHolder_[name]; }

  TH1F* find(std::string name) {return histHolder_[name];}
  TH2F* find2(std::string name) {return histHolder2_[name];}
  
private :
  std::map< std::string, TH1F*> histHolder_;
  std::map< std::string, TH2F*> histHolder2_;

  //std::map< std::pair< std::string, TFile*>, TH1F* > histHolderP_;

  std::string select_;
};

HistHolder::HistHolder() :
  select_("")
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

// void
// HistHolder::load(TString name, const TFile* file) {

//   TH1F* hist = ((TH1F*)file->Get(name));
//   histHolderP_[make_pair(std::string(name.Data()),file)] = hist;

// }

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
    if ( select_ == "*" || i->first.find(select_)!=std::string::npos) 
      i->second->Write();
  }
  
}

void HistHolder::Sumw2() {

  for ( std::map< std::string, TH1F*>::const_iterator i = histHolder_.begin() ; i!= histHolder_.end() ; ++i ) {
    if ( select_ == "*" || i->first.find(select_)!=std::string::npos) 
      i->second->Sumw2();
  } 
}

void HistHolder::SetMaximum( Double_t max) {
  for ( std::map< std::string, TH1F*>::const_iterator i = histHolder_.begin() ; i!= histHolder_.end() ; ++i ) {
    if ( select_ == "*" || i->first.find(select_)!=std::string::npos) 
      i->second->SetMaximum(max);
  }

}
