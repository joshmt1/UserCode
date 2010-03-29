/*
to do:
all sort of protections are needed...

also, i should make sure that all interfaces accept all 3 possibilites:
const char* ("hello world")
std::string
TString
*/
#ifndef _HistHolder_h_
#define _HistHolder_h_

#include <map>
#include <string>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>

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
  Int_t load(TString name, TFile* file); //return 1=1D, 2=2D, 0=failure

  //add a histogram...

  void put(std::string name, TH1F* h) {histHolder_[name] = h;}
  void put(    TString name, TH1F* h) {histHolder_[std::string(name.Data())] = h;}


 //select a group of histograms for the next operation
  void select(const std::string str) {select_=str;}
  void reject(const std::string str) {reject_=str;}

  //these methods operate on the selected histograms
  //use * to work with all (default)

  //this group only operates on 1D histograms
  void normalize(); //scale to unit area
  void Write();
  void Sumw2();
  void SetMinimum(Double_t min);
  void SetMaximum(Double_t max);
  Float_t GetMaximum();
  TString getMaximumName();
  //only for 2D histograms
  //to be implemented!
  //  void ProjectionX(); //should add the custom name argument


  TH1F* operator[](std::string name) {return histHolder_[name]; }
  TH1F* operator[](const char* name) {return histHolder_[std::string(name)]; }
  TH1F* operator[](TString name) {return histHolder_[std::string(name.Data())]; }

  TH1F* find(std::string name) {return histHolder_[name];}
  TH1F* find(std::string name, TFile* file) {return histHolderP_[make_pair(name,file)];}
  TH2F* find2(std::string name) {return histHolder2_[name];}
  TH2F* find2(std::string name, TFile* file) {return histHolderP2_[make_pair(name,file)];} 

  void Print(TString opt="");
  void reset();

  //these are just to allow use of char* or TString
  TH1F* find(TString name) {return find(std::string(name.Data()));}
  TH1F* find(TString name, TFile* file) {return find( std::string(name.Data()),file); }
  TH2F* find2(TString name) {return find2( std::string(name.Data()));} 
  TH2F* find2(TString name, TFile* file) {return find2( std::string(name.Data()),file); }

  TH1F* find(const char* name) {return find(std::string(name));}
  TH1F* find(const char* name, TFile* file) {return find( std::string(name),file); }
  TH2F* find2(const char* name) {return find2( std::string(name));} 
  TH2F* find2(const char* name, TFile* file) {return find2( std::string(name),file); }
 
private :
  bool passesFilter(const std::string mystr, TFile* filep=0);
  bool passesFilter( std::pair<std::string , TFile*> mypair) {return passesFilter( mypair.first,mypair.second);}

  std::map< std::string, TH1F*> histHolder_;
  std::map< std::string, TH2F*> histHolder2_;

  std::map< std::pair< std::string, TFile*>, TH1F* > histHolderP_;
  std::map< std::pair< std::string, TFile*>, TH2F* > histHolderP2_;

  std::string select_;
  std::string reject_;
};

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

#endif
