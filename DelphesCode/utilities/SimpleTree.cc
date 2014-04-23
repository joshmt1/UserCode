#include "SimpleTree.hh"

#include <iostream>

using namespace std;

SimpleTree::SimpleTree(const TString & outfilename) :
  outputFile_(0),
  tree_(0)
{

  outputFile_ = new TFile(outfilename,"RECREATE");
  tree_ = new TTree("simpleTree","tree for analysis");
  tree_->SetMaxTreeSize(30000000000LL); //max size of ~30GB
}

SimpleTree::~SimpleTree() {
  outputFile_->Write();
  outputFile_->Close();

  //clean up arrays
  for (map<TString,float*>::iterator iter = arrays_.begin();iter!=arrays_.end(); ++iter) {
    delete [] iter->second;
  }

  //there's a crash when i try to clean up properly
  //  delete tree_;
  //  delete outputFile_;
}

void SimpleTree::Set(const TString & name, float val, bool accumulate) {
  map<TString,float>::iterator iter = vars_.find(name);
  if (iter != vars_.end() ) {
    if (accumulate)  iter->second += val;
    else             iter->second = val;
  }

}

void SimpleTree::SetDouble(const TString & name, double val) {
  map<TString,double>::iterator iter = doubles_.find(name);
  if (iter != doubles_.end() )     iter->second = val;
}

void SimpleTree::SetBool(const TString & name, bool val) {
  map<TString,bool>::iterator iter = bools_.find(name);
  if (iter != bools_.end() )     iter->second = val;
}

void SimpleTree::SetInt(const TString & name, int val,bool accumulate) {
  map<TString,int>::iterator iter = ints_.find(name);
  if (iter != ints_.end() ) {
    if (accumulate) iter->second += val;
    else            iter->second = val;
  }
}

void SimpleTree::Set(const TString & name,int index, float val) { //for arrays

  map<TString,float*>::iterator iter = arrays_.find(name);
  if (iter != arrays_.end() ) {
    iter->second[index] = val;
  }

}

float SimpleTree::Get(const TString & name) {

  map<TString,float>::iterator iter = vars_.find(name);
  if (iter != vars_.end() ) {
    return iter->second;
  }
  cout<<"Warning in Get -- not found: "<<name<<endl;
  return -1;

}

int SimpleTree::GetInt(const TString & name) {
  map<TString,int>::iterator iter = ints_.find(name);
  if (iter != ints_.end() ) {
    return iter->second;
  }
  cout<<"Warning in GetInt -- not found: "<<name<<endl;
  return -1;
}

void SimpleTree::AddVariable(const TString & name, float dummyVal) {

  if ( vars_.count(name) != 0) {
    cout<<"Variable with name "<<name<<" is already in the tree!"<<endl;
    return;
  }

  cout<<"Adding variable "<<name<<" to tree"<<endl;
  vars_[name] = dummyVal;
  vars_dummyVals_[name]=dummyVal;

  TString decoratedName = name;
  decoratedName += "/F";
  tree_->Branch(name,&vars_[name],decoratedName);

}

//i'm sure this can be done better with templating but what the hell
void SimpleTree::AddDouble(const TString & name, double dummyVal) {

  if ( doubles_.count(name) != 0) {
    cout<<"Variable with name "<<name<<" is already in the tree!"<<endl;
    return;
  }

  cout<<"Adding double "<<name<<" to tree"<<endl;
  doubles_[name] = dummyVal;
  doubles_dummyVals_[name]=dummyVal;

  TString decoratedName = name;
  decoratedName += "/D";
  tree_->Branch(name,&doubles_[name],decoratedName);

}

void SimpleTree::AddInt(const TString & name, int dummyVal) {

  if ( ints_.count(name) != 0) {
    cout<<"Variable with name "<<name<<" is already in the tree!"<<endl;
    return;
  }

  cout<<"Adding int "<<name<<" to tree"<<endl;
  ints_[name] = dummyVal;
  ints_dummyVals_[name]=dummyVal;

  TString decoratedName = name;
  decoratedName += "/I";
  tree_->Branch(name,&ints_[name],decoratedName);

}

void SimpleTree::AddBool(const TString & name, bool dummyVal) {

  if ( bools_.count(name) != 0) {
    cout<<"Variable with name "<<name<<" is already in the tree!"<<endl;
    return;
  }

  cout<<"Adding bool "<<name<<" to tree"<<endl;
  bools_[name] = dummyVal;
  bools_dummyVals_[name]=dummyVal;

  TString decoratedName = name;
  decoratedName += "/O";
  tree_->Branch(name,&bools_[name],decoratedName);

}


void SimpleTree::AddArray(const TString & name, int length, float dummyVal) {

  if ( arrays_.count(name) != 0) {
    cout<<"Variable with name "<<name<<" is already in the tree!"<<endl;
    return;
  }

  cout<<"Adding array "<<name<<" to tree"<<endl;
  arrays_[name] = new float[length];

  for (int ii = 0;ii<length;++ii) arrays_[name][ii]=dummyVal;
  arrays_dummyVals_[name]=dummyVal;
  arrays_length_[name]=length;

  TString decoratedName;
  decoratedName.Form( "%s[%d]/F",name.Data(),length);
  tree_->Branch(name,arrays_[name],decoratedName);

}

void SimpleTree::Fill() {
  //not only fills the tree but also resets to dummy vars

  tree_->Fill();

  for ( map<TString,float>::iterator ii = vars_.begin(); ii!=vars_.end(); ++ii) {
    ii->second = vars_dummyVals_[ii->first];
  }

  for ( map<TString,double>::iterator ii = doubles_.begin(); ii!=doubles_.end(); ++ii) {
    ii->second = doubles_dummyVals_[ii->first];
  }

  for ( map<TString,bool>::iterator ii = bools_.begin(); ii!=bools_.end(); ++ii) {
    ii->second = bools_dummyVals_[ii->first];
  }

  for ( map<TString,float*>::iterator ii = arrays_.begin(); ii!=arrays_.end(); ++ii) {
    for (int jj = 0; jj<arrays_length_[ii->first]; ++jj) {
      ii->second[jj] = arrays_dummyVals_[ii->first];
    }
  }

}
