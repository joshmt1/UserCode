// -*- C++ -*-
//container for flat root tree

#ifndef SIMPLETREE_H
#define SIMPLETREE_H

#include <map>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"

class SimpleTree {

public:
  SimpleTree(const TString & outfilename);
  ~SimpleTree();

  //handles float,double,bool and fixed array of float
  //that's not everything but it is good enough for just about anything in practice

  void AddVariable(const TString & name, float dummyVal = -99);
  void AddDouble(const TString & name, double dummyVal = -99);
  void AddBool(const TString & name, bool dummyVal=false);
  void AddInt(const TString & name, int dummyVal=-99);
  void AddArray(const TString & name, int length, float dummyVal = -99);
  void Fill();
  void Set(const TString & name, float val);
  void SetDouble(const TString & name, double val);
  void SetBool(const TString & name, bool val);
  void SetInt(const TString & name, int val);
  void Set(const TString & name,int index, float val); //for arrays
  float Get(const TString & name);


private:
  TFile* outputFile_;
  TTree* tree_;
  std::map<TString, float> vars_;
  std::map<TString,float> vars_dummyVals_;  //keep track of dummy values for vars
  std::map<TString, double> doubles_;
  std::map<TString,double> doubles_dummyVals_;  //keep track of dummy values for doubles
  std::map<TString, bool> bools_;
  std::map<TString,bool> bools_dummyVals_;  //keep track of dummy values for bools
  std::map<TString, int> ints_;
  std::map<TString,int> ints_dummyVals_;  //keep track of dummy values for ints
  std::map<TString, float*> arrays_;
  std::map<TString,float> arrays_dummyVals_; //again, keep track of dummy vals
  std::map<TString,int> arrays_length_; //also lengths

  //new feature idea -- add check at Fill() that the variable has been set since the last Fill().
  //this could help avoid mistakes.
  //  std::map<TString,bool> isVarSet_;
  //this would also allow us to catch, at Add time, whether a variable name is already taken.

};

#endif
