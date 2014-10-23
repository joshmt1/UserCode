// -*- C++ -*-
#include "ConfigurationDescriptions.h"

#include "TObjArray.h"
#include <iostream>
#include <utility>
#include <cassert>

ConfigurationDescriptions::ConfigurationDescriptions() : 
  default_(""),corrected_("") { }
ConfigurationDescriptions::~ConfigurationDescriptions() {}

TString ConfigurationDescriptions::at(const unsigned int i) {
  if (i == 0) return getDefault();
  else  if (i == 1) return getCorrected();
  else {
    const unsigned int j=i-2;
    unsigned int k=0;
    for (std::map<TString, std::pair<TString,TString> >::iterator iconfig=variationPairs.begin(); iconfig!=variationPairs.end(); ++iconfig) {
      if (j == k) return iconfig->second.first;
      else if (j == k+1) return iconfig->second.second;
      k+=2;
    }
  }

  std::cout<<"WARNING in ConfigurationDescriptions::at() -- asked for element "<<i<<" when there are only "<<this->size()<<" elements!"<<std::endl;
  return "";
}

unsigned int ConfigurationDescriptions::size() {
  unsigned int s=0;
  if (default_ != "") ++s;
  if (corrected_ != "") ++s;

  s+= 2*variationPairs.size();

  return s;
}

TString ConfigurationDescriptions::getVariedSubstring(const TString & currentVariation) {

  TObjArray* baseline = corrected_.Tokenize("_");

  TObjArray* mine = currentVariation.Tokenize("_");

  TString output="";
  for (int i=0; i<baseline->GetEntries(); i++) {
    TString b=baseline->At(i)->GetName();
    TString m=mine->At(i)->GetName();
    if (b!=m) {
      output+=b;
    }
  }
  return output;
}

void ConfigurationDescriptions::addVariation(const TString & description1, const TString & description2) {

  TString var1=getVariedSubstring(description1);
  TString var2=getVariedSubstring(description2);

  assert(var1 == var2);

  variationPairs[var1] = std::make_pair(description1,description2);

}
