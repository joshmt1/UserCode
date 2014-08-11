#include "CrossSections.hh"

#include <iostream>
#include <cassert>

#include "TRegexp.h"

CrossSections::CrossSections(const TString samplename) :
  xs_(1),
  proc_(kNone)
{
  // from https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsWG/Phase2UpgradeStudies
  if (samplename.Contains("B-4p-0-1-v1510_14TEV")) xs_=200944.681290000;
  else if (samplename.Contains("BB-4p-0-300-v1510_14TEV")) xs_=249.977100000;
  else if (samplename.Contains("BB-4p-300-700-v1510_14TEV")) xs_=35.230620000;
  else if (samplename.Contains("BB-4p-700-1300-v1510_14TEV")) xs_=4.137430000;
  else if (samplename.Contains("BB-4p-1300-2100-v1510_14TEV")) xs_=0.417020000;
  else if (samplename.Contains("BB-4p-2100-100000-v1510_14TEV")) xs_=0.047700000;
  else if (samplename.Contains("BBB-4p-0-600-v1510_14TEV")) xs_=2.573040000;
  else if (samplename.Contains("BBB-4p-600-1300-v1510_14TEV")) xs_=0.149350000;
  else if (samplename.Contains("BBB-4p-1300-100000-v1510_14TEV")) xs_=0.012740000;
  else if (samplename.Contains("Bj-4p-0-300-v1510_14TEV")) xs_=34409.923390000;
  else if (samplename.Contains("Bj-4p-300-600-v1510_14TEV")) xs_=2642.853090000;
  else if (samplename.Contains("Bj-4p-600-1100-v1510_14TEV")) xs_=294.123110000;
  else if (samplename.Contains("Bj-4p-1100-1800-v1510_14TEV")) xs_=25.950000000;
  else if (samplename.Contains("Bj-4p-1800-2700-v1510_14TEV")) xs_=2.421110000;
  else if (samplename.Contains("Bj-4p-2700-3700-v1510_14TEV")) xs_=0.226900000;
  else if (samplename.Contains("Bj-4p-3700-100000-v1510_14TEV")) xs_=0.027670000;
  else if (samplename.Contains("Bjj-vbf-4p-0-700-v1510_14TEV")) xs_=86.456040000;
  else if (samplename.Contains("Bjj-vbf-4p-700-1400-v1510_14TEV")) xs_=4.348690000;
  else if (samplename.Contains("Bjj-vbf-4p-1400-2300-v1510_14TEV")) xs_=0.324650000;
  else if (samplename.Contains("Bjj-vbf-4p-2300-3400-v1510_14TEV")) xs_=0.030320000;
  else if (samplename.Contains("Bjj-vbf-4p-3400-100000-v1510_14TEV")) xs_=0.003130000;
  else if (samplename.Contains("H-4p-0-300-v1510_14TEV")) xs_=21.559900000;
  else if (samplename.Contains("H-4p-300-800-v1510_14TEV")) xs_=1.112820000;
  else if (samplename.Contains("H-4p-800-1500-v1510_14TEV")) xs_=0.091880000;
  else if (samplename.Contains("H-4p-1500-100000-v1510_14TEV")) xs_=0.010090000;
  else if (samplename.Contains("LL-4p-0-100-v1510_14TEV")) xs_=1341.369230000;
  else if (samplename.Contains("LL-4p-100-200-v1510_14TEV")) xs_=156.295340000;
  else if (samplename.Contains("LL-4p-200-500-v1510_14TEV")) xs_=42.401320000;
  else if (samplename.Contains("LL-4p-500-900-v1510_14TEV")) xs_=2.843730000;
  else if (samplename.Contains("LL-4p-900-1400-v1510_14TEV")) xs_=0.209140000;
  else if (samplename.Contains("LL-4p-1400-100000-v1510_14TEV")) xs_=0.028910000;
  else if (samplename.Contains("LLB-4p-0-400-v1510_14TEV")) xs_=2.973800000;
  else if (samplename.Contains("LLB-4p-400-900-v1510_14TEV")) xs_=0.228540000;
  else if (samplename.Contains("LLB-4p-900-100000-v1510_14TEV")) xs_=0.020800000;
  else if (samplename.Contains("tB-4p-0-500-v1510_14TEV")) xs_=63.889230000;
  else if (samplename.Contains("tB-4p-500-900-v1510_14TEV")) xs_=7.121720000;
  else if (samplename.Contains("tB-4p-900-1500-v1510_14TEV")) xs_=0.980300000;
  else if (samplename.Contains("tB-4p-1500-2200-v1510_14TEV")) xs_=0.083910000;
  else if (samplename.Contains("tB-4p-2200-100000-v1510_14TEV")) xs_=0.009530000;
  else if (samplename.Contains("tj-4p-0-500-v1510_14TEV")) xs_=109.736020000;
  else if (samplename.Contains("tj-4p-500-1000-v1510_14TEV")) xs_=5.993250000;
  else if (samplename.Contains("tj-4p-1000-1600-v1510_14TEV")) xs_=0.376800000;
  else if (samplename.Contains("tj-4p-1600-2400-v1510_14TEV")) xs_=0.034620000;
  else if (samplename.Contains("tj-4p-2400-100000-v1510_14TEV")) xs_=0.003120000;
  else if (samplename.Contains("tt-4p-0-600-v1510_14TEV")) xs_=530.893580000;
  else if (samplename.Contains("tt-4p-600-1100-v1510_14TEV")) xs_=42.553510000;
  else if (samplename.Contains("tt-4p-1100-1700-v1510_14TEV")) xs_=4.482090000;
  else if (samplename.Contains("tt-4p-1700-2500-v1510_14TEV")) xs_=0.527950000;
  else if (samplename.Contains("tt-4p-2500-100000-v1510_14TEV")) xs_=0.054490000;
  else if (samplename.Contains("ttB-4p-0-900-v1510_14TEV")) xs_=2.667300000;
  else if (samplename.Contains("ttB-4p-900-1600-v1510_14TEV")) xs_=0.250469000;
  else if (samplename.Contains("ttB-4p-1600-2500-v1510_14TEV")) xs_=0.023744100;
  else if (samplename.Contains("ttB-4p-2500-100000-v1510_14TEV")) xs_=0.002088160;
  else if (samplename.Contains("stoc")) xs_=0.009635;
  else {
    //try regexp for more complicated cases
    /*
in particular, the example use case is:

./FlatTree "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/SUSY_SIGNAL/140PileUp/naturalModel//1/wildcard.root" -O ${workdir}/${outfilename} -N 1 1

so we're looking for naturalModel/digit but need to allow for extra / in between
     */
    TRegexp reg1("naturalModel/+[0-9]");
    TString match = samplename(reg1);
    if (match.Length()>0) { //look at the last digit
      int lastdigit=TString(match(match.Length()-1)).Atoi();
      //cross-sections updated according to 2 July email from Karim
      if      (lastdigit==1) xs_=0.07; //Natural Model 1 aka Scenario 1
      else if (lastdigit==2) xs_=0.05; //Natural Model 2
      else if (lastdigit==3) xs_=1.098;  //Natural Model 3
      else 
	std::cout<<" WARNING -- this sample has no known cross-section. Using value of "<<xs_<<std::endl;
    }
    else
      std::cout<<" WARNING -- this sample has no known cross-section. Using value of "<<xs_<<std::endl;
  }

  //now set proc
  SetProc(samplename);
}

CrossSections::~CrossSections() {
 
}

void CrossSections::SetProc(TString name) {

  //set proc_
  if (name.Contains("susy")) proc_=kSignal;
  else if (name.Contains("scenario")) proc_=kSignal;
  else if (name.Contains("stoc")) proc_=kSignal;
  else if (name.Contains("natural")) proc_=kSignal;
  else if (name.Contains("BBB-")) proc_ = kRare;
  else if (name.Contains("Bjj-")) proc_=kBoson;
  else if (name.Contains("ttB-")) proc_=kRare;
  else if (name.Contains("LLB-")) proc_=kRare;
  else if (name.Contains("LL-")) proc_=kBoson;
  else if (name.Contains("Bj-")) proc_=kBoson;
  else if (name.Contains("BB-")) proc_ = kRare;
  else if (name.Contains("tB-")) proc_=kTop;
  else if (name.Contains("tj-")) proc_=kTop;
  else if (name.Contains("tt-")) proc_=kTop;
  else if (name.Contains("H-")) proc_=kRare;
  else if (name.Contains("B-")) proc_=kBoson;
  else proc_=kNone;

  std::cout<<"SetProc: "<<name<<" "<<proc_<<std::endl;
  if (proc_==kNone) assert(0);

}
