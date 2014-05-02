#include "CrossSections.hh"

#include <iostream>
#include <cassert>

CrossSections::CrossSections(const TString samplename) :
  xs_(1)
  //  samplename_(samplename)
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
  else {
    std::cout<<" WARNING -- this sample has no known cross-section. Using value of "<<xs_<<std::endl;
  }
}

CrossSections::~CrossSections() {
 
}

