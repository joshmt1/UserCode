/*
e.g.
 $ root -b -l -q 'EdgeFit_runbulk.C(3000,999999)'
*/
void EdgeFit_runbulk(float l,int max) {

  gROOT->ProcessLine(".L TSelectorMultiDraw.C+");
  gROOT->ProcessLine(".L CrossSectionTable.cxx+");
  gROOT->ProcessLine(".L ConfigurationDescriptions.cxx+");
  
  gROOT->ProcessLine(".L RooEdge.cxx+");
  gROOT->ProcessLine(".L RooEdgeFlavSym.cxx+");
  gROOT->ProcessLine(".L RooTopPairProductionSpline.cxx+");
  gROOT->ProcessLine(".L RooCruijff.cxx+");
  gROOT->ProcessLine(".L RooBifurCB.cxx+");
  gROOT->ProcessLine(".L RooFallingSpectrum.cxx+");
  
  gROOT->ProcessLine(".L EdgeFit.C+");

  full_fit_subtracted_bulk(l,max);

}

