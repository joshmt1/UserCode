/*
e.g.
 $ root -b -l -q 'EdgeFit_runbulkunbinned.C(3000,-1)' >&! DelphesEdge/run_bulk_unbinned_3000.log &
*/
void EdgeFit_runbulkunbinned(float l,float constraint) {

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

  full_fit_unbinned_bulk(l,constraint);

}

