/*
e.g.
 $ root -b -l -q 'EdgeFit_runoneunbinned.C(3000,"DelphesEdge/datasets_3000_109000.root",-1)' >&! DelphesEdge/run_one_unbinned_3000_109000.log &
*/
void EdgeFit_runoneunbinned(float lumi,TString filename,float constraint) {

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

  cout<<"Runing fit in fixed mass mode - 70 GeV"<<endl;
  full_fit(lumi,filename,constraint,false,70);


}

