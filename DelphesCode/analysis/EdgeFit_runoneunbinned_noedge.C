/*
e.g.
 $ root -b -l -q 'EdgeFit_runoneunbinned.C(3000,"DelphesEdge/datasets_3000_109000.root",-1)' >&! DelphesEdge/run_one_unbinned_3000_109000.log &
*/
void EdgeFit_runoneunbinned_noedge(float lumi,Long_t seed,float constraint) {

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

  TString filename;
  filename.Form("DelphesEdge/NoEdge/%d/datasets_%d_%d.root",
		seed%1000,TMath::Nint(lumi),seed);

  full_fit(lumi,filename,constraint,true,70); //last argument fixes the edge mass


}

