#include "drawDelphesBase.C"
/*
--- to compile ---

.L TSelectorMultiDraw.C+
.L CrossSectionTable.cxx+
.L ConfigurationDescriptions.cxx+

.L drawDelphes_Edge.C+

*/


void  drawDelphes_Edge(TString plotsToMake="all") {

  useNewStyle_=true;

  initSamples("tt bj nm1 skimmed");
  setOutputDirectory("DelphesEdge"); 

  treatAllAsSM_=true; //don't treat signal as signal

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false;  

  //  stackSignal_=false;
  leg_x1 = 0.68;

  stackSignal_=true;
  setStackMode(true,false,false); //stack,norm,label override
  removeSample("naturalModel1");
  addSample("naturalModel1:leptonsMatchChi2ToChi1!=1",kRed-9,"Other SUSY");
  addSample("naturalModel1:leptonsMatchChi2ToChi1==1",kRed,"True Edge Signal");

  //edge selection cuts
  TCut dileptons="mll>20";
  TCut sf = "isSF==1";
  TCut noZ = "mll<85|| mll>95";

  TCut dileptons_loose="mll_loose>20";
  TCut sf_loose = "isSF_loose==1";
  TCut noZ_loose = "mll_loose<85|| mll_loose>95";

  TCut jetsloose = "njets40eta3p0>=2";
  TCut jetstight = "njets40eta3p0>=3";
  TCut metloose = "MET>100";
  TCut mettight = "MET>150";


  //first the 8 TeV Edge, with Z veto added
  selection_ = dileptons && sf && jetsloose && mettight && noZ; //sel1

  nbins=40; low=0; high=200;
  var="mll"; xtitle="m_{l+l-} (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_sel1_mll",0,"GeV");

  //first the 8 TeV Edge, with Z veto added
  selection_ = dileptons && sf && jetstight && metloose && noZ; //sel2

  nbins=40; low=0; high=200;
  var="mll"; xtitle="m_{l+l-} (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_sel2_mll",0,"GeV");

  //then tighten
  selection_ = dileptons && sf && TCut("njets40eta3p0>=5") && TCut("HT>1000") && TCut("MET>200") && noZ && TCut("nbjets40tight>=2");
  var="mll"; xtitle="m_{l+l-} (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new1_mll",0,"GeV");

  //try without b-tags
  selection_ = dileptons && sf && TCut("njets40eta3p0>=5") && TCut("HT>1000") && TCut("MET>200") && noZ;
  var="mll"; xtitle="m_{l+l-} (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("002")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new1nob_mll",0,"GeV");

  //plot n-btag distribution
  selection_ = dileptons && sf && TCut("njets40eta3p0>=5") && TCut("HT>1000") && TCut("MET>200") && noZ;
  nbins=8; low=0; high=8;
  var="nbjets40tight"; xtitle="tight b tags";
  stackSignal_=false;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new1nob_nbT",0,"GeV");

  var="nbjets40loose"; xtitle="loose b tags";
  stackSignal_=false;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new1nob_nbL",0,"GeV");

  //add the Z pieces back in
  selection_ = dileptons && sf && TCut("njets40eta3p0>=5") && TCut("HT>1000") && TCut("MET>200");

  var="nbjets40tight"; xtitle="tight b tags";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new1nob-trueEdge_nbT",0,"GeV");

  var="nbjets40loose"; xtitle="loose b tags";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new1nob-trueEdge_nbL",0,"GeV");

  //
  nbins=40; low=0; high=800;
  var="MET"; xtitle="MET";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("004")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new1nob-trueEdge_MET",0,"GeV");

  nbins=40; low=0; high=3000;
  var="HT"; xtitle="HT";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("004")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new1nob-trueEdge_HT",0,"GeV");

  //new2
  selection_ = dileptons && sf && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1");
  nbins=40; low=0; high=200;
  var="mll"; xtitle="m_{l+l-} (GeV)";
  setPlotMaximum(600); //instead of Z veto and upper cut on mll
  if (plotsToMake.Contains("all")||plotsToMake.Contains("005")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2_mll",0,"GeV");

  selection_ = dileptons && sf && TCut("njets40eta3p0>=4") && TCut("MET>400") && TCut("nbjets40tight>=1") && TCut("mll_maxEta<1.4");
  nbins=40; low=0; high=200;
  var="mll"; xtitle="m_{l+l-} (GeV)";
  setPlotMaximum(600); //instead of Z veto and upper cut on mll
  if (plotsToMake.Contains("all")||plotsToMake.Contains("005")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2-central_mll",0,"GeV");

  removeSample("naturalModel1:leptonsMatchChi2ToChi1!=1");
  removeSample("naturalModel1:leptonsMatchChi2ToChi1==1");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose!=1",kRed-9,"Other SUSY");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1",kRed,"True Edge Signal");
  selection_ = dileptons_loose && sf_loose && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1");
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  setPlotMaximum(800); //instead of Z veto and upper cut on mll
  if (plotsToMake.Contains("all")||plotsToMake.Contains("005")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2_mllLoose",0,"GeV");

  //now the mll signal region
  TCut mll_loose_signal = "mll_loose>=20 && mll_loose<70";
  selection_ = mll_loose_signal && sf_loose && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1");
  var = "HT"; xtitle=var;
  nbins=30; low = 0; high = 3000;
  resetPlotMaximum();
  if (plotsToMake.Contains("all")||plotsToMake.Contains("006")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2mllLooseLowMass_HT",0,"GeV");

  removeSample("naturalModel1:leptonsMatchChi2ToChi1_loose!=1");
  removeSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1");
  addSample("naturalModel1:leptonsMatchChi2ToChi1!=1",kRed-9,"Other SUSY");
  addSample("naturalModel1:leptonsMatchChi2ToChi1==1",kRed,"True Edge Signal");

  TCut mll_signal = "mll>=20 && mll<70";
  selection_ = mll_signal && sf && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1");
  var = "HT"; xtitle=var;
  resetPlotMaximum();
  if (plotsToMake.Contains("all")||plotsToMake.Contains("006")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2mllLowMass_HT",0,"GeV");

  selection_ = mll_signal && sf && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1")&& TCut("mll_maxEta<1.4");
  var = "HT"; xtitle=var;
  resetPlotMaximum();
  if (plotsToMake.Contains("all")||plotsToMake.Contains("007")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2mllLowMass-central_HT",0,"GeV");

  selection_ =dileptons && sf && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1")&& TCut("mll_maxEta<1.4")&&TCut("HT>1500");
  var = "mll"; xtitle=var;
  nbins=40; low=0; high=200;
  setPlotMaximum(200); //instead of Z veto and upper cut on mll
  if (plotsToMake.Contains("all")||plotsToMake.Contains("008")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500-central_mll",0,"GeV");

  removeSample("naturalModel1:leptonsMatchChi2ToChi1!=1");
  removeSample("naturalModel1:leptonsMatchChi2ToChi1==1");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose!=1&&leptonsMatchChi4ToChi1_loose!=1",kRed-9,"Other NM1");
  addSample("naturalModel1:leptonsMatchChi4ToChi1_loose==1",kGreen+2,"#tilde{#chi}_{4}^{0} #rightarrow #tilde{l}l #rightarrow l^{+}l^{-} #tilde{#chi}_{1}^{0}");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1",kRed,"#tilde{#chi}_{2}^{0} #rightarrow #tilde{l}l #rightarrow l^{+}l^{-} #tilde{#chi}_{1}^{0}");

  setPadDimensions(800,600);

  isPreliminary_=true;
  selection_ = dileptons_loose && sf_loose && TCut("njets40>=6") && TCut("MET>450")  && TCut("nbjets40medium>=1")&&TCut("HT>1250"); //updated NOMINAL SELECTION
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  resetPlotMaximum();
  doOverflowAddition(false);//jeff's request
  if (plotsToMake.Contains("all")||plotsToMake.Contains("009")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new7_mllLoose",0,"GeV");


  //now plot OF (but otherwise the same cuts)
  selection_ = dileptons_loose && !sf_loose && TCut("njets40>=6") && TCut("MET>450")  && TCut("nbjets40medium>=1")&&TCut("HT>1250");
  resetPlotMaximum();
  if (plotsToMake.Contains("all")||plotsToMake.Contains("009")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new7_mllLoose_OF",0,"GeV");

  resetPadDimensions();
  doOverflowAddition(true);

  //some old plots with very fine binning
  nbins=100; low=0; high=200;

  //for SF-OF
  TH1D* total_sf=0;
  TH1D* total_of = 0;
  TH1D* true_edge = 0;

  selection_ = dileptons_loose && sf_loose && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1")&&TCut("HT>1500");
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  setPlotMaximum(300); //instead of Z veto and upper cut on mll
  if (plotsToMake.Contains("all")||plotsToMake.Contains("010")) {
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_mllLoose-finebins",0,"GeV");
    total_sf = (TH1D*)  totalsmsusy->Clone("total_sf");
    true_edge =  (TH1D*)getHist("naturalModel1:leptonsMatchChi2ToChi1_loose==1")->Clone("true_edge");
  }

  //now plot OF (but otherwise the same cuts)
  selection_ = dileptons_loose && !sf_loose && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1")&&TCut("HT>1500");
  resetPlotMaximum();
  if (plotsToMake.Contains("all")||plotsToMake.Contains("010")) {
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_mllLoose_OF-finebins",0,"GeV");
    total_of = (TH1D*)  totalsmsusy->Clone("total_of");
  }

  if (total_sf!=0 && total_of!=0) { //do OF subtraction
    TH1D* sf_minus_of = (TH1D*)total_sf->Clone("sf_minus_of");
    sf_minus_of->Add(total_of,-1); //subtraction

    renewCanvas();
    true_edge->Draw("hist");

    double max = findOverallMax(true_edge);
    if ( findOverallMax(sf_minus_of)>max) max = findOverallMax(sf_minus_of);
    double min = findOverallMin(true_edge);
    if ( findOverallMin(sf_minus_of)<min) min = findOverallMin(sf_minus_of);
    true_edge->SetMinimum(min);
    true_edge->SetMaximum(max);
    sf_minus_of->SetMarkerStyle(4);
    sf_minus_of->Draw("same");

    TFile fout("DelphesEdge/sfMinusOf-fine.root","recreate");
    true_edge->Write();
    total_sf->Write();
    total_of->Write();
    sf_minus_of->Write();
    thecanvas->Write();
    fout.Close();
  }

  resetPlotMaximum();
  //with tight selection, want to look at
  // max eta
  // lepton isolation
  selection_ = dileptons_loose && sf_loose && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1")&&TCut("HT>1500");
  nbins=24; low=0; high=2.4;
  var="mll_maxEta_loose"; xtitle="maximum lepton eta";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("011")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_maxEtaLoose",0,"");

  selection_ = dileptons_loose && sf_loose && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1")&&TCut("HT>1500");
  nbins=30; low=-1.5; high=1.5;
  var="leptonIso1_loose"; xtitle="RelIso lepton 1";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("011")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_leptonIso1Loose",0,"");
  var="leptonIso2_loose"; xtitle="RelIso lepton 2";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("011")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_leptonIso2Loose",0,"");

  setStackMode(false,true,false); //stack,norm,label override
  var="leptonIso1_loose"; xtitle="RelIso lepton 1";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("011")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_leptonIso1Loose",0,"");
  var="leptonIso2_loose"; xtitle="RelIso lepton 2";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("011")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_leptonIso2Loose",0,"");

}

void print_zbi() {     //let's try to spit out Zbi 
  //assumes histos are already drawn, and assumes the name of the edge signal

    double ntotal = totalsm->Integral();//contains SUSY because we see treatAllAsSM_ flag
    double ntrueedge = getIntegral("naturalModel1:leptonsMatchChi2ToChi1_loose==1");
    double nbackground = ntotal - ntrueedge;

    //assume that all background is FS, and that error on FS background is driven by sqrt(N) on OF sample
    //plus a term for a systematic, which we will take as 4%
    //first term is sqrt(N) squared
    double err = sqrt( nbackground + pow(0.04*nbackground,2));

    double zbi = jmt::zbi( ntotal, nbackground, err);
    cout<<"Zbi = "<<zbi<<endl;

}

void drawDelphes_Edge_CutAndCount(TString plotsToMake="all") {

  initSamples("tt bj nm1 skimmed");
  setOutputDirectory("DelphesEdge"); 

  treatAllAsSM_=true; //don't treat signal as signal

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false;  

  stackSignal_=false;
  setStackMode(true,false,false); //stack,norm,label override
  removeSample("naturalModel1");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose!=1&&leptonsMatchChi4ToChi1_loose!=1",kRed-9,"Other SUSY");
  addSample("naturalModel1:leptonsMatchChi4ToChi1_loose==1",kGreen+2,"True 2nd Edge");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1",kRed,"True Edge Signal");

  //edge selection cuts
  TCut dileptons_loose="mll_loose>20";
  TCut sf_loose = "isSF_loose==1";
  TCut cutncount = "mll_loose<70";

  //NOMINAL (pre-appoval) CUT AND COUNT SELECTION
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1")&&TCut("HT>1500");
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="mll";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) {
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_cutncount1",0,"");
    print_zbi();
  }

  //Same as above but using medium b-tags
 selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40eta3p0>=4") && TCut("MET>400")  && TCut("nbjets40medium>=1")&&TCut("HT>1500");
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="mll";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) {
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_cutncount1b",0,"");
    print_zbi();

  }

 //Same as above but using central jets
 selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=4") && TCut("MET>400")  && TCut("nbjets40medium>=1")&&TCut("HT>1500");
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="mll";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) {
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_cutncount1c",0,"");
    print_zbi();
  }


  //NOMINAL (pre-appoval) CUT AND COUNT SELECTION but with central jets
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1")&&TCut("HT>1500");
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="mll";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) {
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_cutncount1d",0,"");
    print_zbi();
  }


  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1")&&TCut("HT>1500") &&TCut("mll_maxEta_loose<1.4");
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="mll";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("002")) 
    {  drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new2HT1500_cutncount2",0,""); print_zbi();}

  //plot HT 
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=4") && TCut("MET>400")  && TCut("nbjets40tight>=1")&&TCut("HT>800");
  nbins=48; low=800; high=2000;
  var="HT"; xtitle="HT";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003")) 
  drawPlots(var,nbins,low,high,xtitle,"Events", "edge_cutncount_test1",0,"");

  //plot MET
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=4") && TCut("MET>200")  && TCut("nbjets40tight>=1")&&TCut("HT>1400");
  nbins=30; low=200; high=800;
  var="MET"; xtitle="MET";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003")) 
  drawPlots(var,nbins,low,high,xtitle,"Events", "edge_cutncount_test2",0,"");

  //plot njets
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=2") && TCut("MET>400")  && TCut("nbjets40tight>=1")&&TCut("HT>1500");
  nbins=8; low=2; high=10;
  var="njets40"; xtitle="njets40";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("004")) 
  drawPlots(var,nbins,low,high,xtitle,"Events", "edge_cutncount_test3",0,"");

  //nbjets after tighter njets
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=6") && TCut("MET>300")  &&TCut("HT>1200");
  nbins=8; low=0; high=8;
  var="nbjets40medium"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("008")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_cutncount_test6",0,"");

  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=6") && TCut("MET>300")  &&TCut("HT>1200");
  nbins=8; low=0; high=8;
  var="nbjets40tight"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("008")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_cutncount_test6b",0,"");

  //redo with tighter njets
  //plot HT 
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=6") && TCut("MET>300")  && TCut("nbjets40tight>=1")&&TCut("HT>800");
  nbins=48; low=800; high=2000;
  var="HT"; xtitle="HT";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("005")) 
  drawPlots(var,nbins,low,high,xtitle,"Events", "edge_cutncount_test4",0,"");

  //plot MET
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=6") && TCut("MET>200")  && TCut("nbjets40tight>=1")&&TCut("HT>1200");
  nbins=30; low=200; high=800;
  var="MET"; xtitle="MET";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("005")) 
  drawPlots(var,nbins,low,high,xtitle,"Events", "edge_cutncount_test5",0,"");

  // -- new cut and count region -- 
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=6") && TCut("MET>450")  && TCut("nbjets40tight>=1")&&TCut("HT>1250");
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="mll";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("006")) {
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new3_cutncount1",0,"");
    print_zbi();} //NEW SAMPLES -- this is good! Zbi of 4!

  // try very tight MET -- also promising
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=6") && TCut("MET>650")  && TCut("nbjets40tight>=1")&&TCut("HT>1200");
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="mll";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("008")) {
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new5_cutncount1",0,"");
    print_zbi();} 
  
  // try very tight btags
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=6") && TCut("MET>300")  && TCut("nbjets40medium>=3")&&TCut("HT>1200");
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="mll";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("009")) {
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new6_cutncount1",0,"");
    print_zbi();} 

 // so-called new7 -- i.e. post pre-approval retuning
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=6") && TCut("MET>450")  && TCut("nbjets40medium>=1")&&TCut("HT>1250");
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="mll";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("010")) {
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new7_cutncount1",0,"");
    print_zbi();} 


  if (plotsToMake.Contains("cutflowtest")) {
    savePlots_=false;
    assert( TString(gSystem->Getenv("NBJETSCUT"))!="" && TString(gSystem->Getenv("NJETSCUT"))!="");
    cout<<"Cutflow "<<gSystem->Getenv("NJETSCUT")<<" "<<gSystem->Getenv("NBJETSCUT")<<endl;
    for (int htcut = 1000; htcut<2250; htcut+=250) {
      for (int metcut = 250; metcut<700; metcut+=50) {
	TString metcutstring,htcutstring,njetsstring,nbstring;
	metcutstring.Form("MET> %d",metcut);
	htcutstring.Form("HT> %d",htcut);
	njetsstring.Form("njets40>= %d",TString(gSystem->Getenv("NJETSCUT")).Atoi() );
	nbstring.Form("nbjets40medium>= %d",TString(gSystem->Getenv("NBJETSCUT")).Atoi() );
	selection_ = dileptons_loose && sf_loose && cutncount && TCut(njetsstring.Data()) && TCut(metcutstring.Data())  && TCut(nbstring.Data())&&TCut(htcutstring.Data());
	drawPlots(var,nbins,low,high,xtitle,"Events", "blah",0,"");
	cout<<"Cutflow "<<htcut<<" "<<metcut<<" "; print_zbi();
      }
    }
    savePlots_=true;
  }

  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=6") && TCut("MET>450")  && TCut("nbjets40tight>=1")&&TCut("HT>1250")&&TCut("mll_maxEta_loose<1.4");
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="mll";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("006")) {
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new3_cutncount1-central",0,"");
    print_zbi();}

  // -- new cut and count region 2 --
  selection_ = dileptons_loose && sf_loose && cutncount && TCut("njets40>=7") && TCut("MET>350")  && TCut("nbjets40tight>=1")&&TCut("HT>1200");
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="mll";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("007")) {
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new4_cutncount1",0,"");
    print_zbi();
  }



}

void ewkino(TString plotsToMake="all") {

  initSamples("skimmed tt bj nm1");
  setOutputDirectory("DelphesEdgeEwkino");

  treatAllAsSM_=true; //don't treat signal as signal

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false;  

  stackSignal_=true;
  setStackMode(true,false,false); //stack,norm,label override
  removeSample("naturalModel1");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose!=1",kRed-9,"Other SUSY");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1",kRed,"True Edge Signal");

  //edge selection cuts
  TCut dileptons_loose="mll_loose>10";
  TCut sf_loose = "isSF_loose==1";
  TCut lowmass = "mll_loose<80"; //gotta get rid of the Z

  TCut bveto = "nbjets40loose==0";
  //  TCut tau="nTaus>0";
  TCut l3="MT_l3MET_loose>=0";
  //  TCut ht="HT<500";
  TCut loosemet = "MET>50";

  //basic trilepton selection 
  selection_ =dileptons_loose && sf_loose && lowmass && l3 &&loosemet;

  nbins=60; low=0; high=120;
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l_mll",0,"GeV");

  //plot nb
  setStackMode(false,true,false); //stack,norm,label override
  nbins=5; low=0; high=5;
  var="nbjets40loose"; xtitle="nb loose 40";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l_nb40loose",0,"");

  //plot met
  setStackMode(false,true,false); //stack,norm,label override
  nbins=50; low=0; high=500;
  var="MET"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l_met",0,"GeV");

  //plot mt
  setStackMode(false,true,false); //stack,norm,label override
  nbins=50; low=0; high=500;
  var="MT_l3MET_loose"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l_mtl3",0,"GeV");

  // lessons:
  //b veto, high met
  selection_ =dileptons_loose && sf_loose && lowmass && l3 && bveto && TCut("MET>200") && TCut("MT_l3MET_loose>150");
  setStackMode(true,false,false); //stack,norm,label override
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("002")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-2_mll",0,"GeV");

  //let in more background for some studies
  selection_ =dileptons_loose && sf_loose && lowmass && l3 && bveto && TCut("MET>150") && TCut("MT_l3MET_loose>130");
  setStackMode(false,true,false); //stack,norm,label override

  nbins=50; low=0; high=1000;
  var="HT"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-3_ht",0,"GeV");

  nbins=50; low=0; high=300;
  var="leptonPt1_loose"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-3_leppt1",0,"GeV");

  nbins=50; low=0; high=200;
  var="leptonPt2_loose"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-3_leppt2",0,"GeV");

  nbins=8; low=0; high=8;
  var="njets40eta3p0"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-3_njets",0,"GeV");

  //try complete jet veto
  selection_ =dileptons_loose && sf_loose && lowmass && l3 && TCut("njets40eta3p0==0") && TCut("MET>150") && TCut("MT_l3MET_loose>140");
  setStackMode(true,false,false); //stack,norm,label override
  nbins=60; low=0; high=120;
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("004")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-4_mll",0,"GeV");

  setStackMode(false,false,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("005")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-4_mll",0,"GeV");

  setStackMode(true,false,false); //stack,norm,label override
  removeSample("naturalModel1:leptonsMatchChi2ToChi1_loose!=1");
  removeSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1");
  addSample("naturalModel1:isSF_loose==1&&leptonsMatchChi2ToChi1_loose==0&&(Chi2ToChi1Code==1||Chi2ToChi1Code==2)",kRed-9,"SF (edge misreco)");
  addSample("naturalModel1:isSF_loose==1&&leptonsMatchChi2ToChi1_loose==0&&Chi2ToChi1Code>=10",kMagenta,"SF (N2 #rightarrow N1 with tau)");
  addSample("naturalModel1:isSF_loose==1&&leptonsMatchChi2ToChi1_loose==0&&(!(Chi2ToChi1Code==1||Chi2ToChi1Code==2||Chi2ToChi1Code>=10))",kTeal,"SF (no edge)");
  addSample("naturalModel1:isSF_loose==1&&leptonsMatchChi2ToChi1_loose==1",kRed,"SF (edge reco ok)");
  if (plotsToMake.Contains("all")||plotsToMake.Contains("006")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-4-split_mll",0,"GeV");

  nbins=7; low=0; high=7;
  var="nTrueElMu"; xtitle=var;
 if (plotsToMake.Contains("all")||plotsToMake.Contains("006")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-4-split_ntrueElMu",0,"");


 //another tack. no jet cuts except b veto, but tight in MT and MET
 // lessons:
  //b veto, high met
  selection_ =dileptons_loose && sf_loose && lowmass && l3 && bveto && TCut("MET>250") && TCut("MT_l3MET_loose>200");
  setStackMode(true,false,false); //stack,norm,label override
  nbins=60; low=0; high=120;
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("007")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-5_mll",0,"GeV");

}


void ewkino_signalonly(TString plotsToMake="all") {

  TString signalName="naturalModel1";
  initSamples(signalName);
  setOutputDirectory("DelphesSignal");

  treatAllAsSM_=true; //don't treat signal as signal

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false;  

  stackSignal_=true;
  setStackMode(true,false,false); //stack,norm,label override
  /*
  removeSample("naturalModel1");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose!=1",kRed-9,"Other SUSY");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1",kRed,"True Edge Signal");
  */

  //edge selection cuts
  TCut dileptons_loose="mll_loose>0";
  TCut sf_loose = "isSF_loose==1";
  TCut lowmass = "mll_loose<200";

  TCut bveto = "nbjets40loose==0";
  TCut l3="MT_l3MET_loose>=0";
  TCut loosemet = "MET>50";
  selection_ =dileptons_loose && sf_loose && lowmass && l3 && TCut("njets40eta3p0==0") && TCut("MET>150") && TCut("MT_l3MET_loose>140");

  clearSamples();
  addSample(signalName+":SusyProductionMode==2000000",kYellow+1,"slepton");
  addSample(signalName+":SusyProductionMode==200000",kOrange+1,"EWKino");
  addSample(signalName+":SusyProductionMode==2000",kBlue,"stop");
  addSample(signalName+":SusyProductionMode==200",kGreen,"sbottom");
  addSample(signalName+":SusyProductionMode==20000",kMagenta,"gluino");
  addSample(signalName+":SusyProductionMode==10010||SusyProductionMode==20",kRed,"other");
  setStackMode(true,false,false); //stack,norm,label override
  nbins=50; low=0; high=100;
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-4_prod_mll",0,"GeV");
  setStackMode(false,false,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-4_prod_mll",0,"GeV");


  clearSamples();
  addSample(signalName+":isSF_loose==0",kRed,"OF");
  addSample(signalName+":isSF_loose==1&&leptonsMatchChi2ToChi1_loose==0&&(Chi2ToChi1Code==1||Chi2ToChi1Code==2)",kOrange,"SF (edge misreco)");
  addSample(signalName+":isSF_loose==1&&leptonsMatchChi2ToChi1_loose==0&&Chi2ToChi1Code>=10",kMagenta,"SF (N2 #rightarrow N1 with tau)");
  addSample(signalName+":isSF_loose==1&&leptonsMatchChi2ToChi1_loose==0&&(!(Chi2ToChi1Code==1||Chi2ToChi1Code==2||Chi2ToChi1Code>=10))",kTeal,"SF (no edge)");
  addSample(signalName+":isSF_loose==1&&leptonsMatchChi2ToChi1_loose==1",kBlue,"SF (edge reco ok)");
  setStackMode(false,false,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-4_truth_mll",0,"GeV");

  //no OF for stack
  clearSamples();
  addSample(signalName+":isSF_loose==1&&leptonsMatchChi2ToChi1_loose==0&&(Chi2ToChi1Code==1||Chi2ToChi1Code==2)",kOrange,"SF (edge misreco)");
  addSample(signalName+":isSF_loose==1&&leptonsMatchChi2ToChi1_loose==0&&Chi2ToChi1Code>=10",kMagenta,"SF (N2 #rightarrow N1 with tau)");
  addSample(signalName+":isSF_loose==1&&leptonsMatchChi2ToChi1_loose==0&&(!(Chi2ToChi1Code==1||Chi2ToChi1Code==2||Chi2ToChi1Code>=10))",kTeal,"SF (no edge)");
  addSample(signalName+":isSF_loose==1&&leptonsMatchChi2ToChi1_loose==1",kBlue,"SF (edge reco ok)");
  setStackMode(true,false,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_3l-4_truth_mll",0,"GeV");

}

void preselection(TString plotsToMake="all") {
  setOutputDirectory("DelphesEdge");

  initSamples("combinesm");
  setOutputDirectory("DelphesEdge");

  treatAllAsSM_=true; //don't treat signal as signal

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false;  

  //  stackSignal_=false;

  setStackMode(false,true,false); //stack,norm,label override
  removeSample("naturalModel1");
  addSample("naturalModel1:leptonsMatchChi2ToChi1!=1",kRed-9,"Other SUSY");
  addSample("naturalModel1:leptonsMatchChi2ToChi1==1",kRed,"True Edge Signal");

  //edge selection cuts
  TCut dileptons="mll>20";
  TCut dileptons_loose="mll_loose>20";
  TCut sf = "isSF==1";
  TCut sf_loose = "isSF_loose==1";
  TCut jetsloose = "njets40eta3p0>=2";
  TCut jetstight = "njets40eta3p0>=3";
  TCut metloose = "MET>100";
  TCut mettight = "MET>150";
  //  TCut noZ = "mll<85|| mll>95";

  //use 8 TeV edge as a preselection

  selection_ = dileptons && sf && jetstight && metloose ;

  nbins=35; low=100; high=800;
  var="MET"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_presel_met",0,"GeV");
 
  nbins=100; low=0; high=3000;
  var="HT"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_presel_ht",0,"GeV");

  nbins=12; low=0; high=12;
  var="njets40eta3p0"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_presel_njets40",0,"GeV");

  nbins=6; low=0; high=6;
  var="nbjets40loose"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_presel_nbjets40loose",0,"GeV");

  nbins=6; low=0; high=6;
  var="nbjets40tight"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_presel_nbjets40tight",0,"GeV");

  nbins=40; low=0; high=200;
  var="mll"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_presel_mll",0,"GeV");

  selection_ = dileptons_loose && sf_loose && jetstight && metloose ;
  nbins=40; low=0; high=200;
  var="mll_loose"; xtitle=var;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("002")) 
    drawPlots(var,nbins,low,high,xtitle,"Events", "edge_presel_mll_loose",0,"GeV");


}


void compareOFSF(TString option="combinesm") { //
  //compare SF/OF
  useNewStyle_=true;

  setOutputDirectory("DelphesEdge");

  if (option.Contains("combinesm")) {
    initSamples("combinesm skimmed");
    
    clearSamples();
    addSample("tt-4p-0-600-v1510_14TEV:isSF_loose==1",kRed,"SF (tt)");
    addSample("tt-4p-0-600-v1510_14TEV:isSF_loose==0",kBlue,"OF (tt)");

  }
  else if (option.Contains("all")) {
    initSamples("combinesm skimmed");
    chainSamples("tt-4p-0-600-v1510_14TEV","naturalModel1");//combine SM and SUSY into one sample
    clearSamples();
    addSample("tt-4p-0-600-v1510_14TEV:isSF_loose==1",kRed ,"SF");
    addSample("tt-4p-0-600-v1510_14TEV:isSF_loose==0",kBlue,"OF");
  }
  else if (option.Contains("signal")) {
    initSamples("signal");
    
    clearSamples();
    addSample("naturalModel1:isSF_loose==1&&leptonsMatchChi2ToChi1_loose!=1",kRed,"SF (signal; no edges)");
    addSample("naturalModel1:isSF_loose==0",kBlue,"OF (signal)");
  }
  //  else if (option=="total") {
  //initSamples();
    
  //clearSamples();
    //my tricks don't work here
    //    addSample("naturalModel1:isSF==1",kRed,"SF (signal)");
    //    addSample("naturalModel1:isSF==0",kBlue,"OF (signal)");
  //}
  else assert(0);

  doRatio_=true;  
  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);

  stackSignal_=false;

  setStackMode(false,false,false); //stack,norm,label override

  //edge selection cuts
  TCut dileptons="mll_loose>20";
  TCut btag="nbjets40medium>=1";
  TCut rejectedge="leptonsMatchChi2ToChi1_loose!=1";
  //tightened selection without b-tags (remove SF cut)
  selection_ = dileptons && TCut("njets40>=6") && TCut("HT>1250") && TCut("MET>450") &&rejectedge&&btag;
  nbins=60; low=20; high=200;
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  doRatio_=true;  
  if (option.Contains("all"))  drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new7_mllloose_OFSFr_"+jmt::fortranize(option),0,"GeV");
  doRatio_=false;  
  if (option.Contains("all"))  drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new7_mllloose_OFSF_"+jmt::fortranize(option),0,"GeV");

  //plot only dilepton ttbar events
  selection_ = dileptons && TCut("njets40>=6") && TCut("HT>1250") && TCut("MET>450") &&TCut("ttbarDecayCode==11||ttbarDecayCode==12||ttbarDecayCode==13")&&btag;
  nbins=12; low=20; high=200;
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  doRatio_=true;  ratioMin=0.0; ratioMax = 2.0;
  if (option.Contains("combinesm"))  drawPlots(var,nbins,low,high,xtitle,"Events", "edge_new7-genTt2l_mllloose_OFSFr_"+jmt::fortranize(option),0,"GeV");

}

void edge_eff() {
  //plot eff of true edge signal as a function of mll
  setOutputDirectory("DelphesEdge");

  initSamples("nm1");
  
  doRatio_=false;  
  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);

  //edge selection cuts
  TCut dileptons="mll_loose>20";
 
  TCut ee = "abs(leptonFlavor1_loose)==11";
  TCut mm = "abs(leptonFlavor1_loose)==13";

  TCut ee_gencut = "abs(genLepFlavor[0])==11";
  TCut mm_gencut = "abs(genLepFlavor[0])==13";

  TCut theselection= TCut("njets40>=4") && TCut("HT>1250") && TCut("MET>350");//not as tight as my current requirements
  TCut isSF = "isSF_loose==1&&leptonsMatchChi2ToChi1_loose==1";

  savePlots_=false;
  stackSignal_=false;
  setStackMode(false,false,false); //stack,norm,label override
  nbins=14; low=20; high=160;

  //SF ee reco
  selection_ = dileptons && isSF && theselection && ee;
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0,"GeV");
  TH1D* ee_reco = (TH1D*) totalsmsusy->Clone("ee_reco");

  //mm reco
  selection_ = dileptons && isSF && theselection && mm;
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0,"GeV");
  TH1D* mm_reco = (TH1D*) totalsmsusy->Clone("mm_reco");

  //ee gen
  selection_ = TCut("genEdgeMll1>20") && theselection && ee_gencut;
  var="genEdgeMll1"; xtitle="m_{l+l-} (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0,"GeV");
  TH1D* ee_gen = (TH1D*) totalsmsusy->Clone("ee_gen");

  //mm gen
  selection_ = TCut("genEdgeMll1>20") && theselection && mm_gencut;
  var="genEdgeMll1"; xtitle="m_{l+l-} (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0,"GeV");
  TH1D* mm_gen = (TH1D*) totalsmsusy->Clone("mm_gen");

  TGraphAsymmErrors * ee_eff = new TGraphAsymmErrors();
  ee_eff->BayesDivide(ee_reco,ee_gen);

  TGraphAsymmErrors * mm_eff = new TGraphAsymmErrors();
  mm_eff->BayesDivide(mm_reco,mm_gen);

  // TH1D* r_me = (TH1D*) em_gen->Clone("r_me");
  // TH1D* R_SFOF = (TH1D*)em_gen->Clone("R_SFOF");
  // r_me->Reset();
  // R_SFOF->Reset();
  // r_me->Divide(mm_reco,ee_reco);
  // // r_mu,e = sqrt(N_mm / N_ee)
  // for (int ibin=1;ibin<=r_me->GetNbinsX();++ibin) {
  //   double err2
  //     = 0.25 * (1.0 / ee_reco->GetBinContent(ibin)) * (1.0 / mm_reco->GetBinContent(ibin)) * pow(mm_reco->GetBinError(ibin),2)
  //     + 0.25 * mm_reco->GetBinContent(ibin) * pow(ee_reco->GetBinContent(ibin),-3) * pow(ee_reco->GetBinError(ibin),2);
  //   r_me->SetBinContent(ibin, sqrt(r_me->GetBinContent(ibin)));
  //   r_me->SetBinError(ibin, sqrt(err2) );

  //   double rme=r_me->GetBinContent(ibin);
  //   R_SFOF->SetBinContent(ibin, 0.5* (rme+(1.0/rme)));
  //   R_SFOF->SetBinError(ibin, 0.5*(1- pow(rme,-2))*r_me->GetBinError(ibin));
  // }

  mm_eff->SetName("mm_eff");
  ee_eff->SetName("ee_eff");

  mm_eff->SetLineColor(kRed);
  ee_eff->SetLineColor(kBlue);

  mm_eff->SetMarkerColor(kRed);
  ee_eff->SetMarkerColor(kBlue);

  TFile fout("DelphesEdge/EdgeEff.root","recreate");
  ee_eff->Write();
  mm_eff->Write();
  ee_reco->Write();
  ee_gen->Write();
  mm_reco->Write();
  mm_gen->Write();
  // r_me->Write();
  // R_SFOF->Write();
  fout.Close();


}

void compareOFSF_eff() {

  //plot the absolute eff for 2lep ttbar events as a function of mll
  setOutputDirectory("DelphesEdge");

  initSamples("tt");
  
  clearSamples();
  addSample("tt-4p-0-600-v1510_14TEV",kRed,"tt");

  doRatio_=false;  
  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);

  //edge selection cuts
  TCut dileptons="mll_loose>20";
 
  TCut ee="ttbarDecayCode==11";
  TCut em="ttbarDecayCode==12";
  TCut mm="ttbarDecayCode==13";

  //  TCut theselection= TCut("njets40>=4") && TCut("HT>700") && TCut("MET>250");
  //  TCut theselection= TCut("njets40>=6") && TCut("HT>700") && TCut("MET>250");
  TCut theselection= TCut("njets40>=4") && TCut("HT>700") && TCut("MET>450");

  TCut isSF = "isSF_loose==1";

  savePlots_=false;

  stackSignal_=false;
  setStackMode(false,false,false); //stack,norm,label override
  //plot only dilepton ttbar events
  nbins=14; low=20; high=160;

  //SF ee reco
  selection_ = dileptons && isSF && theselection && ee;
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0,"GeV");
  TH1D* ee_reco = (TH1D*) totalsm->Clone("ee_reco");

  //em reco
  selection_ = dileptons && !isSF && theselection && em;
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0,"GeV");
  TH1D* em_reco = (TH1D*) totalsm->Clone("em_reco");

  //mm reco
  selection_ = dileptons && isSF && theselection && mm;
  var="mll_loose"; xtitle="m_{l+l-} (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0,"GeV");
  TH1D* mm_reco = (TH1D*) totalsm->Clone("mm_reco");

  //ee gen
  selection_ = TCut("ttbarGenMll>20") && theselection && ee;
  var="ttbarGenMll"; xtitle="m_{l+l-} (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0,"GeV");
  TH1D* ee_gen = (TH1D*) totalsm->Clone("ee_gen");

  //mm gen
  selection_ = TCut("ttbarGenMll>20") && theselection && mm;
  var="ttbarGenMll"; xtitle="m_{l+l-} (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0,"GeV");
  TH1D* mm_gen = (TH1D*) totalsm->Clone("mm_gen");

  //em gen
  selection_ = TCut("ttbarGenMll>20") && theselection && em;
  var="ttbarGenMll"; xtitle="m_{l+l-} (GeV)";
  drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0,"GeV");
  TH1D* em_gen = (TH1D*) totalsm->Clone("em_gen");

  TGraphAsymmErrors * ee_eff = new TGraphAsymmErrors();
  ee_eff->BayesDivide(ee_reco,ee_gen);

  TGraphAsymmErrors * em_eff = new TGraphAsymmErrors();
  em_eff->BayesDivide(em_reco,em_gen);

  TGraphAsymmErrors * mm_eff = new TGraphAsymmErrors();
  mm_eff->BayesDivide(mm_reco,mm_gen);

  TH1D* r_me = (TH1D*) em_gen->Clone("r_me");
  TH1D* R_SFOF = (TH1D*)em_gen->Clone("R_SFOF");
  r_me->Reset();
  R_SFOF->Reset();
  r_me->Divide(mm_reco,ee_reco);
  // r_mu,e = sqrt(N_mm / N_ee)
  for (int ibin=1;ibin<=r_me->GetNbinsX();++ibin) {
    double err2
      = 0.25 * (1.0 / ee_reco->GetBinContent(ibin)) * (1.0 / mm_reco->GetBinContent(ibin)) * pow(mm_reco->GetBinError(ibin),2)
      + 0.25 * mm_reco->GetBinContent(ibin) * pow(ee_reco->GetBinContent(ibin),-3) * pow(ee_reco->GetBinError(ibin),2);
    r_me->SetBinContent(ibin, sqrt(r_me->GetBinContent(ibin)));
    r_me->SetBinError(ibin, sqrt(err2) );

    double rme=r_me->GetBinContent(ibin);
    R_SFOF->SetBinContent(ibin, 0.5* (rme+(1.0/rme)));
    R_SFOF->SetBinError(ibin, 0.5*(1- pow(rme,-2))*r_me->GetBinError(ibin));
  }

  mm_eff->SetName("mm_eff");
  em_eff->SetName("em_eff");
  ee_eff->SetName("ee_eff");

  mm_eff->SetLineColor(kRed);
  ee_eff->SetLineColor(kBlue);
  em_eff->SetLineColor(kBlack);

  mm_eff->SetMarkerColor(kRed);
  ee_eff->SetMarkerColor(kBlue);
  em_eff->SetMarkerColor(kBlack);

  TFile fout("DelphesEdge/ttbar2leff.root","recreate");
  r_me->Write();
  R_SFOF->Write();
  ee_eff->Write();
  em_eff->Write();
  mm_eff->Write();
  ee_reco->Write();
  ee_gen->Write();
  em_reco->Write();
  em_gen->Write();
  mm_reco->Write();
  mm_gen->Write();
  fout.Close();

}

void drawSignalByProductionMode() {
 lumiScale_=3000e3;

  initSamples("nm1 skimmed");
  setOutputDirectory("DelphesEdge");
  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(false);
  doRatio_=false;  
  treatAllAsSM_=true; //don't treat signal as signal

  setStackMode(true,false,false); //stack,norm,label override

  //edge selection cuts
  TCut dileptons="mll_loose>20";
  TCut sf = "isSF_loose==1";
  TCut jets = "njets40eta3p0>=4";
  TCut met = "MET>400";
  TCut ht = "HT>1500";  

  TCut btags = "nbjets40tight>=1";
  TCut bveto = "nbjets40tight==0";

  // does not include sf cut
  TCut analysis_selection = dileptons && jets && met && ht && btags && sf;
  nbins = 40; low = 20; high=100;
  var = "mll_loose"; xtitle="mll";

  clearSamples();
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose!=1",kRed-9,"Non-edge SUSY");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1&&(SusyProductionMode==2000000||SusyProductionMode==200000)",kYellow+1,"l~+EWKino edge");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1&&(SusyProductionMode==2000||SusyProductionMode==200)",kBlue,"t~t~+b~b~ edge");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1&&SusyProductionMode==20000",kMagenta,"g~g~ edge");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1&&SusyProductionMode==10010",kCyan,"g~q~ edge");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1&&SusyProductionMode==20",kGreen,"q~q~ edge");
  addSample("naturalModel1:leptonsMatchChi2ToChi1_loose==1&&!(SusyProductionMode==2000000||SusyProductionMode==200000||SusyProductionMode==2000||SusyProductionMode==200||SusyProductionMode==20000||SusyProductionMode==10010||SusyProductionMode==20)",kRed,"other edge");
  

  selection_ = analysis_selection;
  drawPlots(var,nbins,low,high,xtitle,"Events", "edge_signal_by_prod",0,"GeV");
 
}

void generateHistosForFit(bool vetoChi4=false) {
  lumiScale_=3000e3;

  initSamples("tt bj nm1 skimmed");
  setOutputDirectory("DelphesEdge");
  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false;  

  if (vetoChi4) {
    removeSample("naturalModel1");
    addSample("naturalModel1:leptonsMatchChi4ToChi1_loose!=1",kRed,"NM1 (no Chi4 Edge)");
  }

  setStackMode(true,false,false); //stack,norm,label override
  stackSignal_=false; //the logic of this is a bit broken at the moment, so just don't stack it

  //edge selection cuts
  TCut dileptons="mll_loose>20";
  TCut sf = "isSF_loose==1";
  TCut jets = "njets40>=6";
  TCut met = "MET>450";
  TCut ht = "HT>1250";  

  TCut btags = "nbjets40medium>=1";
  TCut bveto = "nbjets40medium==0";

  TCut ee = "abs(leptonFlavor1_loose)==11";
  TCut mm = "abs(leptonFlavor1_loose)==13";

  // does not include sf cut
  TCut analysis_selection = dileptons && jets && met && ht && btags;
  //could split by e and mu
  TCut dy_selection = dileptons && sf && jets && ht && TCut("MET<40") && bveto; 

  nbins = 280; low = 20; high=300;
  var = "mll_loose"; xtitle="mll";

  selection_ = analysis_selection && sf && ee;
  drawPlots(var,nbins,low,high,xtitle,"Events", "templates_new7_mll_ee",0,"GeV");
  TH1D* smsusy_mll_ee = (TH1D*)  totalsmsusy->Clone("smsusy_mll_ee");
  TH1D* sm_mll_ee = (TH1D*)  totalsm->Clone("sm_mll_ee");

  selection_ = analysis_selection && sf && mm;
  drawPlots(var,nbins,low,high,xtitle,"Events", "templates_new7_mll_mm",0,"GeV");
  TH1D* smsusy_mll_mm = (TH1D*)  totalsmsusy->Clone("smsusy_mll_mm");
  TH1D* sm_mll_mm = (TH1D*)  totalsm->Clone("sm_mll_mm");

  selection_ = analysis_selection && !sf;
  drawPlots(var,nbins,low,high,xtitle,"Events", "templates_new7_mll_OF",0,"GeV");
  TH1D* smsusy_mll_OF = (TH1D*)  totalsmsusy->Clone("smsusy_mll_OF");
  TH1D* sm_mll_OF = (TH1D*)  totalsm->Clone("sm_mll_OF");

  selection_ = dy_selection;
  nbins = 280*5; low = 20; high=300;
  drawPlots(var,nbins,low,high,xtitle,"Events", "templates_new7_mll_DY",0,"GeV");
  TH1D* smsusy_mll_DY = (TH1D*)  totalsmsusy->Clone("smsusy_mll_DY");
  TH1D* sm_mll_DY = (TH1D*)  totalsm->Clone("sm_mll_DY");

  // new histograms: include sm and non-edge susy
  //i.e. the "background only" hypothesis when fitting for edge signal

  
  if (vetoChi4) {
    removeSample("naturalModel1:leptonsMatchChi4ToChi1_loose!=1");
    addSample("naturalModel1:leptonsMatchChi4ToChi1_loose!=1&&leptonsMatchChi2ToChi1_loose!=1",kRed-9,"No Chi4 or Chi2");
  }
  else {
    removeSample("naturalModel1");
    addSample("naturalModel1:leptonsMatchChi2ToChi1_loose!=1",kRed-9,"Other SUSY");
  }
  nbins = 280; low = 20; high=300;
  selection_ = analysis_selection && sf && ee;
  drawPlots(var,nbins,low,high,xtitle,"Events", "templates_new7_mll_ee_noedge",0,"GeV");
  TH1D* smsusynoedge_mll_ee = (TH1D*)  totalsmsusy->Clone("smsusynoedge_mll_ee");

  selection_ = analysis_selection && sf && mm;
  drawPlots(var,nbins,low,high,xtitle,"Events", "templates_new7_mll_mm_noedge",0,"GeV");
  TH1D* smsusynoedge_mll_mm = (TH1D*)  totalsmsusy->Clone("smsusynoedge_mll_mm");

  selection_ = analysis_selection && !sf ;
  drawPlots(var,nbins,low,high,xtitle,"Events", "templates_new7_mll_OF_noedge",0,"GeV");
  TH1D* smsusynoedge_mll_OF = (TH1D*)  totalsmsusy->Clone("smsusynoedge_mll_OF");
 
  TString outfile = vetoChi4 ? "DelphesEdge/templates_noChi4.root" : "DelphesEdge/templates.root";
  TFile fout(outfile,"recreate");
  smsusy_mll_ee->Write();
  smsusy_mll_mm->Write();
  smsusy_mll_OF->Write();
  smsusy_mll_DY->Write();
  smsusynoedge_mll_ee->Write();
  smsusynoedge_mll_mm->Write();
  smsusynoedge_mll_OF->Write();
  sm_mll_ee->Write();
  sm_mll_mm->Write();
  sm_mll_OF->Write();
  sm_mll_DY->Write();
  fout.Close();

}


//copy/paste from MT2 code, then modified for edge
void make_cutflowtable_fancy() {

  TString samples="bj tt nm1";
  initSamples(samples);
  setOutputDirectory("DelphesEdge");
  lumiScale_=3000e3; //request to use 3000

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false;  

  //  stackSignal_=false;
  setQuiet(true);
  stackSignal_=false;
  setStackMode(true,false,false); //stack,norm,label override

  removeSample("naturalModel1");
  TString othersignal = "naturalModel1:leptonsMatchChi2ToChi1_loose!=1&&leptonsMatchChi4ToChi1_loose!=1";
  TString chi4signal ="naturalModel1:leptonsMatchChi4ToChi1_loose==1";
  TString chi2signal="naturalModel1:leptonsMatchChi2ToChi1_loose==1";
  addSample(othersignal,kRed-9,"Other NM1");
  addSample(chi4signal,kGreen+2,"#tilde{#chi}_{4}^{0} #rightarrow #tilde{l}l #rightarrow l^{+}l^{-} #tilde{#chi}_{1}^{0}");
  addSample(chi2signal,kRed,"#tilde{#chi}_{2}^{0} #rightarrow #tilde{l}l #rightarrow l^{+}l^{-} #tilde{#chi}_{1}^{0}");

  setPadDimensions(800,600);

  //  selection_ = dileptons_loose && sf_loose && TCut("njets40>=6") && TCut("MET>450")  && TCut("nbjets40medium>=1")&&TCut("HT>1250"); //updated NOMINAL SELECTION

  TCut dileptons_loose = "mll_loose>20 && mll_loose<70";
  TCut sf_loose = "isSF_loose==1";

  TCut btags="nbjets40medium>=1";
  TCut jets="njets40>=6";

  TCut htcut = "HT>1250";
  TCut metcut="MET>450";

  std::vector<TCut> cut_list;
  cut_list.push_back(jets );
  cut_list.push_back(jets && htcut);
  cut_list.push_back(jets && htcut &&metcut);
  cut_list.push_back(jets && htcut &&metcut &&btags);
  cut_list.push_back(jets && htcut &&metcut &&btags && dileptons_loose && sf_loose);

  const  double dB=0.5;

  savePlots_=false;
  cout<<"Fractional error on background assumed to be "<<dB<<endl;
  for (int k=(int)samples_.size()-1;k>=0;k--) if (isSampleSM(samples_.at(k))) cout<<" & "; //one divider per SM sample
  cout<<"   & \\multicolumn{c}{S/\\sqrt{S+B+\\deltaB^2}} & \\multicolumn{c}{Zbi} \\\\"<<endl;
  for (int k=(int)samples_.size()-1;k>=0;k--) if (isSampleSM(samples_.at(k))) cout<<samples_.at(k)<<" & "; //one divider per SM sample
  cout<<"  SM & Other & Chi4 & Chi2 & Other & Chi4 & Chi2 & Other & Chi4 & Chi2 \\\\"<<endl;
  for (unsigned int k=0; k<cut_list.size() ; k++) {
    selection_ = cut_list.at(k);
    cout<<selection_<<" & ";

    nbins=1; low=0; high=1e9;
    var="MT2"; xtitle="MT2 (GeV)";
    drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0,"");
    for (int isample=(int)samples_.size()-1;isample>=0;isample--) {
      if (isSampleSM(samples_.at(isample))) {
	TString o;
	o.Form("%.2f &",getIntegral(samples_.at(isample)));
	cout<<o;
      }
    }
    double B = getIntegral("totalsm");
    double S1 = getIntegral(othersignal);
    double S2 = getIntegral(chi4signal);
    double S3 = getIntegral(chi2signal);

    TString output;
    output.Form(" %.2f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f &%.1f & %.1f \\\\",
		///		selection_.GetTitle(),
		B,
		S1,S2,S3,
		//S1/sqrt(B),S2/sqrt(B),S3/sqrt(B),
		S1/sqrt(S1 + B + dB*dB*B*B),S2/sqrt(S2 + B + dB*dB*B*B),S3/sqrt(S3 + B + dB*dB*B*B) ,
		jmt::zbi(S1+B,B,dB*B),jmt::zbi(S2+B,B,dB*B),jmt::zbi(S3+B,B,dB*B)
		);
    cout<<output<<endl;
  }


}
