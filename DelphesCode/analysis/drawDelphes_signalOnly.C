#include "drawDelphesBase.C"
/*
--- to compile ---

.L TSelectorMultiDraw.C+
.L CrossSectionTable.cxx+
.L ConfigurationDescriptions.cxx+

.L drawDelphes_signalOnly.C+

*/

void drawDelphes_compareSignals(TString plotsToMake="all") {
  initSamples("nm1 nm2 nm3 stoc");
  setOutputDirectory("DelphesSignal");

  treatAllAsSM_=true; //don't treat signal as signal

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false;  

  stackSignal_=false;

  setStackMode(false,false,false); //stack,norm,label override
  selection_ = "(1)"; //no cuts!
  nbins=8; low=1; high=9;
  //slepton, ewkino, stop, sbottom, gluino, gluino+squark, sq/sq
  setPlotMaximum(20000);
  var="(SusyProductionMode>=2000000)*1+(SusyProductionMode>=200000&&SusyProductionMode<2000000)*2+(SusyProductionMode==2000)*3+(SusyProductionMode==200)*4+(SusyProductionMode==20000)*5+(SusyProductionMode==10010)*6+(SusyProductionMode==20)*7"; xtitle="Susy Production Mode";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001"))  drawPlots(var,nbins,low,high,xtitle,"Events", "compareSignals_prod",0,"");
  resetPlotMaximum();

  selection_ = "SusyProductionMode==20000"; //require gluino-gluino production
  nbins=5; low=0; high=5;
  var="nElectrons+nMuons"; xtitle="n e+#mu";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001"))  drawPlots(var,nbins,low,high,xtitle,"Events", "compareSignals_nLep",0,"");
  selection_ = "SusyProductionMode==20000"; //require gluino-gluino production
  nbins=16; low=0; high=16;
  var="njets40"; xtitle="Jet Multiplicity (40 GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001"))  drawPlots(var,nbins,low,high,xtitle,"Events", "compareSignals_njets40",0,"");

  selection_ = "SusyProductionMode==20000 && nElectrons==0 && nMuons==0"; //require gluino-gluino production -- and no leptons
  nbins=100; low=0; high=3000;
  var="HT"; xtitle="HT (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001"))  drawPlots(var,nbins,low,high,xtitle,"Events", "compareSignals_HT",0,"");

  selection_ = "SusyProductionMode==20000"; //require gluino-gluino production
  nbins=7; low=0; high=7;
  var="nbjets40medium"; xtitle="Number of medium b tags";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001"))  drawPlots(var,nbins,low,high,xtitle,"Events", "compareSignals_btagsloose",0,"");

  selection_ = "SusyProductionMode==20000"; //require gluino-gluino production
  nbins=100; low=0; high=1000;
  var="MT2"; xtitle="MT2 (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001"))  drawPlots(var,nbins,low,high,xtitle,"Events", "compareSignals_MT2",0,"");

  setStackMode(false,true,false); //stack,norm,label override
  TCut cleanup="VSPT<175 && minDeltaPhi>0.3";
  TCut noleptons = "nElectrons==0 && nMuons==0";
  TCut btags="nbjets40medium>=3";
  TCut tighterjets="njets40>=7";
  TCut tightht="HT>1500";
  selection_ = cleanup&&noleptons&&btags&&tighterjets&&tightht && TCut("MT2>550");
  nbins=8; low=1; high=9;
  //slepton, ewkino, stop, sbottom, gluino, gluino+squark, sq/sq
  var="(SusyProductionMode>=2000000)*1+(SusyProductionMode>=200000&&SusyProductionMode<2000000)*2+(SusyProductionMode==2000)*3+(SusyProductionMode==200)*4+(SusyProductionMode==20000)*5+(SusyProductionMode==10010)*6+(SusyProductionMode==20)*7"; xtitle="Susy Production Mode";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("006"))  drawPlots(var,nbins,low,high,xtitle,"Events", "compareSignals_prod_tight_hadronic_cuts",0,"");

}

void  drawDelphes_signalOnly(TString plotsToMake="all",TString signalToDraw="nm1") {

  TString signalName="naturalModel";
  if (signalToDraw=="nm1") signalName+="1";
  else if (signalToDraw=="nm2") signalName+="2";
  else if (signalToDraw=="nm3") signalName+="3";
  else assert(0);

  initSamples(signalToDraw);
  setOutputDirectory("DelphesSignal");

  treatAllAsSM_=true; //don't treat signal as signal

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false;  

  stackSignal_=false;

  setStackMode(false,false,false); //stack,norm,label override
  selection_ = "(1)"; //no cuts!
  nbins=9; low=1; high=10;
  //slepton, ewkino, stop, sbottom, gluino, gluino+squark, sq/sq, other
  var="(SusyProductionMode>=2000000)*1+(SusyProductionMode>=200000&&SusyProductionMode<2000000)*2+(SusyProductionMode==2000)*3+(SusyProductionMode==200)*4+(SusyProductionMode==20000)*5+(SusyProductionMode==10010)*6+(SusyProductionMode==20)*7"; xtitle="Susy Production Mode";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001"))  drawPlots(var,nbins,low,high,xtitle,"Events", signalToDraw+"_prod",0,"");
  setPlotMaximum(12e3);
  if (plotsToMake.Contains("all")||plotsToMake.Contains("001"))  drawPlots(var,nbins,low,high,xtitle,"Events", signalToDraw+"_prod_zoom",0,"");
  resetPlotMaximum();

  setStackMode(false,true,false); //stack,norm,label override
  //split by SUSY mode
  clearSamples();
  addSample(signalName+":SusyProductionMode>=2000000",kYellow+1,"slepton");
  addSample(signalName+":SusyProductionMode>=200000&&SusyProductionMode<2000000",kOrange+1,"EWKino");
  addSample(signalName+":SusyProductionMode==2000",kBlue,"stop");
  addSample(signalName+":SusyProductionMode==200",kGreen,"sbottom");
  addSample(signalName+":SusyProductionMode==20000",kMagenta,"gluino");
  addSample(signalName+":SusyProductionMode==10010||SusyProductionMode==20",kRed,"other");

  selection_ = "(1)"; //no cuts!
  nbins=50; low=0; high=800;
  var="MET"; xtitle="MET";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("002"))  drawPlots(var,nbins,low,high,xtitle,"Events", signalToDraw+"_byMode_MET",0,"");

  nbins=50; low=0; high=3000;
  var="HT"; xtitle="HT";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("002"))  drawPlots(var,nbins,low,high,xtitle,"Events", signalToDraw+"_byMode_HT",0,"");

  nbins=12; low=0; high=12;
  var="njets30eta3p0"; xtitle="njets30eta3p0";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("002"))  drawPlots(var,nbins,low,high,xtitle,"Events", signalToDraw+"_byMode_njets30eta3p0",0,"");

  nbins=6; low=0; high=6;
  var="nbjets40tight"; xtitle="nbjets40tight";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("002"))  drawPlots(var,nbins,low,high,xtitle,"Events", signalToDraw+"_byMode_nbjets40tight",0,"");

  selection_ = "mll>0 && isSF==1"; //loosest possible dilepton selection
  nbins=40; low=0; high=200;
  var="mll"; xtitle="m_{l+l-}";

  setStackMode(true,false,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("002"))  drawPlots(var,nbins,low,high,xtitle,"Events", signalToDraw+"_byMode_mll",0,"");
  setStackMode(false,true,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("002"))  drawPlots(var,nbins,low,high,xtitle,"Events", signalToDraw+"_byMode_mll",0,"");


  //split instead by Chi2ToChi1Code
 setStackMode(false,true,false); //stack,norm,label override
  //split by SUSY mode
  clearSamples();
  addSample(signalName+":Chi2ToChi1Code==0",kYellow+1,"No N2 #rightarrow N1");
  addSample(signalName+":Chi2ToChi1Code==1",kRed,"N2 #rightarrow N1 via e/#mu");
  addSample(signalName+":Chi2ToChi1Code==2",kOrange,"N2 #rightarrow N1 via e/#mu x 2");
  addSample(signalName+":Chi2ToChi1Code>=10",kBlue,"N2 #rightarrow N1 via #tau");

  selection_ = "(1)"; //no cuts!
  nbins=8; low=1; high=9;
  //slepton, ewkino, stop, sbottom, gluino, gluino+squark, sq/sq
  var="(SusyProductionMode>=2000000)*1+(SusyProductionMode>=200000&&SusyProductionMode<2000000)*2+(SusyProductionMode==2000)*3+(SusyProductionMode==200)*4+(SusyProductionMode==20000)*5+(SusyProductionMode==10010)*6+(SusyProductionMode==20)*7"; xtitle="Susy Production Mode";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003"))  drawPlots(var,nbins,low,high,xtitle,"Events", signalToDraw+"_prod_byN2toN1",0,"");

  setStackMode(false,false,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003"))  drawPlots(var,nbins,low,high,xtitle,"Events", signalToDraw+"_prod_byN2toN1",0,"");

  selection_ = "mll>0 && isSF==1"; //loosest possible dilepton selection
  nbins=40; low=0; high=200;
  var="mll"; xtitle="m_{l+l-}";

  setStackMode(false,false,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003"))  drawPlots(var,nbins,low,high,xtitle,"Events", signalToDraw+"_byN2toN1_mll",0,"");
  setStackMode(false,true,false); //stack,norm,label override
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003"))  drawPlots(var,nbins,low,high,xtitle,"Events", signalToDraw+"_byN2toN1_mll",0,"");

  //now, forget gen-level (or not completely)
  clearSamples();
  addSample(signalName+":isSF==0",kRed,"OF");
  addSample(signalName+":isSF==1&&leptonsMatchChi2ToChi1==0&&(Chi2ToChi1Code==1||Chi2ToChi1Code==2)",kOrange,"SF (edge misreco)");
  addSample(signalName+":isSF==1&&leptonsMatchChi2ToChi1==0&&Chi2ToChi1Code>=10",kMagenta,"SF (N2 #rightarrow N1 with tau)");
  addSample(signalName+":isSF==1&&leptonsMatchChi2ToChi1==0&&(!(Chi2ToChi1Code==1||Chi2ToChi1Code==2||Chi2ToChi1Code>=10))",kTeal,"SF (no edge)");
  //  addSample(signalName+":isSF==1&&leptonsMatchChi2ToChi1==0&&nZFromSusy==0",kTeal,"SF (no reco match; no Z)");
  //  addSample(signalName+":isSF==1&&leptonsMatchChi2ToChi1==0&&nZFromSusy>0",kOrange,"SF (no reco match; with Z)");
  addSample(signalName+":isSF==1&&leptonsMatchChi2ToChi1==1",kBlue,"SF (edge reco ok)");

  //leptonsMatchChi2ToChi1

  setStackMode(false,false,false); //stack,norm,label override 
  selection_ = "mll>0"; //loosest possible dilepton selection
  nbins=40; low=0; high=200;
  var="mll"; xtitle="m_{l+l-}";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("004"))  drawPlots(var,nbins,low,high,xtitle,"Events", signalToDraw+"_bySF_mll",0,"");


  clearSamples(); //now just plot SF part for stack
  addSample(signalName+":isSF==1&&leptonsMatchChi2ToChi1==0&&(Chi2ToChi1Code==1||Chi2ToChi1Code==2)",kOrange,"SF (edge misreco)");
  addSample(signalName+":isSF==1&&leptonsMatchChi2ToChi1==0&&Chi2ToChi1Code>=10",kMagenta,"SF (N2 #rightarrow N1 with tau)");
  addSample(signalName+":isSF==1&&leptonsMatchChi2ToChi1==0&&(!(Chi2ToChi1Code==1||Chi2ToChi1Code==2||Chi2ToChi1Code>=10))",kTeal,"SF (no edge)");
  //  addSample(signalName+":isSF==1&&leptonsMatchChi2ToChi1==0&&nZFromSusy==0",kTeal,"SF (no reco match; no Z)");
  //  addSample(signalName+":isSF==1&&leptonsMatchChi2ToChi1==0&&nZFromSusy>0",kOrange,"SF (no reco match; with Z)");
  addSample(signalName+":isSF==1&&leptonsMatchChi2ToChi1==1",kBlue,"SF (edge reco ok)");
  setStackMode(true,false,false); //stack,norm,label override 
  if (plotsToMake.Contains("all")||plotsToMake.Contains("004"))  drawPlots(var,nbins,low,high,xtitle,"Events", signalToDraw+"_bySF_mll",0,"");

  //finally, let's go for the actual cut of the Edge analysis
  clearSamples();
  addSample(signalName+":mll_maxEta<1.4",kRed,signalName+" (central)");
  addSample(signalName+":mll_maxEta>1.4",kBlue,signalName+" (forward)");

  setStackMode(false,false,false); //stack,norm,label override 
  nbins=40; low=0; high=200;
  var="mll"; xtitle="m_{l+l-}";

  selection_="mll>20&&isSF==1 && njets40eta3p0>=2 && MET>150"; //one search region
  if (plotsToMake.Contains("all")||plotsToMake.Contains("005"))  drawPlots(var,nbins,low,high,xtitle,"Events", signalToDraw+"_sel1_mll",0,"");

  selection_="mll>20&&isSF==1 && njets40eta3p0>=3 && MET>100"; //the other search region
  if (plotsToMake.Contains("all")||plotsToMake.Contains("005"))  drawPlots(var,nbins,low,high,xtitle,"Events", signalToDraw+"_sel2_mll",0,"");

  //do it again, but split by SUSY production mode and stack ; integrate over maxEta
  clearSamples();
  addSample(signalName+":SusyProductionMode>=2000000",kYellow+1,"slepton");
  addSample(signalName+":SusyProductionMode>=200000&&SusyProductionMode<2000000",kOrange+1,"EWKino");
  addSample(signalName+":SusyProductionMode==2000",kBlue,"stop");
  addSample(signalName+":SusyProductionMode==200",kGreen,"sbottom");
  addSample(signalName+":SusyProductionMode==20000",kMagenta,"gluino");
  addSample(signalName+":SusyProductionMode==10010||SusyProductionMode==20",kRed,"other");

  setStackMode(true,false,false); //stack,norm,label override 
  nbins=40; low=0; high=200;
  var="mll"; xtitle="m_{l+l-}";

  selection_="mll>20&&isSF==1 && njets40eta3p0>=2 && MET>150"; //one search region
  if (plotsToMake.Contains("all")||plotsToMake.Contains("005"))  drawPlots(var,nbins,low,high,xtitle,"Events", signalToDraw+"_sel1_byMode_mll",0,"");

  selection_="mll>20&&isSF==1 && njets40eta3p0>=3 && MET>100"; //the other search region
  if (plotsToMake.Contains("all")||plotsToMake.Contains("005"))  drawPlots(var,nbins,low,high,xtitle,"Events", signalToDraw+"_sel2_byMode_mll",0,"");


  setStackMode(false,true,false); //stack,norm,label override
  TCut cleanup="VSPT<175 && minDeltaPhi>0.3";
  TCut noleptons = "nElectrons==0 && nMuons==0";
  TCut btags="nbjets40medium>=3";
  TCut tighterjets="njets40>=7";
  TCut tightht="HT>1500";
  selection_ = cleanup&&noleptons&&btags&&tighterjets&&tightht && TCut("MT2>550");
  nbins=8; low=1; high=9;
  //slepton, ewkino, stop, sbottom, gluino, gluino+squark, sq/sq
  var="(SusyProductionMode>=2000000)*1+(SusyProductionMode>=200000&&SusyProductionMode<2000000)*2+(SusyProductionMode==2000)*3+(SusyProductionMode==200)*4+(SusyProductionMode==20000)*5+(SusyProductionMode==10010)*6+(SusyProductionMode==20)*7"; xtitle="Susy Production Mode";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("006"))  drawPlots(var,nbins,low,high,xtitle,"Events", signalToDraw+"_prod_tight_hadronic_cuts",0,"");

}
