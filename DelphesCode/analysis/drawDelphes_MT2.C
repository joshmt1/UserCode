#include "drawDelphesBase.C"
/*
--- to compile ---

.L TSelectorMultiDraw.C+
.L CrossSectionTable.cxx+
.L ConfigurationDescriptions.cxx+

.L drawDelphes_MT2.C+

*/

void drawVSPT() {
  //to study VSPT shape in ttbar
  initSamples("tt");
  setOutputDirectory("DelphesMT2");
  lumiScale_=300e3; //discovery point


  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false;  

  //  stackSignal_=false;

  stackSignal_=false;
  setStackMode(false,true,false); //stack,norm,label override

  selection_ = "ttbarDecayCode==2"; //hadronic ttbar as a proxy for qcd

  clearSamples();
  addSample("tt-4p-0-600-v1510_14TEV:VSPT<35",kBlack,"0<VSPT<35");
  addSample("tt-4p-0-600-v1510_14TEV:VSPT>=35&&VSPT<70",kRed,"35<VSPT<70");
  addSample("tt-4p-0-600-v1510_14TEV:VSPT>=70&&VSPT<105",kGreen,"70<VSPT<105");
  addSample("tt-4p-0-600-v1510_14TEV:VSPT>=105&&VSPT<140",kBlue,"105<VSPT<140");
  addSample("tt-4p-0-600-v1510_14TEV:VSPT>=140&&VSPT<175",kMagenta,"140<VSPT<175");
  addSample("tt-4p-0-600-v1510_14TEV:VSPT>=175",kCyan,"VSPT>175");

  setLogY(true);
  setPlotMinimum(1e-5);
 //nbjets
  nbins=100; low=0; high=500;
  var="MT2"; xtitle=var;
  drawPlots(var,nbins,low,high,xtitle,"Events", "hadttbar_VSPTstudy",0,"");

}

void  drawDelphes_MT2(TString plotsToMake="all") {
  useNewStyle_=true;  

  initSamples("nm1 nm2 nm3 bj tt htskim");
  setOutputDirectory("DelphesMT2");
  lumiScale_=3000e3; //request is to use 3000fb-1
  //  treatAllAsSM_=true; //don't treat signal as signal

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false;  

  //  stackSignal_=false;

  stackSignal_=false;
  setStackMode(true,false,false); //stack,norm,label override

  TCut cleanup="VSPT<175 && minDeltaPhi>0.3";

  //tighter further!

  TCut noleptons = "nElectrons==0 && nMuons==0";
  TCut tightermt2 = "MT2>500";
  TCut tighterjets = "njets40>=8";
  TCut ht2000="HT>2000";

  //need to investigate the VSPT and minDeltaPhi distributions
  selection_ = noleptons && ht2000 && tightermt2 && tighterjets ;
  nbins=50; low=0; high=500;
  var="VSPT"; xtitle="VSPT (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("002"))   drawPlots(var,nbins,low,high,xtitle,"Events", "hadronic_mt2tightNoCleanup_VSPT",0,"");
  nbins=50; low=0; high=3.5;
  var="minDeltaPhi"; xtitle="minDeltaPhi";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("002"))   drawPlots(var,nbins,low,high,xtitle,"Events", "hadronic_mt2tightNoCleanup_minDeltaPhi",0,"");

  //now the full selection
  setPadDimensions(800,600);
  isPreliminary_=true;
  selection_ = cleanup && noleptons && ht2000 && tightermt2 && tighterjets ;
 //nbjets
  nbins=8; low=0; high=8;
  var="nbjets40medium"; xtitle="b tag multiplicity";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003"))   drawPlots(var,nbins,low,high,xtitle,"Events", "hadronic_mt2tight1_nbjets",0,"");

  relPosX+=0.005;
  TCut loosermt2 = "MT2>=200";
  TCut btags = "nbjets40medium>=3";
  selection_ = cleanup && noleptons && ht2000 && loosermt2 && tighterjets&&btags;
  nbins=24; low=200; high=800;
  var="MT2"; xtitle="M_{T2} (GeV)";
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003"))   drawPlots(var,nbins,low,high,xtitle,"Events", "hadronic_mt2-3b_mt2",0,"GeV");

  //request: remake this plot with longer MT2 range
  nbins=20; low=200; high=1200;
  if (plotsToMake.Contains("all")||plotsToMake.Contains("003"))   drawPlots(var,nbins,low,high,xtitle,"Events", "hadronic_mt2-3b_mt2_wide",0,"GeV");
}

void make_cutflowtable(bool useTauVeto,float htCutValue,int nbtags,int njets=8) {


  initSamples("all htskim");
  setOutputDirectory("DelphesMT2");
  lumiScale_=300e3; //discovery point
  //  treatAllAsSM_=true; //don't treat signal as signal

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false;  

  //  stackSignal_=false;
  setQuiet(true);
  stackSignal_=false;
  setStackMode(true,false,false); //stack,norm,label override

  TCut cleanup="VSPT<175 && minDeltaPhi>0.3";

  //tighter further!
  TCut noleptons = "nElectrons==0 && nMuons==0";
  TCut tauveto = "nTaus==0";

  TString nbstring;
  nbstring.Form("nbjets40medium>=%d",nbtags);
  TCut btags=nbstring.Data();

  TString njstring;
  njstring.Form("njets40>=%d",njets);
  TCut tighterjets=njstring.Data();

  TCut baseline=cleanup&&noleptons&&tighterjets&&btags;

  if (useTauVeto) baseline = baseline&&tauveto;

  TString htstring;
  htstring.Form("HT>%.0f",htCutValue);
  baseline = baseline && TCut( htstring.Data());


  /*
scan over mt2 cut values; compute S/sqrt(B) or S/sqrt(S+B+ dB*dB)
try:
with and without tau veto
for:
HT>1000
HT>1500
HT>2000
HT>2500
  */

  const  double dB=0.5;

  savePlots_=false;
  cout<<baseline.GetTitle()<<endl;
  cout<<"Fractional error on background assumed to be "<<dB<<endl;
  //  cout<<"     &    & \\multicolumn{c}{S/\\sqrt{B++\\deltaB^2}} & \\multicolumn{c}{S/\\sqrt{S+B+\\deltaB^2}} \\\\"<<endl;
  cout<<"     & SM & NM1 & NM2 & NM3 & NM1 & NM2 & NM3  \\\\"<<endl;
  for (float mt2CutValue=350; mt2CutValue<=650; mt2CutValue+=50) {
    TString mt2string;
    mt2string.Form("MT2>%.0f",mt2CutValue);
    selection_ = baseline && TCut( mt2string.Data());

    nbins=1; low=0; high=1e9;
    var="MT2"; xtitle="MT2 (GeV)";
    drawPlots(var,nbins,low,high,xtitle,"Events", "dummy",0,"");
    double B = getIntegral("totalsm");
    double S1 = getIntegral("naturalModel1");
    double S2 = getIntegral("naturalModel2");
    double S3 = getIntegral("naturalModel3");

    TString output;
    output.Form(" MT2 > %.0f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f  \\\\",
		mt2CutValue,
		B,
		S1,S2,S3,
		/*		S1/sqrt(B),S2/sqrt(B),S3/sqrt(B),*/
		//		S1/sqrt( B + dB*dB*B*B),S2/sqrt( B + dB*dB*B*B),S3/sqrt( B + dB*dB*B*B),
		jmt::zbi(S1+B,B,dB*B),jmt::zbi(S2+B,B,dB*B),jmt::zbi(S3+B,B,dB*B) );
    cout<<output<<endl;
  }


}



void make_cutflowtable_fancy(bool useTauVeto=false,float htCutValue=2000,int nbtags=3,int njets=8,bool PhaseII=true) {

  TString samples="bj tt nm1 nm2 nm3";
  if (!PhaseII) samples+= " PhaseI_50";
  initSamples(samples);
  //initSamples("signal");
  setOutputDirectory("DelphesMT2");
  lumiScale_=3000e3; //request to use 3000
  if (!PhaseII) lumiScale_ = 300e3; //except in the PhaseI case

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=false;  

  //  stackSignal_=false;
  setQuiet(true);
  stackSignal_=false;
  setStackMode(true,false,false); //stack,norm,label override

  TCut cleanup="VSPT<175 && minDeltaPhi>0.3";

  TCut color="SusyProductionMode==20000 || SusyProductionMode==10010";

  //tighter further!
  TCut noleptons = "nElectrons==0 && nMuons==0";
  TCut tauveto = "nTaus==0";
  if (useTauVeto) noleptons = noleptons&&tauveto;

  TString nbstring;
  nbstring.Form("nbjets40medium>=%d",nbtags);
  TCut btags=nbstring.Data();

  TString njstring;
  njstring.Form("njets40>=%d",njets);
  TCut tighterjets=njstring.Data();

  //  TCut baseline=cleanup&&noleptons&&tighterjets&&btags;
  TString htstring;
  htstring.Form("HT>%.0f",htCutValue);

  std::vector<TCut> cut_list;
  //  cut_list.push_back( color);
  //  cut_list.push_back( color&&noleptons);
  //  cut_list.push_back( color&&cleanup && noleptons);
  cut_list.push_back( cleanup && noleptons && tighterjets);//removed 'color' cut
  cut_list.push_back( cleanup && noleptons && tighterjets && TCut(htstring.Data()));
  cut_list.push_back( cleanup && noleptons && tighterjets && TCut(htstring.Data()) && btags);
  //cut_list.push_back( cleanup && noleptons && tighterjets && TCut(htstring.Data()) && btags);
  cut_list.push_back( cleanup && noleptons && tighterjets && TCut(htstring.Data()) && btags && TCut("MT2>400"));
  cut_list.push_back( cleanup && noleptons && tighterjets && TCut(htstring.Data()) && btags && TCut("MT2>500"));
  cut_list.push_back( cleanup && noleptons && tighterjets && TCut(htstring.Data()) && btags && TCut("MT2>600"));
  cut_list.push_back( cleanup && noleptons && tighterjets && TCut(htstring.Data()) && btags && TCut("MT2>700"));
  cut_list.push_back( cleanup && noleptons && tighterjets && TCut(htstring.Data()) && btags && TCut("MT2>800"));

  const  double dB=0.5;

  savePlots_=false;
  cout<<"Fractional error on background assumed to be "<<dB<<endl;
  for (int k=(int)samples_.size()-1;k>=0;k--) if (isSampleSM(samples_.at(k))) cout<<" & "; //one divider per SM sample
  cout<<"   & \\multicolumn{c}{S/\\sqrt{S+B+\\deltaB^2}} & \\multicolumn{c}{Zbi} \\\\"<<endl;
  for (int k=(int)samples_.size()-1;k>=0;k--) if (isSampleSM(samples_.at(k))) cout<<samples_.at(k)<<" & "; //one divider per SM sample
  cout<<"  SM & NM1 & NM2 & NM3 & NM1 & NM2 & NM3 & NM1 & NM2 & NM3 \\\\"<<endl;
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
    double S1 = getIntegral("naturalModel1");
    double S2 = getIntegral("naturalModel2");
    double S3 = getIntegral("naturalModel3");

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

void leptonEta() {


  initSamples("combinesm htskim");
  setOutputDirectory("DelphesMT2");
  lumiScale_=300e3; //discovery point
  //  treatAllAsSM_=true; //don't treat signal as signal

  int nbins;
  float low,high;
  TString var,xtitle;

  doOverflowAddition(true);
  doRatio_=true;  

  //  stackSignal_=false;

  stackSignal_=false;
  setStackMode(false,false,false); //stack,norm,label override

  TCut cleanup="VSPT<175 && minDeltaPhi>0.3";

  //tighter further!

  TCut noleptons = "nElectrons==0 && nMuons==0";
  TCut tightermt2 = "MT2>400";
  TCut tighterjets = "njets40>=7";
  TCut ht2000="HT>2000";

  clearSamples();
  addSample("tt-4p-0-600-v1510_14TEV:nElectrons+nMuons==0",kRed,"0 e+#mu");
  addSample("tt-4p-0-600-v1510_14TEV",kBlack,"No lepton veto");

  //now the full selection
  TCut ht1000="HT>1000";
   TCut loosermt2 = "MT2>=200";
  TCut btags = "nbjets40medium>=3";
  selection_ = cleanup  && ht1000 && loosermt2 && tighterjets;//&&btags;
  nbins=60; low=-6; high=6;
  var="genLepEta[0]*(abs(genLepEta[0])<abs(genLepEta[1]))+genLepEta[1]*(abs(genLepEta[0])>=abs(genLepEta[1]))"; xtitle="Gen Lepton #eta";
  drawPlots(var,nbins,low,high,xtitle,"Events", "hadronic_mt2-nob_lostLeptonEta",0,"");

  selection_ = cleanup  && ht1000 && loosermt2 && tighterjets&&btags;
  nbins=60; low=-6; high=6;
  var="genLepEta[0]*(abs(genLepEta[0])<abs(genLepEta[1]))+genLepEta[1]*(abs(genLepEta[0])>=abs(genLepEta[1]))"; xtitle="Gen Lepton #eta";
  drawPlots(var,nbins,low,high,xtitle,"Events", "hadronic_mt2_lostLeptonEta",0,"");


}
