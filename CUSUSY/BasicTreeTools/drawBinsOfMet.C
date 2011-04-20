//void plotSlices( const TCut btag = "nbjets>=1")
{

  const TCut btag = "nbjets>=2";

  //  const TString var ="minDeltaPhi";
  const TString var ="bestTopMass";

  const bool drawLSB = true;
  const bool drawMSB = false;
  const bool drawSB = true;
  const bool drawSIG = false;

  const float customMax = 0.4;

  //updated to use reducedTrees

  bool addOverflow = false;
  if (var=="minDeltaPhi") addOverflow=false;

  //  gROOT->SetStyle("CMS");
  gStyle->SetOptStat(0);

  TCut baseline ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passInconsistentMuon==1 && passBadPFMuon==1  && weight<1000"; //no mindp, no MET
  TCut minDPcut = "cutDeltaPhi==1";
  if (var!="minDeltaPhi") baseline = baseline&&minDPcut;
  TCut LSB = "MET<50";
  TCut MSB = "MET>=50 && MET<100";
  TCut SB = "MET>=100 && MET <150";
  TCut SIG = "MET>150";

  TCut cut1=baseline && LSB && btag;
  TCut cut2=baseline && MSB && btag;
  TCut cut3=baseline && SB && btag;
  TCut cut4=baseline && SIG && btag;

  TString selection1 = TString("weight*(")+cut1.GetTitle()+")";
  TString selection2 = TString("weight*(")+cut2.GetTitle()+")";
  TString selection3 = TString("weight*(")+cut3.GetTitle()+")";
  TString selection4 = TString("weight*(")+cut4.GetTitle()+")";

  TChain madgraph("reducedTree");
  madgraph.Add("/cu2/joshmt/V00-03-01_2/reducedTree.Baseline0_PF_JERbias_pfMEThigh_PFLepRA20e0mu_minDP_MuonCleaning.QCD.root");

  TChain pypu("reducedTree");
  pypu.Add("/cu2/joshmt/V00-03-01_2/reducedTree.Baseline0_PF_JERbias_pfMEThigh_PFLepRA20e0mu_minDP_MuonCleaning.PythiaPUQCD.root");

  TChain py("reducedTree");
  py.Add("/cu2/joshmt/V00-03-01_2/reducedTree.Baseline0_PF_JERbias_pfMEThigh_PFLepRA20e0mu_minDP_MuonCleaning.PythiaQCD.root");

  TChain pypuflat("reducedTree");
  pypuflat.Add("/cu2/joshmt/V00-03-01_2/reducedTree.Baseline0_PF_pfMEThigh_PFLepRA20e0mu_minDP_MuonCleaning.PythiaPUQCDFlat.root");

  TH1::SetDefaultSumw2(); //trick to turn on Sumw2 for all histos

  //const int nbins=6;
  //  const double varbins[]={0.,160.,180.,260.,400.,800.,2000};
  int nbins=40;
  float min=0;
  float max=TMath::Pi();
  if (var=="bestTopMass") max=800;
  
  TH1D Hlow("Hlow","",nbins,min,max);
  TH1D Hmed("Hmed","",nbins,min,max);
  TH1D Hmh("Hmh","",nbins,min,max);
  TH1D Hhigh("Hhigh","",nbins,min,max);
  /*
  TH1D Hlow("Hlow","",nbins,varbins);
  TH1D Hmed("Hmed","",nbins,varbins);
  TH1D Hmh("Hmh","",nbins,varbins);
  TH1D Hhigh("Hhigh","",nbins,varbins);
  */
  Hlow.SetLineColor(kBlue);
  Hmed.SetLineColor(1);
  Hmh.SetLineColor(kRed);
  Hhigh.SetLineColor(6);

  float width=2.5;
  Hlow.SetLineWidth(width);
  Hmed.SetLineWidth(width);
  Hmh.SetLineWidth(width);
  Hhigh.SetLineWidth(width);

  //was minDeltaPhi, bestTopMass
  if (drawLSB) pypu.Draw(var+">>Hlow",selection1);
  if (drawMSB) pypu.Draw(var+">>Hmed",selection2);
  if (drawSB) pypu.Draw(var+">>Hmh",selection3);
  if (drawSIG) pypu.Draw(var+">>Hhigh",selection4);

  if (addOverflow) {
    Hlow.SetBinContent(nbins, Hlow.GetBinContent(nbins)+Hlow.GetBinContent(nbins+1));
    Hmed.SetBinContent(nbins, Hmed.GetBinContent(nbins)+Hmed.GetBinContent(nbins+1));
    Hmh.SetBinContent(nbins, Hmh.GetBinContent(nbins)+Hmh.GetBinContent(nbins+1));
    Hhigh.SetBinContent(nbins, Hhigh.GetBinContent(nbins)+Hhigh.GetBinContent(nbins+1));
  }

  if (Hlow.Integral()>0) Hlow.Scale( 1.0 / Hlow.Integral());
  if (Hmed.Integral()>0) Hmed.Scale( 1.0 / Hmed.Integral());
  if (Hmh.Integral()>0) Hmh.Scale( 1.0 / Hmh.Integral());
  if (Hhigh.Integral()>0) Hhigh.Scale( 1.0 / Hhigh.Integral());
  
  if (var=="minDeltaPhi")  Hhigh.SetXTitle("#Delta #phi_{min}");
  else if (var=="bestTopMass")  {
    TString title="best 3-jet mass (GeV)";
    Hhigh.SetXTitle(title);
    Hmh.SetXTitle(title);
    Hmed.SetXTitle(title);
    Hlow.SetXTitle(title);
  }
  TString ytitle="Arbitrary units";
  Hhigh.SetYTitle(ytitle);
  Hmh.SetYTitle(ytitle);
  Hmed.SetYTitle(ytitle);
  Hlow.SetYTitle(ytitle);
  TString thetitle=btag.GetTitle();
  Hhigh.SetTitle(thetitle);
  Hmh.SetTitle(thetitle);
  Hmed.SetTitle(thetitle);
  Hlow.SetTitle(thetitle);


  TString drawopt="hist e";
  if (drawSIG)  {Hhigh.Draw(drawopt); drawopt="hist e SAME"; if (customMax>0) Hhigh.SetMaximum(customMax);}
  if (drawSB)   {Hmh.Draw(drawopt); drawopt="hist e SAME";   if (customMax>0) Hmh.SetMaximum(customMax);}
  if (drawMSB)  {Hmed.Draw(drawopt); drawopt="hist e SAME";  if (customMax>0) Hmed.SetMaximum(customMax);}
  if (drawLSB)  {Hlow.Draw(drawopt); drawopt="hist e SAME";  if (customMax>0) Hlow.SetMaximum(customMax);}

  if (var=="minDeltaPhi") {
    Hhigh.SetMinimum(0);
    Hhigh.GetXaxis()->SetRangeUser(0,2);
  }

  TLegend leg(0.55,0.55,0.85,0.85);
  leg.SetFillColor(0);
  if (drawLSB)  leg.AddEntry(&Hlow,"E_{T}^{miss} < 50 GeV (LSB)");
  if (drawMSB)  leg.AddEntry(&Hmed,"50 < E_{T}^{miss} < 100 GeV");
  if (drawSB)   leg.AddEntry(&Hmh,"100 < E_{T}^{miss} < 150 GeV (SB)");
  if (drawSIG)  leg.AddEntry(&Hhigh,"E_{T}^{miss} > 150 GeV");
  leg.Draw();

//   double chi2=0;
//   for (int i=1; i<=nbins; i++) {
//     double denom=Hlow.GetBinError(i)*Hlow.GetBinError(i) + Hmh.GetBinError(i)*Hmh.GetBinError(i);
//     double c2= denom>0 ? pow( Hlow.GetBinContent(i) - Hmh.GetBinContent(i) ,2) / denom : 0;
//     chi2+=c2;
//   }

//   cout<<"Hand chi^2 = "<<chi2<<endl;
//   Hlow.Chi2Test(&Hmh,"WW p");

}

// plotAll() {

//   TCanvas cplots("cplots","plots",900,600);
//   cplots.Divide(1,3);

//   cplots.cd(1);
//   plotSlices("nbjets==1");

//   cplots.cd(2);
//   plotSlices("nbjets>=1");

//   cplots.cd(3);
//   plotSlices("nbjets>=2");

//   cplots.SaveAs("compM3jsbvslsb.eps");

// }
