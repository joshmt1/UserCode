//void plotMETDeltaPhiCorrelation()
{
  //updated to use reducedTrees
  //  gROOT->SetStyle("CMS");
  gStyle->SetOptStat(0);

  TCut baseline ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passInconsistentMuon==1 && passBadPFMuon==1  && weight<1000"; //&&cutDeltaPhi==1
  TCut LSB = "MET<50";
  TCut MSB = "MET>=50 && MET<100";
  TCut SB = "MET>=100 && MET <150";
  TCut SIG = "MET>150";
  TCut btag = "nbjets==1";

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

  TH1::SetDefaultSumw2(); //trick to turn on Sumw2 for all histos

  //const int nbins=6;
  //  const double varbins[]={0.,160.,180.,260.,400.,800.,2000};
  int nbins=50;
  float min=0;
  float max=TMath::Pi();
    //float max=800;
  
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
  Hlow.SetLineColor(1);
  Hmed.SetLineColor(2);
  Hmh.SetLineColor(4);
  Hhigh.SetLineColor(6);

  float width=2.5;
  Hlow.SetLineWidth(width);
  Hmed.SetLineWidth(width);
  Hmh.SetLineWidth(width);
  Hhigh.SetLineWidth(width);

  //was minDeltaPhi, bestTopMass
  pypu.Draw("minDeltaPhi>>Hlow",selection1);
  pypu.Draw("minDeltaPhi>>Hmed",selection2);
  pypu.Draw("minDeltaPhi>>Hmh",selection3);
  pypu.Draw("minDeltaPhi>>Hhigh",selection4);

  Hlow.SetBinContent(nbins, Hlow.GetBinContent(nbins)+Hlow.GetBinContent(nbins+1));
  Hmed.SetBinContent(nbins, Hmed.GetBinContent(nbins)+Hmed.GetBinContent(nbins+1));
  Hmh.SetBinContent(nbins, Hmh.GetBinContent(nbins)+Hmh.GetBinContent(nbins+1));
  Hhigh.SetBinContent(nbins, Hhigh.GetBinContent(nbins)+Hhigh.GetBinContent(nbins+1));

  if (Hlow.Integral()>0) Hlow.Scale( 1.0 / Hlow.Integral());
  if (Hmed.Integral()>0) Hmed.Scale( 1.0 / Hmed.Integral());
  if (Hmh.Integral()>0) Hmh.Scale( 1.0 / Hmh.Integral());
  if (Hhigh.Integral()>0) Hhigh.Scale( 1.0 / Hhigh.Integral());
  
  Hhigh.SetXTitle("#Delta #phi_{min}");
  Hhigh.SetYTitle("Arbitrary units");

  Hhigh.Draw();
  Hmh.Draw("same");
  Hmed.Draw("same");
  Hlow.Draw("same");

  Hhigh.SetMinimum(0);
  Hhigh.GetXaxis()->SetRangeUser(0,2);

  TLegend leg(0.4,0.4,0.85,0.8);
  leg.SetFillColor(0);
  leg.AddEntry(&Hlow,"E_{T}^{miss} < 50 GeV");
  leg.AddEntry(&Hmed,"50 < E_{T}^{miss} < 100 GeV");
  leg.AddEntry(&Hmh,"100 < E_{T}^{miss} < 150 GeV");
  leg.AddEntry(&Hhigh,"E_{T}^{miss} > 150 GeV");
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
