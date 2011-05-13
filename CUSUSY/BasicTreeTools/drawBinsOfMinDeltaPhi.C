//void plotSlices( const TCut btag = "nbjets>=1")
{

  const TCut btag = "nbjets>=1";
  //const TCut btag = "nbjetsSSV0>=1";

  //  const TString var ="minDeltaPhi";
  //  const TString var ="bestTopMass";
  const TString var ="MET";
  //const TString var ="njets";

  const bool drawSB=false;

  int nbins=10;
  float min=0;
  float max=TMath::Pi();
  if (var=="bestTopMass") max=800;
  else if (var=="MET") {min=100; max=200;}
  else if (var=="njets") {min=2; max=10; nbins=8;}

  const float customMax = 1;
  const bool doRatio=true; //hack for a specific ratio plot

  bool doNJetReweighting=false;

  //updated to use reducedTrees

  bool addOverflow = true;

  TCut baseline ="cutHT==1 && cutPV==1 && cutTrigger==1 && cut3Jets==1 && cutEleVeto==1 && cutMuVeto==1 && passInconsistentMuon==1 && passBadPFMuon==1  && weight<1000 && njets==3"; //no mindp, no MET
  //  TCut LSB = "MET<50";
  TCut MSB = "MET>=50 && MET<100";
  TCut SB = "MET>=100 && MET <150";
  TCut SIG = "MET>100";

  TCut passMinDp = "cutDeltaPhi==1";
  TCut failMinDp = "cutDeltaPhi==0";

  TCut cut1=baseline && SB && btag && passMinDp;
  TCut cut2=baseline && SB && btag && failMinDp;
  TCut cut3=baseline && SIG && btag && passMinDp;
  TCut cut4=baseline && SIG && btag && failMinDp;

  TString selection1 = TString("weight*(")+cut1.GetTitle()+")";
  TString selection2 = TString("weight*(")+cut2.GetTitle()+")";
  TString selection3 = TString("weight*(")+cut3.GetTitle()+")";
  TString selection4 = TString("weight*(")+cut4.GetTitle()+")";

  if (doNJetReweighting) {
    TString reweight = "*((njets==3)*0.66 + (njets==4)*1.29 + (njets==5)*3.81 + (njets==6)*3.44 + (njets==7)*7.01 + (njets==8)*5.9 + (njets>=9)*10.13)";
    //    selection1 += reweight;
    //    selection2 += reweight;
    //    selection3 += reweight;
    selection4 += reweight;
  }

  TChain madgraph("reducedTree");
  madgraph.Add("/cu2/joshmt/V00-03-01_6/reducedTree.Baseline0_PF_JERbias6_pfMEThigh_PFLepRA20e0mu_minDP_MuonCleaning.QCD.root");

  TChain pypu("reducedTree");
  pypu.Add("/cu2/joshmt/V00-03-01_6/reducedTree.Baseline0_PF_JERbias6_pfMEThigh_PFLepRA20e0mu_minDP_MuonCleaning.PythiaPUQCD.root");

  TChain py("reducedTree");
  py.Add("/cu2/joshmt/V00-03-01_6/reducedTree.Baseline0_PF_JERbias6_pfMEThigh_PFLepRA20e0mu_minDP_MuonCleaning.PythiaQCD.root");

  TH1::SetDefaultSumw2(); //trick to turn on Sumw2 for all histos

  //const int nbins=6;
  //  const double varbins[]={0.,160.,180.,260.,400.,800.,2000};

  int height= doRatio ? 800 : 600;
  TCanvas * thecanvas= new TCanvas("thecanvas","the canvas",700,height);
  if (doRatio) {
    thecanvas->Divide(1,2);
    const float padding=0.01; const float ydivide=0.2;
    thecanvas->GetPad(1)->SetPad( padding, ydivide + padding, 1-padding, 1-padding);
    thecanvas->GetPad(2)->SetPad( padding, padding, 1-padding, ydivide-padding);
    thecanvas->cd(1);
  }
  
  TH1D Hfail("Hfail","",nbins,min,max);
  TH1D Hpass("Hpass","",nbins,min,max);
  TH1D Hratio("Hratio","",nbins,min,max);

  TH1D Hfail1("Hfail1","",nbins,min,max);
  TH1D Hpass1("Hpass1","",nbins,min,max);
  /*
  TH1D Hlow("Hlow","",nbins,varbins);
  TH1D Hmed("Hmed","",nbins,varbins);
  TH1D Hmh("Hmh","",nbins,varbins);
  TH1D Hhigh("Hhigh","",nbins,varbins);
  */
  Hfail.SetLineColor(kBlue);
  Hpass.SetLineColor(kMagenta);

  Hfail1.SetLineColor(kBlue-9);
  Hpass1.SetLineColor(kMagenta-9);

  float width=2.5;
  Hfail.SetLineWidth(width);
  Hpass.SetLineWidth(width);

  Hfail1.SetLineWidth(width);
  Hpass1.SetLineWidth(width);

  pypu.Draw(var+">>Hfail",selection4);
  pypu.Draw(var+">>Hpass",selection3);

  if (drawSB) {
    pypu.Draw(var+">>Hfail1",selection2);
    pypu.Draw(var+">>Hpass1",selection1);
  }

  if (Hfail.Integral()>0) Hfail.Scale( 1.0 / Hfail.Integral());
  if (Hpass.Integral()>0) Hpass.Scale( 1.0 / Hpass.Integral());

  if (drawSB) {
    if (Hfail1.Integral()>0) Hfail1.Scale( 1.0 / Hfail1.Integral());
    if (Hpass1.Integral()>0) Hpass1.Scale( 1.0 / Hpass1.Integral());
  }

  if (var=="minDeltaPhi")  Hpass.SetXTitle("#Delta #phi_{min}");
  else if (var=="bestTopMass")  {
    TString title="best 3-jet mass (GeV)";
    Hpass.SetXTitle(title);
    Hfail.SetXTitle(title);
  }
  else {
    Hpass.SetXTitle(var);
    Hfail.SetXTitle(var);
    Hpass1.SetXTitle(var);
    Hfail1.SetXTitle(var);
  }

  TString ytitle="Arbitrary units";
  Hpass.SetYTitle(ytitle);
  Hfail.SetYTitle(ytitle);
  Hpass1.SetYTitle(ytitle);
  Hfail1.SetYTitle(ytitle);
  TString thetitle=btag.GetTitle();
  Hpass.SetTitle(thetitle);
  Hfail.SetTitle(thetitle);
  Hpass1.SetTitle(thetitle);
  Hfail1.SetTitle(thetitle);

  TString drawopt="hist e";
  Hpass.Draw(drawopt); drawopt="hist e SAME"; if (customMax>0) Hpass.SetMaximum(customMax);
  Hfail.Draw(drawopt); drawopt="hist e SAME";  if (customMax>0) Hfail.SetMaximum(customMax);
  if (drawSB) {
    Hpass1.Draw(drawopt); drawopt="hist e SAME"; if (customMax>0) Hpass1.SetMaximum(customMax);
    Hfail1.Draw(drawopt); drawopt="hist e SAME";  if (customMax>0) Hfail1.SetMaximum(customMax);
  }

  if (var=="minDeltaPhi") {
    Hpass.SetMinimum(0);
    Hpass.GetXaxis()->SetRangeUser(0,2);
  }

  TLegend leg(0.55,0.55,0.85,0.85);
  leg.SetFillColor(0);
  leg.AddEntry(&Hfail,"fail minDeltaPhi cut");
  leg.AddEntry(&Hpass,"pass minDeltaPhi cut");
  if (drawSB) {
    leg.AddEntry(&Hfail1,"[low MET] fail minDeltaPhi cut");
    leg.AddEntry(&Hpass1,"[low MET] pass minDeltaPhi cut");
  }
  leg.Draw();

  if (doRatio) {
    Hratio.Divide(&Hpass,&Hfail);
    Hratio.SetLineWidth(2);
    thecanvas->cd(2);
    Hratio.Draw();
    Hratio.SetMinimum(0);
    Hratio.SetMaximum(3);
    TLine* l1 = new TLine(min,1,max,1);
    l1->SetLineColor(kMagenta);
    l1->SetLineWidth(2);
    l1->Draw();
    cout<<"KS test results = "<<Hpass.KolmogorovTest(&Hfail)<<endl;
  }

//   double chi2=0;
//   for (int i=1; i<=nbins; i++) {
//     double denom=Hfail.GetBinError(i)*Hfail.GetBinError(i) + Hmh.GetBinError(i)*Hmh.GetBinError(i);
//     double c2= denom>0 ? pow( Hfail.GetBinContent(i) - Hmh.GetBinContent(i) ,2) / denom : 0;
//     chi2+=c2;
//   }

//   cout<<"Hand chi^2 = "<<chi2<<endl;
//   Hfail.Chi2Test(&Hmh,"WW p");

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
