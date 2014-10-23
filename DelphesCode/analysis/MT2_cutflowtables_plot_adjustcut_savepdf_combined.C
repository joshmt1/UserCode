//void drawme()
{

  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();

  gROOT->LoadMacro("CMS_lumi.C");

  writeExtraText = true;   //need Preliminary
  cmsTextSize *= 1.5;
  lumiTextSize = cmsTextSize*0.9;
  cmsText     = "CMS Phase I/II Delphes Simulation";

  int iPeriod = 14;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 14= PU=140,14TeV 

  TGraph * gr[3];
  TLegend theleg(0.7,0.15,0.9,0.4);
  TCanvas thecanvas("thecanvas","",800,600);
  TString drawopt="AL";
  int colors[3]={2,4,6};//colors request by ken
  float max=0;
  //now load my graphs, add decorations, and save as pdf
  for (Long_t model = 1; model<=3; model++) {

    TString modelname = "NM";
    modelname+=model;

    //    TString fname = TString("DelphesMT2/MT2_significance_NM")+model; fname.Append(".root");
    //Update -- use Phase1+Phase2 curves
    TString fname = TString("DelphesMT2/MT2_significance_Phase12_NM")+model; fname.Append(".root");
    TFile fin(fname);
    gr[model-1] = (TGraph*) fin.Get("mt2_550_0p50");
    gr[model-1]->Draw(drawopt);
    gr[model-1]->SetLineColor(colors[model-1]);
    gr[model-1] ->GetHistogram()->GetXaxis()->SetRangeUser(0,3000);
    theleg.AddEntry(gr[model-1],modelname);
    drawopt="L";

    gr[model-1]->SetFillColor(0);

    float thismax= gr[model-1]->GetHistogram()->GetMaximum();
    if (thismax>max) max=thismax;

  }

  if (max >12) max=12;//ken's request
  gr[0]->GetHistogram()->SetMaximum(max);

  gr[0]->GetHistogram()->SetXTitle("Luminosity (fb^{-1})");
  gr[0]->GetHistogram()->SetYTitle("Expected significance (#sigma)");

  int font=42;
  gr[0]->GetHistogram()->GetYaxis()->SetLabelFont(font);
  gr[0]->GetHistogram()->GetYaxis()->SetTitleOffset(0.9);
  gr[0]->GetHistogram()->GetYaxis()->SetTitleFont(font);
  gr[0]->GetHistogram()->GetXaxis()->SetLabelFont(font);
  gr[0]->GetHistogram()->GetXaxis()->SetTitleFont(font);

  theleg.SetBorderSize(0);
  theleg.SetLineStyle(0);
  theleg.SetTextFont(42);
  theleg.SetFillStyle(0);
  theleg.Draw();
  thecanvas.GetPad(0)->SetRightMargin(0.07);
    CMS_lumi( &thecanvas, iPeriod,11 );
    TString    fname = TString("DelphesMT2/MT2_significance_0p50.pdf");

    thecanvas->SaveAs(fname);

}

