//void drawme()
{

  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();

  gROOT->LoadMacro("CMS_lumi.C");

  writeExtraText = false;       // if extra text
  cmsTextSize = 0.53;
  lumiTextSize = cmsTextSize*0.9;

  int iPeriod = 14;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 14= PU=140,14TeV 

  //now load my TCanvases, add decorations, and save as pdf
  for (Long_t model = 1; model<=3; model++) {

    TString modelname = "NM";
    modelname+=model;

    TString fname = TString("DelphesMT2/MT2_significance_NM")+model; fname.Append(".root");
    TFile fin(fname);
    TCanvas* thecanvas = (TCanvas*) fin.Get("thecanvas");
    TCanvas* thecanvas2 = (TCanvas*) fin.Get("thecanvas2");
    thecanvas->Draw();
    CMS_lumi( thecanvas, iPeriod,11 );
    TText modelbox(0.15,0.75,modelname.Data());
    modelbox.SetNDC();
    modelbox.Draw();

    fname.ReplaceAll(".root",".pdf");
    thecanvas->SaveAs(fname);
    thecanvas2->Draw();
    fname = TString("DelphesMT2/MT2_significance_bestcut_NM")+model; fname.Append(".pdf");
    thecanvas2->SaveAs(fname);
  }

}

