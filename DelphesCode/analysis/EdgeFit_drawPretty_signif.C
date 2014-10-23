{
  //set this string to choose what plot the script makes
  TString which ="fracyielderror" ; //or "fracyielderror", "masserror", "significance"
  
  TString summaryfile,graphname,ytitle,outputname;
  float titleoffset=1;
  float ymaxscalefactor=1;
  int headerAlignment=11;
  if (which=="significance") {
    summaryfile = "DelphesEdge/unbinned0p05fixed70_fits_summary.root";
    graphname = "significance";
    ytitle="Expected significance (#sigma)";
    outputname="significanceLee.pdf";
    titleoffset=1;
  }
  else if (which=="fracyielderror") {
    summaryfile = "DelphesEdge/unbinned0p05_fits_summary.root";
    graphname="yield_fracerror";
    ytitle="Fractional uncertainty on edge yield";
    outputname="EdgeYieldFracError.pdf";
    titleoffset=1.45;
    ymaxscalefactor=1.1;
    headerAlignment=33;
  }
  else if (which=="masserror") {
    summaryfile = "DelphesEdge/unbinned0p05_fits_summary.root";
    graphname="mass_sigma";
    ytitle="Uncertainty on edge mass (GeV)";
    outputname="MassGausSigma.pdf";
    titleoffset=1.3;
    ymaxscalefactor=1.06;
    headerAlignment=33;
  }
  else assert(0);

  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();

  gROOT->LoadMacro("CMS_lumi.C");
  lumi_14TeV="PU=140";
  writeExtraText = true;       // if extra text
  cmsTextSize *= 1.4;
  lumiTextSize = cmsTextSize*0.9;

  int iPeriod = 14;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 14= PU=140,14TeV 

  TFile f(summaryfile);
  TGraph * thegraph = (TGraph*) f.Get(graphname);

  thegraph->GetHistogram()->SetYTitle(ytitle);
  thegraph->GetHistogram()->SetXTitle("Luminosity (fb^{-1})");


  TCanvas * thecanvas=new TCanvas("thecanvas","",800,600);
  thegraph->Draw("PAL");
  thegraph->SetMarkerStyle(kCircle);
  thegraph->SetMarkerSize(1.5);
  thegraph->SetLineWidth(2);
  thegraph->GetHistogram()->SetMinimum(0); //y axis minimum
  float max=  thegraph->GetHistogram()->GetMaximum(); 
  thegraph->GetHistogram()->SetMaximum(ymaxscalefactor*max); //y axis minimum

  const float titlescale=1.3;
  float size0 = thegraph->GetHistogram()->GetYaxis()->GetTitleSize();
  thegraph->GetHistogram()->GetYaxis()->SetTitleSize( size0 * titlescale);
  size0 = thegraph->GetHistogram()->GetXaxis()->GetTitleSize();
  thegraph->GetHistogram()->GetXaxis()->SetTitleSize( size0 * titlescale);

  const float labelscale=1.15;
  size0 = thegraph->GetHistogram()->GetXaxis()->GetLabelSize();
  thegraph->GetHistogram()->GetXaxis()->SetLabelSize( size0 * labelscale);
  size0 = thegraph->GetHistogram()->GetYaxis()->GetLabelSize();
  thegraph->GetHistogram()->GetYaxis()->SetLabelSize( size0 * labelscale);

  int font=42;
  thegraph->GetHistogram()->GetXaxis()->SetTitleFont(font);
  thegraph->GetHistogram()->GetYaxis()->SetTitleFont(font);
  thegraph->GetHistogram()->GetXaxis()->SetLabelFont(font);
  thegraph->GetHistogram()->GetYaxis()->SetLabelFont(font);
  thegraph->GetHistogram()->GetYaxis()->SetTitleOffset(titleoffset);
  CMS_lumi( thecanvas, iPeriod,headerAlignment );

  summaryfile.ReplaceAll("summary.root",outputname);
  thecanvas->SaveAs(summaryfile);


}
