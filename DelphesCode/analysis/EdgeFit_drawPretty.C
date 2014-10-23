//void drawme()
{
  //for unbinned
  Long_t seed=106986;
  Long_t lumi=3000;

  //for 3000fb-1 binned fit
  //seed=987;
  //  lumi=3000;

  //for 500fb-1 binned fit
  //seed=102000;
  //lumi=500;

  int plotToDraw=1; //1 is SF, 2 is OF, 3 is binned
  TString plotName,histName,saveStub,totalPdfName;
  TString path="DelphesEdge/"; 
  TString file;

  if (plotToDraw==1) {
    plotName="fullfitframe";
    histName="h_data_Cut[sample==sample::SF]";
    saveStub="_SF.pdf";
    totalPdfName="pdf_total_Norm[mll]";
    path+=lumi; path+="/";
    file="unbinnedWithConstraint0p05_3000_";
    file+=seed;
  }
  else if (plotToDraw==2) {
    plotName="fullfitframe_of";
    histName="h_data_Cut[sample==sample::OF]";
    saveStub="_OF.pdf";
    totalPdfName="pdf_total_of_Norm[mll]";
    path+=lumi; path+="/";
    file="unbinnedWithConstraint0p05_3000_";
    file+=seed;
  }
  else if (plotToDraw==3) {
    if (lumi==3000) {
      plotName=";";
      plotName+=seed;
    }
    else if (lumi==500) {
      plotName="binnedfitplot_";
      plotName+=lumi;
      plotName+="_";
      plotName+=seed;
    }
    else assert(0);
    file="binnedfitplots_";
    file+=lumi;
    histName="h_data";
    totalPdfName="pdf_total_Norm[mll]";
    saveStub.Form("_%d.pdf",(int)seed);
  }
  else assert(0);
  file+=".root";

  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();

  gROOT->LoadMacro("CMS_lumi.C");
  lumi_14TeV.Form("%d",(int)lumi);
  lumi_14TeV+=" fb^{-1}, PU=140";

  writeExtraText = true;       // if extra text
  //  cmsTextSize *= 1.1;
  //  lumiTextSize = cmsTextSize*0.9;

  int iPeriod = 14;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 14= PU=140,14TeV 

  //load the file

  TString filepath=path+file;

  TFile f(filepath);

  TCanvas *c=new TCanvas("c","",800,600);
  RooPlot* fullfit = (RooPlot*) f.Get(plotName);
  fullfit->Draw();

  int font=42;
  fullfit->GetXaxis()->SetTitleFont(font);
  fullfit->GetYaxis()->SetTitleFont(font);

  float sz=0.05;
  fullfit->GetXaxis()->SetTitleSize(sz);
  fullfit->GetYaxis()->SetTitleSize(sz);

  fullfit->GetXaxis()->SetTitle("m_{l+l-} (GeV)");
  TString thetitle=  fullfit->GetYaxis()->GetTitle();
  thetitle.ReplaceAll("(","");
  thetitle.ReplaceAll(")","");
  thetitle+="GeV";
  fullfit->GetYaxis()->SetTitle(thetitle.Data());

  fullfit->GetXaxis()->SetLabelFont(font);
  fullfit->GetYaxis()->SetLabelFont(font);

  RooHist* dh=fullfit->getHist(histName);
  RooCurve* pdfback= plotToDraw==1 ? fullfit->getCurve("pdf_total_Norm[mll]_Comp[pdf_fs,pdf_dy]") : 0;

  RooCurve* pdftotal=fullfit->getCurve(totalPdfName);
  RooCurve* pdffs= plotToDraw==1 ? fullfit->getCurve("pdf_total_Norm[mll]_Comp[pdf_fs]") : 0;

  dh->SetFillColor(0);
  if (plotToDraw==1) {
    pdfback->SetFillColor(0);
    pdffs->SetFillColor(0);
  }
  pdftotal->SetFillColor(0);

  TLegend theleg(0.7,0.55,0.94,0.8);
  theleg.AddEntry(dh,"Pseudodata","eplf");
  if (plotToDraw==1) {
    theleg.AddEntry(pdftotal,"Total PDF","L");
    theleg.AddEntry(pdfback,"FS+Z PDF","L");
    theleg.AddEntry(pdffs,"FS PDF","L");
  }
  else if (plotToDraw==2) {
    theleg.AddEntry(pdftotal,"FS PDF","L");
  }
  else if (plotToDraw==3) {
    theleg.AddEntry(pdftotal,"Total PDF","L");
  }
  else assert(0);
  theleg.SetBorderSize(0);
  theleg.SetLineStyle(0);
  theleg.SetTextFont(42);
  theleg.SetFillStyle(0);
  theleg.Draw();

  c->cd()->SetTopMargin(0.08);


  CMS_lumi( c, iPeriod,11 );
  fullfit->GetYaxis()->SetTitleOffset(1.25);
  fullfit->GetXaxis()->SetTitleOffset(1.15);

  float minval=0;
  if (plotToDraw==3) minval= fullfit->GetMinimum();
  fullfit->GetYaxis()->SetRangeUser(minval, fullfit->GetMaximum()*1.1);

  TString outfile = filepath;
  outfile.ReplaceAll(".root",saveStub);
  c->SaveAs(outfile);

}
