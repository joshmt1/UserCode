{

  //======
  // to do -- split SUSY into 'SUSY signal' and 'SUSY background' ?
  //======
  //  const TString dir="plots_V00-00-02/";
  const TString dir="";

  //const TString hname="Hnjets";
  //const TString hname="Hnjets_ge2b";
  //  const TString hname="HdeltaPhiMPTMET_ge2b";
  const TString hname="H_MET_ge2b";
  //const TString hname="Hjetpt1_ge2b";
  //const TString hname="HdeltaPhiMPTMET";
   //const TString hname="HminDeltaR_bj_ge2b";
  //const TString hname="Hnjets_nocuts";
  //  const TString xtitle="Number of jets";
  //  const TString xtitle="#Delta #phi (MPT,MET)";
    const TString xtitle="caloMET";
  //const TString xtitle="pT of jet 1";
  //const TString xtitle="event minimum DeltaR(b,non-b)";

  //  const TString hname="HdeltaPhiMPTMET";
  //  const TString hname="HdeltaPhiMPTMET_ge2b";
  //  const TString xtitle="#Delta #Phi (MET,MPT)";
  
  //    const TString cutstring = "RA2METminDP_NoMET_NoDeltaPhi";
  //  const TString cutstring = "RA2wideMETminDP";
  const TString cutstring = "RA2_calo_METwide_minDP";
  //const TString cutstring = "RA2medMET";
  const TString sigsample = "LM9";

  //  const TString customTitle =""; //set to this to not use a custom title
  const TString customTitle="10.8 pb^{-1} (CU group cuts)";

  const bool dodata=true;
  const double MCreweight = (2771.0+167.0+3898+511+3467)/(100e3); //use to adjust the MC luminosity

  const int dataColor=7;
  //range for rescaling
  const bool dorescaling=false;
  const  double lowlimit=130;
  const  double highlimit=170;

  //range for plotting
  const bool customrange=true; //do we want to customize the range?
  const double lowplotlimit=50;
  const double highplotlimit=200;

  //rebin?
  const bool dorebin=false;
  int nrebin=2;

  gStyle->SetOptStat(0); //comment out as desired


  const TString ytitle="Events";

  FileHolder fh;

  TString signalfile = "plots.";
  signalfile+=cutstring;
  signalfile+=".";
  signalfile+=sigsample;
  signalfile+=".root";

  TString filebase=dir;
  filebase+="plots.";
  filebase+=cutstring;

  
  //samples!
  //first get data
  TFile* fdata=0;
  TH1D* hdata=0;
  double norm_data=0;
  if (dodata) {
    fdata=new TFile(filebase+".data.root");
    hdata=(TH1D*) fdata->Get(hname);
    hdata->SetLineColor(dataColor);
    if (dorebin) hdata->Rebin(nrebin);
    char ccc[100];
    sprintf(ccc,"%15s %f","Data",hdata->Integral());
    cout<<ccc<<endl;
    norm_data = hdata->Integral( hdata->FindBin(lowlimit), hdata->FindBin(highlimit));
  }

  //then MC
  TFile fqcd(filebase+".QCD.root"); //combined with hadd
  TFile fttbar(filebase+".TTbarJets.root");
  TFile fwjets(filebase+".WJets.root");
  TFile fzjets(filebase+".ZJets.root");
  TFile fzinvisible(filebase+".Zinvisible.root");
  TFile fsusy(dir+signalfile);

  fh.add(&fqcd);
  fh.add(&fttbar);
  fh.add(&fwjets);
  fh.add(&fzjets);
  fh.add(&fzinvisible);
  fh.add(&fsusy);

  HistHolder hh_njets;
  THStack mystack("mystack","--");
  TLegend leg(0.6,0.6,0.9,0.9);
  leg.SetFillColor(0);

  double norm_qcd=0;
  double norm_bkg=0;

  { //should get rid of this extra brace
    for (int i=0; i<fh.size(); i++) {
      hh_njets.load(hname,fh.at(i));
      hh_njets.find(hname,fh.at(i))->SetFillColor(i+1);
      hh_njets.find(hname,fh.at(i))->SetLineColor(i+1);

      hh_njets.find(hname,fh.at(i))->Scale(MCreweight);
      if (dorebin)  hh_njets.find(hname,fh.at(i))->Rebin(nrebin);

 //      hh_njets.find(hname,fh.at(i))->SetXTitle("Number of Jets");
//       hh_njets.find(hname,fh.at(i))->SetYTitle("Weighted events");
      //      hh_njets.find(hname,fh.at(i))->SetTitle("100 pb^{-1}");
      //hh_njets.find(hname,fh.at(i))->SetLineWidth(3);

      mystack.Add( hh_njets.find(hname,fh.at(i)) );

      //copy title from histograms to histogram stack
      if (i==0) {
	if (customTitle!="") mystack.SetTitle(customTitle);
	else	mystack.SetTitle(hh_njets.find(hname,fh.at(i))->GetTitle());
      }

      TString filename=fh.at(i)->GetName();
      TObjArray* pieces = filename.Tokenize(".");
      //gonna have to be careful to stick with filename convention = plots.blah.sampleName.root
      TString sampleName=pieces->At(2)->GetName();
      leg.AddEntry( hh_njets.find(hname,fh.at(i)), sampleName);

      char ccc[100];
      sprintf(ccc,"%15s %f",sampleName.Data(),hh_njets.find(hname,fh.at(i))->Integral());
      cout<<ccc<<endl;
      if (sampleName == "QCD") {
	norm_qcd=hh_njets.find(hname,fh.at(i))->Integral(hh_njets.find(hname,fh.at(i))->FindBin(lowlimit),hh_njets.find(hname,fh.at(i))->FindBin(highlimit));
      }
      else if (!sampleName.Contains("LM") ) {
	norm_bkg+= hh_njets.find(hname,fh.at(i))->Integral(hh_njets.find(hname,fh.at(i))->FindBin(lowlimit),hh_njets.find(hname,fh.at(i))->FindBin(highlimit));}
    }
  }

  double rescale=1;
  if (dorescaling) {
    cout<<"number of data events in low MET region after subtracting non-QCD background = "<<norm_data - norm_bkg<<endl;
    rescale = (norm_data - norm_bkg)/norm_qcd;
    cout<<"rescale = "<<rescale<<endl;
  }

  for (int i=0; i<fh.size(); i++) {
    TString filename=fh.at(i)->GetName();
    TObjArray* pieces = filename.Tokenize(".");
    //gonna have to be careful to stick with filename convention = plots.blah.sampleName.root
    TString sampleName=pieces->At(2)->GetName();
    if (dorescaling && sampleName=="QCD")    hh_njets.find(hname,fh.at(i))->Scale(rescale);
  }

  //  hh_njets.find(hname,fh.at(0))->Draw("HIST");

  mystack.Draw("HIST");
  //  mystack.Draw("SAME");

  mystack.GetHistogram()->GetXaxis()->SetTitle(xtitle);
  mystack.GetHistogram()->GetYaxis()->SetTitle(ytitle);
  if (customrange)  mystack.GetHistogram()->GetXaxis().SetRangeUser(lowplotlimit,highplotlimit);

  hdata->SetMarkerColor(dataColor);
  hdata->SetLineWidth(2);
  hdata->Draw("SAME");

  leg.Draw();

}
