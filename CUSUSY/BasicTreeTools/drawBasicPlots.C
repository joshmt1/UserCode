{

  //======
  // to do -- split SUSY into 'SUSY signal' and 'SUSY background' ?
  //======
  //  const TString dir="plots_V00-00-02/";
  const TString dir="";

  //const TString hname="Hnjets";
  //  const TString hname="Hnjets_ge3b";
  //  const TString hname="HdeltaPhiMPTMET_ge2b";
  const TString hname="HminDeltaR_bj_ge2b";
  //const TString hname="Hnjets_nocuts";
  //  const TString xtitle="Number of jets";
  //  const TString xtitle="DeltaPhi(MPT,MET)";
  const TString xtitle="event minimum DeltaR(b,non-b)";

  //  const TString hname="HdeltaPhiMPTMET";
  //  const TString hname="HdeltaPhiMPTMET_ge2b";
  //  const TString xtitle="#Delta #Phi (MET,MPT)";
  
  const TString signalfile = "plots.RA2MET.LM9.root";

  //

  const TString ytitle="Weighted events";

  FileHolder fh;

  //samples!
  TFile fqcd(dir+"plots.RA2MET.QCD.root"); //combined with hadd
  TFile fttbar(dir+"plots.RA2MET.TTbarJets.root");
  TFile fwjets(dir+"plots.RA2MET.WJets.root");
  TFile fzjets(dir+"plots.RA2MET.ZJets.root");
  TFile fzinvisible(dir+"plots.RA2MET.Zinvisible.root");
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

  {
    for (int i=0; i<fh.size(); i++) {
      hh_njets.load(hname,fh.at(i));
      hh_njets.find(hname,fh.at(i))->SetFillColor(i+1);
      hh_njets.find(hname,fh.at(i))->SetLineColor(i+1);

 //      hh_njets.find(hname,fh.at(i))->SetXTitle("Number of Jets");
//       hh_njets.find(hname,fh.at(i))->SetYTitle("Weighted events");
      //      hh_njets.find(hname,fh.at(i))->SetTitle("100 pb^{-1}");
      //hh_njets.find(hname,fh.at(i))->SetLineWidth(3);

      mystack.Add( hh_njets.find(hname,fh.at(i)) );

      //copy title from histograms to histogram stack
      if (i==0) mystack.SetTitle(hh_njets.find(hname,fh.at(i))->GetTitle());

      TString filename=fh.at(i)->GetName();
      TObjArray* pieces = filename.Tokenize(".");
      //gonna have to be careful to stick with filename convention = plots.sampleName.root
      TString sampleName=pieces->At(1)->GetName();
      leg.AddEntry( hh_njets.find(hname,fh.at(i)), sampleName);

      char ccc[100];
      sprintf(ccc,"%15s %f",sampleName.Data(),hh_njets.find(hname,fh.at(i))->Integral());
      cout<<ccc<<endl;
    }
  }

  //  hh_njets.find(hname,fh.at(0))->Draw("HIST");

  mystack.Draw("HIST");
  //  mystack.Draw("SAME");

  mystack.GetHistogram()->GetXaxis()->SetTitle(xtitle);
  mystack.GetHistogram()->GetYaxis()->SetTitle(ytitle);

  leg.Draw();

}
