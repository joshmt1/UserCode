{

  //======
  // to do -- split SUSY into 'SUSY signal' and 'SUSY background' ?
  //======

  //const TString hname="Hnjets";
  const TString hname="Hnjets_ge3b";
  //const TString hname="Hnjets_nocuts";
  const TString xtitle="Number of jets";

  //  const TString hname="HdeltaPhiMPTMET";
  //  const TString hname="HdeltaPhiMPTMET_ge2b";
  //  const TString xtitle="#Delta #Phi (MET,MPT)";
  
  const TString signalfile = "plots.mMSSM.root";

  //

  const TString ytitle="Weighted events";

  FileHolder fh;

  //samples!
  TFile fqcd("plots.QCD.root"); //combined with hadd
  TFile fttbar("plots.TTbarJets.root");
  TFile fwjets("plots.WJets.root");
  TFile fzjets("plots.ZJets.root");
  TFile fzinvisible("plots.Zinvisible.root");
  TFile fsusy(signalfile);

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