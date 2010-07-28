/*

.L basicLoop.C++

*/

{

  gSystem->Load("basicLoop_C.so");

  {
  //ttbar
  TChain* ttbarChain = new TChain("BasicTreeMaker/tree");
  ttbarChain->Add("~/data/BasicNtuples/V00-00-01/TTbarJets/Spring10/BasicNtuple*");

  basicLoop l_ttbar(ttbarChain);
  l_ttbar.Loop();
  }

  {
  //LM9
  TChain* lm9Chain = new TChain("BasicTreeMaker/tree");
  lm9Chain->Add("~/data/BasicNtuples/V00-00-01/LM9/Spring10/BasicNtuple*");

  basicLoop l_lm9(lm9Chain);
  l_lm9.Loop();
  }

  {
  //LM13
  TChain* lm13Chain = new TChain("BasicTreeMaker/tree");
  lm13Chain->Add("~/data/BasicNtuples/V00-00-01/LM13/Spring10/BasicNtuple*");

  basicLoop l_lm13(lm13Chain);
  l_lm13.Loop();
  }

  {
  //MoreMSSM
  TChain* mmssmChain = new TChain("BasicTreeMaker/tree");
  mmssmChain->Add("~/data/BasicNtuples/V00-00-01/MoreMSSM/Spring10/BasicNtuple*");

  basicLoop l_mmssm(mmssmChain);
  l_mmssm.Loop();
  }

  {
  //QCD-Pt1000toInf-madgraph
  TChain* qcd1000Chain = new TChain("BasicTreeMaker/tree");
  qcd1000Chain->Add("~/data/BasicNtuples/V00-00-01/QCD-Pt1000toInf-madgraph/Spring10/BasicNtuple*");

  basicLoop l_qcd1000(qcd1000Chain);
  l_qcd1000.Loop();
  }

  {
  //QCD-Pt100to250-madgraph
  TChain* qcd100Chain = new TChain("BasicTreeMaker/tree");
  qcd100Chain->Add("~/data/BasicNtuples/V00-00-01/QCD-Pt100to250-madgraph/Spring10/BasicNtuple*");

  basicLoop l_qcd100(qcd100Chain);
  l_qcd100.Loop();
  }

  {
  //QCD-Pt250to500-madgraph
  TChain* qcd250Chain = new TChain("BasicTreeMaker/tree");
  qcd250Chain->Add("~/data/BasicNtuples/V00-00-01/QCD-Pt250to500-madgraph/Spring10/BasicNtuple*");

  basicLoop l_qcd250(qcd250Chain);
  l_qcd250.Loop();
  }

  {
  //QCD-Pt500to1000-madgraph
  TChain* qcd500Chain = new TChain("BasicTreeMaker/tree");
  qcd500Chain->Add("~/data/BasicNtuples/V00-00-01/QCD-Pt500to1000-madgraph/Spring10/BasicNtuple*");

  basicLoop l_qcd500(qcd500Chain);
  l_qcd500.Loop();
  }

  {
  //WJets
  TChain* wjetsChain = new TChain("BasicTreeMaker/tree");
  wjetsChain->Add("~/data/BasicNtuples/V00-00-01/WJets/Spring10/BasicNtuple*");

  basicLoop l_wjets(wjetsChain);
  l_wjets.Loop();
  }

  {
  //ZJets
  TChain* zjetsChain = new TChain("BasicTreeMaker/tree");
  zjetsChain->Add("~/data/BasicNtuples/V00-00-01/ZJets/Spring10/BasicNtuple*");

  basicLoop l_zjets(zjetsChain);
  l_zjets.Loop();
  }

  {
  //Zinvisible
  TChain* znunuChain = new TChain("BasicTreeMaker/tree");
  znunuChain->Add("~/data/BasicNtuples/V00-00-01/Zinvisible/Spring10/BasicNtuple*");

  basicLoop l_znunu(znunuChain);
  l_znunu.Loop();
  }

}
