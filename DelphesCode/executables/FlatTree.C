/*
root -l examples/Example1.C\(\"delphes_output.root\"\)

or just:

./FlatTree delphes_output.root

N.B. you must touch FlatTree.cpp in order to convince Make that there is anything to do to build this!
*/

//------------------------------------------------------------------------------

bool isPuJet(Jet * jet,const bool isPhaseII) {
  /* (Dominick)
     It's not in the AN, though I should put it there. Here's the code I use.  The cut values come from:
     https://github.com/sethzenz/Delphes/blob/master/Cards/JetStudies_Phase_II_140PileUp_conf4.tcl
     
     I'm not sure if the values are recommended for general use, but I looked at a few distributions and they seem reasonable.
  */
  // bar: |eta| < 1.5
  // ec: 1.5 <= |eta| < 2.4
  // fwd: 2.4 <= |eta| < 4.0
  // vfwd: |eta| >= 4.0
  
  // values from Jul31
  float cut_beta_bar = 0.13;
  float cut_beta_ec = 0.15;
  float cut_beta_fwd = 0.0;
  float cut_meansqdr_bar = 0.07;
  float cut_meansqdr_ec = 0.07;
  float cut_meansqdr_fwd = 0.07;
  float cut_meansqdr_vfwd = 0.01;
  // apply this only for PhaseII scenario
  if (isPhaseII) {
    cut_beta_fwd = 0.15;
  }

  float eta = jet->Eta;
  float beta = jet->Beta;
  float meansqdr = jet->MeanSqDeltaR;
  if (fabs(eta) < 1.5) {
    if (beta <= cut_beta_bar) return true;
    if (meansqdr >= cut_meansqdr_bar) return true;
  } else if (fabs(eta) < 2.4) {
    if (beta <= cut_beta_ec) return true;
    if (meansqdr >= cut_meansqdr_ec) return true;
  } else if (fabs(eta) < 4.0) {
    if (beta <= cut_beta_fwd) return true;
    if (meansqdr >= cut_meansqdr_fwd) return true;
  } else {
    if (meansqdr >= cut_meansqdr_vfwd) return true;
  }

  return false;

}

void FlatTree(TString inputFile,TString outputFile,const int jobIndex, const int nJobs,const bool doJetLeptonCleaning,const bool usePuJetId)
{
  const bool debugWeights = false;
  //  gSystem->Load("libDelphes");

  DelWeight khWeight;
  khWeight.initialize();
  //inputFile can be something like:
  //"root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/upgrade/PhaseI/Configuration0/NoPileUp/tt-4p-600-1100-v1510_14TEV/tt-4p-600-1100-v1510_14TEV*.root"

  cout<<" *** beginning of FlatTree ***"<<endl;
  cout<<"inputFile is: "<<inputFile<<endl;
  cout<<"Jet/lepton disambiguation is ";
  if (doJetLeptonCleaning) cout<<"enabled"<<endl;
  else cout<<"disabled"<<endl;

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);
  

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchGenParticles = treeReader->UseBranch("Particle");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchMet = usePuJetId ? treeReader->UseBranch("PileUpJetIDMissingET") : treeReader->UseBranch("MissingET");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");

  bool verbose =false;
  //  bool listJetTowers = false;

  //turns out that i don't want this fancy logic.
  //it is easier and more straightforward to handle the naming in the submission script and not here
  //  TString fileEnding;
  //  fileEnding.Form("_%d_%d.root",jobIndex,nJobs);
  //  assert(outputFile.EndsWith(".root"));
  //  outputFile.ReplaceAll(".root",fileEnding);

  SimpleTree tr(outputFile);
  //bookkeeping
  tr.AddDouble("weight");
  tr.AddVariable("kfactor",1); //already included in 'weight' but keep it just for completeness

  //stuff to characterize the mc truth (mostly for signal)
  tr.AddInt("ttbarDecayCode");
  tr.AddVariable("ttbarGenMll");
  tr.AddVariable("genMet",-99);
  tr.AddVariable("genMetPhi",-99);
  tr.AddInt("genNneutrinos",-99);
  tr.AddInt("SusyProductionMode");
  tr.AddInt("Chi2ToChi1Code");
  tr.AddInt("Chi4ToChi1Code");
  tr.AddVariable("genEdgeMll1");
  tr.AddVariable("genEdgeMll2");
  tr.AddVariable("genEdgeLepPt1");
  tr.AddVariable("genEdgeLepPt2");

  tr.AddInt("nZFromSusy");
  tr.AddInt("nbFromSusy");
  tr.AddInt("ntFromSusy");
  tr.AddInt("nTrueElMu");
  tr.AddInt("nTrueTau");
  tr.AddArray("genLepPt",2);
  tr.AddArray("genLepEta",2);
  tr.AddArray("genLepFlavor",2);
  tr.AddArray("genLepRecoIso",2);
  tr.AddBool("leptonsMatchChi2ToChi1");
  tr.AddBool("leptonsMatchChi2ToChi1_veryloose");
  tr.AddBool("leptonsMatchChi2ToChi1_loose");
  tr.AddBool("leptonsMatchChi4ToChi1_loose");

  //jet observables
  tr.AddVariable("HT",0);
  tr.AddVariable("MET");
  tr.AddVariable("METphi");
  tr.AddVariable("MHT");
  tr.AddVariable("MHTphi");
  tr.AddVariable("VSPT");
  tr.AddVariable("MT2");
  tr.AddVariable("minDeltaPhi");
  tr.AddVariable("mll"); //for edge analysis
  tr.AddVariable("mll_maxEta");//max eta of the OS dileptons
  tr.AddVariable("MT_l3MET",-99);
  /*
notes on this variable: right now it uses the same pT threshold as two edge leptons
in reality one would certainly want to allow at least two different thresholds when
doing a 3l analysis. This would require reworking the code.

Also, prob need to store the flavor and charge of this lepton if we're really going to use this
   */
  tr.AddBool("isSF");
  tr.AddVariable("leptonPt1");
  tr.AddVariable("leptonPt2");
  tr.AddVariable("leptonEta1");
  tr.AddVariable("leptonEta2");
  tr.AddInt("leptonFlavor1");
  tr.AddInt("leptonFlavor2");
  tr.AddVariable("leptonIso1");//relative isolation
  tr.AddVariable("leptonIso2");

  tr.AddVariable("mll_veryloose"); //for edge analysis test
  tr.AddBool("isSF_veryloose");
  tr.AddInt("leptonFlavor1_veryloose");
  tr.AddInt("leptonFlavor2_veryloose");
  tr.AddVariable("mll_maxEta_veryloose");//max eta of the OS dileptons
  tr.AddVariable("leptonIso1_veryloose");//relative isolation
  tr.AddVariable("leptonIso2_veryloose");

  tr.AddVariable("mll_loose"); //for edge analysis test
  tr.AddBool("isSF_loose");
  tr.AddVariable("leptonPt1_loose");
  tr.AddVariable("leptonPt2_loose");
  tr.AddInt("leptonFlavor1_loose");
  tr.AddInt("leptonFlavor2_loose");
  tr.AddVariable("mll_maxEta_loose");//max eta of the OS dileptons
  tr.AddVariable("leptonIso1_loose");//relative isolation
  tr.AddVariable("leptonIso2_loose");
  tr.AddVariable("MT_l3MET_loose",-99);

  //llq inv mass; minimized over the lead two jets
  tr.AddVariable("mllqmin");
  //lq invariant masses (using two leading jets)
  tr.AddArray("mlq",4);
  tr.AddInt("nPUjets30",0);
  tr.AddInt("njets30",0);
  tr.AddInt("njets30eta3p0",0);
  tr.AddInt("njets40",0);
  tr.AddInt("njets40eta3p0",0);
  tr.AddInt("nbjets40loose",0);
  tr.AddInt("nbjets40medium",0);
  tr.AddInt("nbjets40tight",0);
  //leptons for MT2
  tr.AddInt("nElectrons",0);
  tr.AddInt("nMuons",0);
  tr.AddInt("nTaus",0);
  const int MAX_njets=20;
  tr.AddArray("jetPt",MAX_njets); //could add eta, etc
  const int MAX_leptons=4;
  tr.AddArray("electronIso",MAX_leptons); //these *are* RelIso, I think
  tr.AddArray("muonIso",MAX_leptons);
  tr.AddArray("electronPt",MAX_leptons); 
  tr.AddArray("muonPt",MAX_leptons);
  //no need for photons or tighter muons, probably
  //MC truth for DY events
  tr.AddVariable("DYgenPt1");
  tr.AddVariable("DYgenPt2");
  tr.AddInt("DYflavor",-99999);
  tr.AddBool("DYgenRecoMatch1");
  tr.AddBool("DYgenRecoMatch2");

  //get a full file name instead of just path/*.root
  TString a_full_file_name =  chain.GetListOfFiles()->At(0)->GetTitle();

  CrossSections cross_section(a_full_file_name);
  //whole chain needs to be available in order to get the normalization right
  const Long64_t n_events_generated = numberOfEntries; 

  //count number of 'bad jets' and number of events with bad jets
  Long64_t nbadjets=0,neventsWithBadjets=0;

  TStopwatch loopTimer; //automatically starts the timer

  McTruthInfo geninfo;
  //set k-factors for signal
  //code from stefan depends on seeing the full path (or really the file name) of the sample
  //so make sure to pass a full path including file name instead of path/*.root
  if (cross_section.GetProcess()==CrossSections::kSignal) {
    geninfo.setNloKFactor_fromfileName(a_full_file_name);
  }


  // Loop over selected events
  Long64_t nperjob = numberOfEntries/nJobs;
  Long64_t firstEvent = (jobIndex-1)*nperjob;
  Long64_t lastEvent = nJobs==jobIndex ? numberOfEntries : jobIndex*nperjob;
  cout<<"Running on events from "<<firstEvent+1<<" through "<<lastEvent<<endl;
  for (Long64_t entry = firstEvent; entry < lastEvent; ++entry) {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
  
    if (verbose||((entry-firstEvent)%5000==0 &&entry>0)) cout << "Event " << entry-firstEvent << " / " << lastEvent-firstEvent << endl;

    //disambiguate leptons and jets
    if (doJetLeptonCleaning)    JetLeptonCleaner(branchJet,branchElectron,branchMuon,branchPhoton);

    //MC truth stuff for signal
    double kfactor=1;
    geninfo.Set(branchGenParticles);
    if (cross_section.GetProcess()==CrossSections::kSignal) {
      int prodcode=geninfo.getSusyProductionProcess();
      kfactor=geninfo.getNloKFactor(prodcode);
      tr.Set("kfactor",(float)kfactor);
      tr.SetInt("SusyProductionMode", prodcode);
      int code =  geninfo.findChi2ToChi1();
      tr.SetInt("Chi2ToChi1Code",  code);
      if (code==1 || code==2) {
	tr.Set("genEdgeMll1",geninfo.getGenMll(1));
	tr.Set("genEdgeMll2",geninfo.getGenMll(2));
	tr.Set("genEdgeLepPt1",geninfo.getGenEdgeLeptonPt(1));
	tr.Set("genEdgeLepPt2",geninfo.getGenEdgeLeptonPt(2));
      }
      int code4 =  geninfo.findChi2ToChi1(4);//finds chi4->chi1
      tr.SetInt("Chi4ToChi1Code",  code4);
      //cout<<"nZ = "<< geninfo.findZinSusy()<<endl;
      tr.SetInt("nZFromSusy",    geninfo.findPinSusy(23));
      tr.SetInt("nbFromSusy",    geninfo.findPinSusy(5));
      tr.SetInt("ntFromSusy",    geninfo.findPinSusy(6));
    }


    //signal uses HepMCEvent; others use LHEFEvent
    HepMCEvent* evt1=0; LHEFEvent* evt2=0;
    double w=0;
    if (branchEvent) {
      evt1 = dynamic_cast<HepMCEvent*>( branchEvent->At(0));
      evt2 = dynamic_cast<LHEFEvent*>( branchEvent->At(0));
      if  (evt1) w=evt1->Weight;
      else if (evt2) w=evt2->Weight;
      else assert(0);
      //      cout<<" w [ntuple] = "<<w<<endl;
    }
    else w=1;
    //here's the story:
    //for SM samples it is supposed to be important to use the evt->Weight
    //then for signal I found that the underlying class is different, hence the
    //dynamic cast above
    //But then I discovered that for signal, the evt->Weight seems to be a bogus constant value
    //so this if statement replaces it with 1.
    //so obviously i could streamline this and get rid of the dynamic cast, but for now
    //let's just leave this duplicate logic in place
    //in case of signal, override the value from the ntuple no matter what
    if (cross_section.GetProcess()==CrossSections::kSignal) w = 1;

    if ( (cross_section.GetProcess()==CrossSections::kTop || cross_section.GetProcess()==CrossSections::kBoson) && (branchEvent==0|| debugWeights)) {
      //this if block is only here for the buggy samples with no branchEvent. keep it for posterity but let's make sure it does not run
      assert(0);

      int code = (cross_section.GetProcess()==CrossSections::kTop) ? 1 : 2;
      //need  std::vector<GenParticle> GenParticles
      vector<GenParticle> gparticles;
      for (int k=0; k<branchGenParticles->GetEntries();++k) {
	GenParticle * c =(GenParticle*) branchGenParticles->At(k);
	if (c!=0) 	  gparticles.push_back( *c);
      }
      if (debugWeights) cout<<"w [ ntuple ] = "<<w<<endl;
      w= khWeight.weight(code,gparticles);
      if (debugWeights) cout<<"w [   ken  ] = "<<w<<endl;
    }    //otherwise, keep the default value of 1 from the else block above

    tr.SetDouble("weight",w * kfactor * cross_section.Get() / n_events_generated); //weight for 1 pb-1
    //geninfo.Dump();
    if (cross_section.GetProcess()==CrossSections::kTop) {
      float ttgenmll;
      const      int ttdc=geninfo.GetTtbarDecayCode(ttgenmll);
      tr.SetInt("ttbarDecayCode", ttdc);
      tr.Set("ttbarGenMll",ttgenmll);
    }


    //more MC truth
    tr.SetInt("nTrueElMu",geninfo.countTrueLeptons(McTruthInfo::kElMu));
    tr.SetInt("nTrueTau", geninfo.countTrueLeptons(McTruthInfo::kTau));
    for (unsigned int ilep = 0; ilep<geninfo.getGenLeptons().size();ilep++) {
      if (ilep>=2) break;

      tr.Set("genLepPt",(int)ilep,(float) geninfo.getGenLeptons().at(ilep)->PT);
      tr.Set("genLepEta",(int)ilep,(float)geninfo.getGenLeptons().at(ilep)->Eta);
      tr.Set("genLepFlavor",(int)ilep,std::abs(geninfo.getGenLeptons().at(ilep)->PID));
      tr.Set("genLepRecoIso",(int)ilep, geninfo.getIsolationOfMatch(ilep,branchElectron,branchMuon));
    }
    std::pair<float,float> gmet = geninfo.getGenMet();
    tr.Set("genMet",gmet.first);
    tr.Set("genMetPhi",gmet.second);
    tr.SetInt("genNneutrinos",geninfo.countNeutrinos());

    //DY truth
    if (cross_section.GetProcess()==CrossSections::kBoson) {
      std::vector< std::pair< TLorentzVector, int> > gendy=geninfo.GetDYTruth() ;
      vector<int> dymatches =  geninfo.MatchDYRecoGen( gendy,branchElectron,branchMuon);
      if ( gendy.size()==2 ) {
	tr.Set("DYgenPt1",gendy[0].first.Pt());
	tr.Set("DYgenPt2",gendy[1].first.Pt());
	tr.SetInt("DYflavor",std::abs(gendy[0].second));
	tr.SetBool("DYgenRecoMatch1",dymatches[0]>=0);
	tr.SetBool("DYgenRecoMatch2",dymatches[1]>=0);
      }
    }

    //store MET
    assert( branchMet->GetEntries() ==1); //sanity
    MissingET * met = (MissingET*) branchMet->At(0);
    tr.Set("MET",met->MET);
    tr.Set("METphi",met->Phi);
    
    float MHTx=0,MHTy=0,MSTx=0,MSTy=0;

    double pa[3]={0,0,0};
    double pb[3]={0,0,0};
    double pmiss[3]={0,0,0};
    pmiss[1] = met->MET*cos(met->Phi);
    pmiss[2] = met->MET*sin(met->Phi);

    //count electrons, muons
    for (int i = 0 ; i < branchElectron->GetEntries() ; i++) {
      Electron *el = (Electron*) branchElectron->At(i);
      if (el->PT < 10 ) continue;
      if ( std::abs(el->Eta) > 2.4) continue; 
      if (el->IsolationVar >0.2) continue; //add isolation cut
      //good electron
      if ( tr.GetInt("nElectrons") < MAX_leptons) {
	tr.Set("electronIso",tr.GetInt("nElectrons"),el->IsolationVar);
	tr.Set("electronPt",tr.GetInt("nElectrons"),el->PT);
      }
      tr.SetInt("nElectrons",1,true);
      MSTx -= el->PT * cos(el->Phi);
      MSTy -= el->PT * sin(el->Phi);
    }
    for (int i = 0 ; i < branchMuon->GetEntries() ; i++) {
      Muon *mu = (Muon*) branchMuon->At(i);
      if (mu->PT < 10 ) continue;
      if ( std::abs(mu->Eta) > 2.4) continue; 
      if (mu->IsolationVar > 0.2) continue;
      //good muon
      if ( tr.GetInt("nMuons") < MAX_leptons) {
	tr.Set("muonIso",tr.GetInt("nMuons"),mu->IsolationVar);
	tr.Set("muonPt",tr.GetInt("nMuons"),mu->PT);
      }
      tr.SetInt("nMuons",1,true);
      MSTx -= mu->PT * cos(mu->Phi);
      MSTy -= mu->PT * sin(mu->Phi);
    }

    //OSDL Edge variables
    MllComputer mllcomp(branchElectron,branchMuon);
    tr.Set("mll",mllcomp.GetMll());
    tr.Set("mll_maxEta",mllcomp.GetMaxEta()); //must be called after GetMll()
    tr.SetBool("isSF",mllcomp.isSF()); //ditto
    TLorentzVector* lepton1 = mllcomp.GetLeptonP4(1);
    TLorentzVector* lepton2 = mllcomp.GetLeptonP4(2);
    if ( mllcomp.GetNExtraLeptons() >0) {
      TLorentzVector pl = mllcomp.GetExtraLeptonP4(0);
      tr.Set("MT_l3MET",sqrt(2*met->MET*pl.Pt()*(1-cos(Util::DeltaPhi(pl.Phi(),met->Phi)))));
    }

    //redo mll with variations
    MllComputer mllcomp_veryloose(branchElectron,branchMuon);
    //very loose == unrealistically loose -- no eta cuts and 5 GeV pT threshold
    mllcomp_veryloose.minpt_ = 5;//i'm nervous about the idea of going to 0
    mllcomp_veryloose.maxetacut_ = 5;
    mllcomp_veryloose.removegap_ = false;
    tr.Set("mll_veryloose",mllcomp_veryloose.GetMll());
    tr.SetBool("isSF_veryloose",mllcomp_veryloose.isSF()); //ditto
    TLorentzVector* lepton1_veryloose = mllcomp_veryloose.GetLeptonP4(1);
    TLorentzVector* lepton2_veryloose = mllcomp_veryloose.GetLeptonP4(2);
    if (lepton1_veryloose!=0)  tr.SetInt("leptonFlavor1_veryloose",mllcomp_veryloose.GetLeptonFlavor(1));
    if (lepton2_veryloose!=0)  tr.SetInt("leptonFlavor2_veryloose",mllcomp_veryloose.GetLeptonFlavor(2));
    if (lepton1_veryloose!=0)  tr.Set("leptonIso1_veryloose",mllcomp_veryloose.GetLeptonIsolation(1));
    if (lepton2_veryloose!=0)  tr.Set("leptonIso2_veryloose",mllcomp_veryloose.GetLeptonIsolation(2));
    tr.Set("mll_maxEta_veryloose",mllcomp_veryloose.GetMaxEta()); //must be called after GetMll()

    MllComputer mllcomp_loose(branchElectron,branchMuon);
    //loose == stanndard eta cuts but loosen pT to 10
    mllcomp_loose.minpt_ = 10;//
    tr.Set("mll_loose",mllcomp_loose.GetMll());
    tr.SetBool("isSF_loose",mllcomp_loose.isSF()); //ditto
    TLorentzVector* lepton1_loose = mllcomp_loose.GetLeptonP4(1);
    TLorentzVector* lepton2_loose = mllcomp_loose.GetLeptonP4(2);
    if (lepton1_loose!=0)  tr.SetInt("leptonFlavor1_loose",mllcomp_loose.GetLeptonFlavor(1));
    if (lepton2_loose!=0)  tr.SetInt("leptonFlavor2_loose",mllcomp_loose.GetLeptonFlavor(2));
    if (lepton1_loose!=0)  tr.Set("leptonPt1_loose",lepton1_loose->Pt());
    if (lepton2_loose!=0)  tr.Set("leptonPt2_loose",lepton2_loose->Pt());
    if (lepton1_loose!=0)  tr.Set("leptonIso1_loose",mllcomp_loose.GetLeptonIsolation(1));
    if (lepton2_loose!=0)  tr.Set("leptonIso2_loose",mllcomp_loose.GetLeptonIsolation(2));
    tr.Set("mll_maxEta_loose",mllcomp_loose.GetMaxEta()); //must be called after GetMll()
    if ( mllcomp_loose.GetNExtraLeptons() >0) {
      TLorentzVector pl = mllcomp_loose.GetExtraLeptonP4(0);
      tr.Set("MT_l3MET_loose",sqrt(2*met->MET*pl.Pt()*(1-cos(Util::DeltaPhi(pl.Phi(),met->Phi)))));
    }

    //back to the regular mll
    if (lepton1!=0)  tr.Set("leptonPt1",lepton1->Pt());
    if (lepton2!=0)  tr.Set("leptonPt2",lepton2->Pt());
    if (lepton1!=0)  tr.Set("leptonEta1",lepton1->Eta());
    if (lepton2!=0)  tr.Set("leptonEta2",lepton2->Eta());
    if (lepton1!=0)  tr.SetInt("leptonFlavor1",mllcomp.GetLeptonFlavor(1));
    if (lepton2!=0)  tr.SetInt("leptonFlavor2",mllcomp.GetLeptonFlavor(2));
    if (lepton1!=0)  tr.Set("leptonIso1",mllcomp.GetLeptonIsolation(1));
    if (lepton2!=0)  tr.Set("leptonIso2",mllcomp.GetLeptonIsolation(2));

    //if reco leptons and Signal
    if (lepton1!=0 && lepton2!=0 && cross_section.GetProcess()==CrossSections::kSignal) {
      //do these match the leptons from N2->l~ l->N1 
      tr.SetBool("leptonsMatchChi2ToChi1", geninfo.matchesChi2ToChi1Gen(*lepton1,*lepton2,mllcomp.GetLeptonFlavor(1),mllcomp.GetLeptonFlavor(2)));
      //could also check match to SUSY Z; would need additional code

    }
    if (lepton1_veryloose!=0 && lepton2_veryloose!=0 && cross_section.GetProcess()==CrossSections::kSignal) {
      //do these match the leptons from N2->l~ l->N1 
      tr.SetBool("leptonsMatchChi2ToChi1_veryloose", geninfo.matchesChi2ToChi1Gen(*lepton1_veryloose,*lepton2_veryloose,mllcomp_veryloose.GetLeptonFlavor(1),mllcomp_veryloose.GetLeptonFlavor(2)));
    }
    if (lepton1_loose!=0 && lepton2_loose!=0 && cross_section.GetProcess()==CrossSections::kSignal) {
      //do these match the leptons from N2->l~ l->N1 
      tr.SetBool("leptonsMatchChi2ToChi1_loose", geninfo.matchesChi2ToChi1Gen(*lepton1_loose,*lepton2_loose,mllcomp_loose.GetLeptonFlavor(1),mllcomp_loose.GetLeptonFlavor(2)));
      //do these match the leptons from N4      ->N1
      tr.SetBool("leptonsMatchChi4ToChi1_loose",geninfo.matchesChi2ToChi1Gen(*lepton1_loose,*lepton2_loose,mllcomp_loose.GetLeptonFlavor(1),mllcomp_loose.GetLeptonFlavor(2),4) );
    }

    //jets to go into the MT2 calculation
    vector<float> mt2jets_px;
    vector<float> mt2jets_py;
    vector<float> mt2jets_pz;
    vector<float> mt2jets_e;

    float mllq1=1e9,mllq2=1e9;
    float minDeltaPhi=1e9;
    bool badjetinthisevent=false;
      // Loop over jets
      for (int i = 0 ; i < branchJet->GetEntries() ; i++) {
	Jet *jet = (Jet*) branchJet->At(i);

	if (jet->PT<30) continue;//do not use jets with pt<30 for anything

	if ( usePuJetId && isPuJet(jet,true) ) { //2nd argument is whether or not the sample is phase II or not
	  tr.SetInt("nPUjets30",1,true);
	  continue;
	}

	if ( jet->Mass < 0) {
	  //	  cout<<" ** Bad jet **, fixing it by setting mass to 0"<<endl;
	  jet->Mass = 0;
	  nbadjets++;
	  if (!badjetinthisevent) { neventsWithBadjets++; badjetinthisevent=true;}
	  //	    cout<<jet->PT<<" "<<jet->Eta<<" "<<jet->Mass<<endl;	    
	}

	if (std::abs(jet->Eta)<2.4) { //mt2 and MHT jets  (pt>30)

	  TLorentzVector jjj = jet->P4();
	  mt2jets_px.push_back( jjj.Px());
	  mt2jets_py.push_back( jjj.Py());
	  mt2jets_pz.push_back( jjj.Pz());
	  mt2jets_e.push_back(  jjj.E());

	  //my nomenclaure -- MHT uses only jets, MST uses jets and leptons
	  MHTx -= jet->PT * cos(jet->Phi);
	  MHTy -= jet->PT * sin(jet->Phi);
	  MSTx -= jet->PT * cos(jet->Phi);
	  MSTy -= jet->PT * sin(jet->Phi);

	  //use 30 GeV threshold for taus as well ; was 20 in the real analysis
	  if ( jet->TauTag == 1)  tr.SetInt("nTaus",1,true);
	}


	//eta 3.0 jets: used for edge analysis and for HT calculation in MT2 analysis
	if (std::abs(jet->Eta)<3) {  //pt>30
	  tr.SetInt("njets30eta3p0",1,true);

	  if (lepton1!=0 && lepton2!=0) {
	    int nj=tr.GetInt("njets30eta3p0");
	    if      (nj == 1) {
	      mllq1 = ((*lepton1)+(*lepton2)+jet->P4()).M();
	      tr.Set("mlq",0, ((*lepton1)+jet->P4()).M());
	      tr.Set("mlq",1, ((*lepton2)+jet->P4()).M());
	    }
	    else if (nj == 2) {
	      mllq2 = ((*lepton1)+(*lepton2)+jet->P4()).M();
	      tr.Set("mlq",2, ((*lepton1)+jet->P4()).M());
	      tr.Set("mlq",3, ((*lepton2)+jet->P4()).M());
	    }
	  }

	  if (jet->PT>40) {
	    tr.SetInt("njets40eta3p0",1,true);
	    if (jet->PT>50) { //HT calculation for MT2
	      tr.Set("HT",jet->PT,true); 
	    }
	  }
	}

	//"all purpose" jets e.g. b-tagging
	if (jet->PT > 40. && std::abs(jet->Eta)<2.4) {
	  if ( tr.GetInt("njets40")<MAX_njets) tr.Set("jetPt",tr.GetInt("njets40"),jet->PT);
	  tr.SetInt("njets40",1,true);

	  //in Jul28 production, switch to loose/medium/tight
	  if ( jet->BTag&1) tr.SetInt("nbjets40loose",1,true);
	  if ( jet->BTag&2) tr.SetInt("nbjets40medium",1,true);
	  if ( jet->BTag&4) tr.SetInt("nbjets40tight",1,true);

	  //use lead 4 jets for mindeltaphi
	  if (tr.GetInt("njets40") <=4 && Util::DeltaPhi(jet->Phi,met->Phi) < minDeltaPhi) minDeltaPhi = Util::DeltaPhi(jet->Phi,met->Phi);


	} //jet pT 40
      } //loop over jets
      // } // verbose 

      float mllq = mllq1<mllq2 ? mllq1 : mllq2;
      tr.Set("mllqmin",mllq);

      tr.Set("MHT",sqrt(MHTx*MHTx + MHTy*MHTy));
      tr.Set("MHTphi",atan2(MHTy,MHTx)); 
      tr.Set("minDeltaPhi",minDeltaPhi);

      //VSPT calculation
      double MST = sqrt(MSTx*MSTx + MSTy*MSTy);
      double MSTphi = atan2(MSTy,MSTx);
      TVector3 METvec( met->MET*cos(met->Phi),met->MET*sin(met->Phi) ,0);
      TVector3 MSTvec( MST*cos(MSTphi),MST*sin(MSTphi),0);
      TVector3 vdiff = METvec-MSTvec;
      tr.Set("VSPT",vdiff.Mag());

      tr.SetInt("njets30",(int)mt2jets_px.size());

      if (mt2jets_px.size()>=2) {
// 	FindHemispheres hemifinder(branchJet);
// 	pa[1] = hemifinder.Px(1);
// 	pa[2] = hemifinder.Py(1);
// 	pb[1] = hemifinder.Px(2);
// 	pb[2] = hemifinder.Py(2);

	Hemisphere ethhemi(mt2jets_px,mt2jets_py,mt2jets_pz,mt2jets_e,2,3); //final arguments are options (inv mass seeds, min Lund distance)
	pa[1] = ethhemi.getAxis1()[0] * ethhemi.getAxis1()[3];
	pa[2] = ethhemi.getAxis1()[1] * ethhemi.getAxis1()[3];
	pb[1] = ethhemi.getAxis2()[0] * ethhemi.getAxis2()[3];
	pb[2] = ethhemi.getAxis2()[1] * ethhemi.getAxis2()[3];
	mt2_bisect::mt2 mt2_event_eth;
	mt2_event_eth.set_momenta( pa, pb, pmiss );
	mt2_event_eth.set_mn( 0 ); //LSP mass
	float  mt2_value_eth = (float) mt2_event_eth.get_mt2();
	tr.Set("MT2", mt2_value_eth);

      }
      else {tr.Set("MT2", -99);}//tr.Set("MT2jmt", -1);}

      tr.Fill();

  } // event

  loopTimer.Stop();

  cout<<"-- Number of 'bad jets' (30 GeV) = "<<nbadjets<<endl;
  cout<<"-- number of events with a 'bad jet' (30 GeV) = "<<neventsWithBadjets<<endl;

  cout<<"events / CPU (Wall): "<<lastEvent-firstEvent <<" / "<<loopTimer.CpuTime()<<" ("<<loopTimer.RealTime()<<") sec = "<<
    (lastEvent-firstEvent)/loopTimer.CpuTime()<<" ("<<(lastEvent-firstEvent)/loopTimer.RealTime()<<") Hz"<<endl;

}

