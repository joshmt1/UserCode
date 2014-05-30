/*
root -l examples/Example1.C\(\"delphes_output.root\"\)

or just:

./FlatTree delphes_output.root

N.B. you must touch FlatTree.cpp in order to convince Make that there is anything to do to build this!
*/

//------------------------------------------------------------------------------

void FlatTree(TString inputFile,TString outputFile,const int jobIndex, const int nJobs)
{
  //  gSystem->Load("libDelphes");

  //inputFile can be something like:
  //"root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/upgrade/PhaseI/Configuration0/NoPileUp/tt-4p-600-1100-v1510_14TEV/tt-4p-600-1100-v1510_14TEV*.root"

  cout<<" *** beginning of FlatTree ***"<<endl;
  cout<<"inputFile is: "<<inputFile<<endl;

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);
  

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  // Get pointers to branches used in this analysis
  TClonesArray *branchGenParticles = treeReader->UseBranch("Particle");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchMet = treeReader->UseBranch("MissingET");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  //  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");

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

  //stuff to characterize the mc truth (mostly for signal)
  tr.AddInt("ttbarDecayCode");
  tr.AddInt("SusyProductionMode");
  tr.AddInt("Chi2ToChi1Code");
  tr.AddVariable("genEdgeMll1");
  tr.AddVariable("genEdgeMll2");
  tr.AddInt("nZFromSusy");
  tr.AddInt("nbFromSusy");
  tr.AddInt("ntFromSusy");
  tr.AddInt("nTrueElMu");
  tr.AddInt("nTrueTau");
  tr.AddBool("leptonsMatchChi2ToChi1");
  tr.AddBool("leptonsMatchChi2ToChi1_loose");
  tr.AddBool("leptonsMatchChi2ToChi1_rand");

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
  tr.AddBool("isSF");
  tr.AddVariable("leptonPt1");
  tr.AddVariable("leptonPt2");
  tr.AddVariable("leptonEta1");
  tr.AddVariable("leptonEta2");
  tr.AddInt("leptonFlavor1");
  tr.AddInt("leptonFlavor2");

  tr.AddVariable("mll_loose"); //for edge analysis test
  tr.AddBool("isSF_loose");
  tr.AddInt("leptonFlavor1_loose");
  tr.AddInt("leptonFlavor2_loose");

  tr.AddVariable("mll_rand"); //for edge analysis test
  tr.AddBool("isSF_rand");
  tr.AddInt("leptonFlavor1_rand");
  tr.AddInt("leptonFlavor2_rand");

  //llq inv mass; minimized over the lead two jets
  tr.AddVariable("mllqmin");
  //lq invariant masses (using two leading jets)
  tr.AddArray("mlq",4);
  tr.AddInt("njets30",0);
  tr.AddInt("njets30eta3p0",0);
  tr.AddInt("njets40",0);
  tr.AddInt("njets40eta3p0",0);
  tr.AddInt("nbjets40loose",0);
  tr.AddInt("nbjets40tight",0);
  //leptons for MT2
  tr.AddInt("nElectrons",0);
  tr.AddInt("nMuons",0);
  tr.AddInt("nTaus",0);
  const int MAX_njets=20;
  tr.AddArray("jetPt",MAX_njets); //could add eta, etc
  const int MAX_leptons=4;
  tr.AddArray("electronIso",MAX_leptons); //*not* relative iso
  tr.AddArray("muonIso",MAX_leptons); //*not* relative iso
  tr.AddArray("electronPt",MAX_leptons);
  tr.AddArray("muonPt",MAX_leptons);
  //no need for photons or tighter muons, probably
  //MC truth for DY events
  tr.AddVariable("DYgenPt1");
  tr.AddVariable("DYgenPt2");
  tr.AddInt("DYflavor",-99999);
  tr.AddBool("DYgenRecoMatch1");
  tr.AddBool("DYgenRecoMatch2");

  CrossSections cross_section(inputFile);
  const Long64_t n_events_generated = numberOfEntries; //for now, enforce that the whole sample is run over in one job!

  //count number of 'bad jets' and number of events with bad jets
  Long64_t nbadjets=0,neventsWithBadjets=0;

  TStopwatch loopTimer; //automatically starts the timer

  McTruthInfo geninfo;

  // Loop over selected events
  Long64_t nperjob = numberOfEntries/nJobs;
  Long64_t firstEvent = (jobIndex-1)*nperjob;
  Long64_t lastEvent = nJobs==jobIndex ? numberOfEntries : jobIndex*nperjob;
  cout<<"Running on events from "<<firstEvent+1<<" through "<<lastEvent<<endl;
  for (Long64_t entry = firstEvent; entry < lastEvent; ++entry) {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
  
    if (verbose||((entry-firstEvent)%5000==0 &&entry>0)) cout << "Event " << entry-firstEvent << " / " << lastEvent-firstEvent << endl;

    assert(branchEvent->GetEntries()==1);
    //signal uses HepMCEvent; others use LHEFEvent
    HepMCEvent* evt1 = dynamic_cast<HepMCEvent*>( branchEvent->At(0));
    LHEFEvent* evt2 = dynamic_cast<LHEFEvent*>( branchEvent->At(0));
    double w=0;
    if  (evt1) w=evt1->Weight;
    else if (evt2) w=evt2->Weight;
    else assert(0);
    //    cout<<w<<endl;
    tr.SetDouble("weight",w * cross_section.Get() / n_events_generated); //weight for 1 pb-1
    geninfo.Set(branchGenParticles);
    if (cross_section.GetProcess()==CrossSections::kTop)  tr.SetInt("ttbarDecayCode", geninfo.GetTtbarDecayCode());

     //for debug
    //          geninfo.Dump();
    //    cout<<"nTrue e+mu tau "<<geninfo.countTrueLeptons(McTruthInfo::kElMu)<<" "<<geninfo.countTrueLeptons(McTruthInfo::kTau)<<endl;
    /*
       vector<int> susymoms=    geninfo.findSusyMoms();
    for (size_t iiii=0;iiii<susymoms.size();iiii++) {
      cout<<"SUSY mom = "<<susymoms[iiii]<<endl;
    }
    cout<<" Susy prod code = "<< geninfo.getSusyProductionProcess()<<" ; "<<endl;
    */

    // production code
    if (cross_section.GetProcess()==CrossSections::kSignal) {
      tr.SetInt("SusyProductionMode",   geninfo.getSusyProductionProcess());
      int code =  geninfo.findChi2ToChi1();
      tr.SetInt("Chi2ToChi1Code",  code);
      if (code==1 || code==2) {
	tr.Set("genEdgeMll1",geninfo.getGenMll(1));
	tr.Set("genEdgeMll2",geninfo.getGenMll(2));
      }
      //    cout<<"Chi2ToChi1Code = "<<geninfo.findChi2ToChi1()<<endl;
      //cout<<"nZ = "<< geninfo.findZinSusy()<<endl;
      tr.SetInt("nZFromSusy",    geninfo.findPinSusy(23));
      tr.SetInt("nbFromSusy",    geninfo.findPinSusy(5));
      tr.SetInt("ntFromSusy",    geninfo.findPinSusy(6));
    }
    tr.SetInt("nTrueElMu",geninfo.countTrueLeptons(McTruthInfo::kElMu));
    tr.SetInt("nTrueTau", geninfo.countTrueLeptons(McTruthInfo::kTau));

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
      //good electron
      //store electron iso
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
      //good muon
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

    //redo mll with variations
    MllComputer mllcomp_loose(branchElectron,branchMuon);
    mllcomp_loose.minpt_ = 5;//i'm nervous about the idea of going to 0
    mllcomp_loose.maxetacut_ = 5;
    mllcomp_loose.removegap_ = false;
    tr.Set("mll_loose",mllcomp_loose.GetMll());
    tr.SetBool("isSF_loose",mllcomp_loose.isSF()); //ditto
    TLorentzVector* lepton1_loose = mllcomp_loose.GetLeptonP4(1);
    TLorentzVector* lepton2_loose = mllcomp_loose.GetLeptonP4(2);
    if (lepton1_loose!=0)  tr.SetInt("leptonFlavor1_loose",mllcomp_loose.GetLeptonFlavor(1));
    if (lepton2_loose!=0)  tr.SetInt("leptonFlavor2_loose",mllcomp_loose.GetLeptonFlavor(2));


    MllComputer mllcomp_rand(branchElectron,branchMuon);
    //    mllcomp_rand.minpt_ = 5;//i'm nervous about the idea of going to 0
    //    mllcomp_rand.maxetacut_ = 5;
    //    mllcomp_rand.removegap_ = false;
    mllcomp_rand.randomizeLeptons_=true;
    mllcomp_rand.seed_ = entry+10001;//use a different seed for each event
    tr.Set("mll_rand",mllcomp_rand.GetMll());
    tr.SetBool("isSF_rand",mllcomp_rand.isSF()); //ditto
    TLorentzVector* lepton1_rand = mllcomp_rand.GetLeptonP4(1);
    TLorentzVector* lepton2_rand = mllcomp_rand.GetLeptonP4(2);
    if (lepton1_rand!=0)  tr.SetInt("leptonFlavor1_rand",mllcomp_rand.GetLeptonFlavor(1));
    if (lepton2_rand!=0)  tr.SetInt("leptonFlavor2_rand",mllcomp_rand.GetLeptonFlavor(2));


    //back to the regular mll
    if (lepton1!=0)  tr.Set("leptonPt1",lepton1->Pt());
    if (lepton2!=0)  tr.Set("leptonPt2",lepton2->Pt());
    if (lepton1!=0)  tr.Set("leptonEta1",lepton1->Eta());
    if (lepton2!=0)  tr.Set("leptonEta2",lepton2->Eta());
    if (lepton1!=0)  tr.SetInt("leptonFlavor1",mllcomp.GetLeptonFlavor(1));
    if (lepton2!=0)  tr.SetInt("leptonFlavor2",mllcomp.GetLeptonFlavor(2));

    //if reco leptons and Signal
    if (lepton1!=0 && lepton2!=0 && cross_section.GetProcess()==CrossSections::kSignal) {
      //do these match the leptons from N2->l~ l->N1 
      tr.SetBool("leptonsMatchChi2ToChi1", geninfo.matchesChi2ToChi1Gen(*lepton1,*lepton2,mllcomp.GetLeptonFlavor(1),mllcomp.GetLeptonFlavor(2)));
      //could also check match to SUSY Z; would need additional code

    }
    if (lepton1_loose!=0 && lepton2_loose!=0 && cross_section.GetProcess()==CrossSections::kSignal) {
      //do these match the leptons from N2->l~ l->N1 
      tr.SetBool("leptonsMatchChi2ToChi1_loose", geninfo.matchesChi2ToChi1Gen(*lepton1_loose,*lepton2_loose,mllcomp_loose.GetLeptonFlavor(1),mllcomp_loose.GetLeptonFlavor(2)));
    }
    if (lepton1_rand!=0 && lepton2_rand!=0 && cross_section.GetProcess()==CrossSections::kSignal) {
      //do these match the leptons from N2->l~ l->N1 
      tr.SetBool("leptonsMatchChi2ToChi1_rand", geninfo.matchesChi2ToChi1Gen(*lepton1_rand,*lepton2_rand,mllcomp_rand.GetLeptonFlavor(1),mllcomp_rand.GetLeptonFlavor(2)));
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

	if (jet->PT>30 && jet->Mass < 0) {
	  //	  cout<<" ** Bad jet **, fixing it by setting mass to 0"<<endl;
	  jet->Mass = 0;
	  nbadjets++;
	  if (!badjetinthisevent) { neventsWithBadjets++; badjetinthisevent=true;}
	  //	    cout<<jet->PT<<" "<<jet->Eta<<" "<<jet->Mass<<endl;	    
	}

	if (jet->PT>30 && std::abs(jet->Eta)<2.4) { //mt2 and MHT jets 

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
	if (jet->PT>30 && std::abs(jet->Eta)<3) { 
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

	  //TODO check that I have tight versus loose correct!
	  if ( jet->BTag&1) tr.SetInt("nbjets40loose",1,true);
	  if ( jet->BTag&2) tr.SetInt("nbjets40tight",1,true);

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
