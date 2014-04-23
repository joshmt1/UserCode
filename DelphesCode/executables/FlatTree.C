/*
root -l examples/Example1.C\(\"delphes_output.root\"\)

or just:

./FlatTree delphes_output.root

N.B. you must touch FlatTree.cpp in order to convince Make that there is anything to do to build this!
*/

//------------------------------------------------------------------------------

void FlatTree(const char *inputFile)
{
  //  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);
  
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  // Get pointers to branches used in this analysis
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchMet = treeReader->UseBranch("MissingET");
  //  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  //  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  //  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");

  bool verbose =false;
  bool listJetTowers = false;

  SimpleTree tr("simpleTree.root");
  tr.AddVariable("HT",0);
  tr.AddVariable("MET");
  tr.AddVariable("METphi");
  tr.AddVariable("MHT");
  tr.AddVariable("MHTphi");
  tr.AddVariable("VSPT");
  tr.AddVariable("MT2");
  tr.AddInt("njets20",0);
  tr.AddInt("njets40",0);
  tr.AddInt("nbjets40loose",0);
  tr.AddInt("nbjets40tight",0);
  tr.AddArray("jetPt",20);

  //i think we veto on these things
  tr.AddInt("nElectrons");
  tr.AddInt("nMuons");
  tr.AddInt("nTaus");
  //no need for photons or tighter muons, probably


  // Loop over all events
  for(Long64_t entry = 0; entry < numberOfEntries; ++entry) {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
  
    if (verbose||entry%5000==0) cout << "Event " << entry << " / " << numberOfEntries << endl;

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

    //TODO count electrons, muons, taus

    //also calculate MST

    //jets to go into the MT2 calculation
    vector<float> mt2jets_px;
    vector<float> mt2jets_py;
    vector<float> mt2jets_pz;
    vector<float> mt2jets_e;

    float minDeltaPhi=1e9;
      // Loop over jets
      for (int i = 0 ; i < branchJet->GetEntries() ; i++) {
	Jet *jet = (Jet*) branchJet->At(i);

	//TODO add eta, quality cuts

	if (jet->PT>20 && std::abs(jet->Eta)<2.4) { //mt2 and MHT jets ; TODO quality cuts
	  mt2jets_px.push_back( jet->PT * cos(jet->Phi));
	  mt2jets_py.push_back( jet->PT * sin(jet->Phi));
	  TLorentzVector jjj = jet->P4(); //let ROOT do the math for me
	  mt2jets_pz.push_back( jjj.Pz());
	  mt2jets_e.push_back(  jjj.E());

	  //my nomenclaure -- MHT uses only jets, MST uses jets and leptons
	  MHTx -= jet->PT * cos(jet->Phi);
	  MHTy -= jet->PT * sin(jet->Phi);
	  MSTx -= jet->PT * cos(jet->Phi);
	  MSTy -= jet->PT * sin(jet->Phi);
	}

	//HT calculation ; TODO quality cuts
	if (jet->PT>50 && std::abs(jet->Eta)<3) {
	  tr.Set("HT",jet->PT,true);
	}

	//"all purpose" jets e.g. b-tagging TODO quality
	if (jet->PT > 40. && std::abs(jet->Eta)<2.4) {
	  if ( tr.GetInt("njets40")<20) tr.Set("jetPt",tr.GetInt("njets40"),jet->PT);
	  tr.SetInt("njets40",1,true);

	  //TODO check that I have tight versus loose correct!
	  if ( jet->BTag&1) tr.SetInt("nbjets40loose",1,true);
	  if ( jet->BTag&2) tr.SetInt("nbjets40tight",1,true);

	  //use lead 4 jets for mindeltaphi
	  if (tr.GetInt("njets40") <=4 && Util::DeltaPhi(jet->Phi,met->Phi) < minDeltaPhi) minDeltaPhi = Util::DeltaPhi(jet->Phi,met->Phi);

/*
	  cout << "  Jet " << i << endl;
	  cout << "    pT: " << jet->PT << endl;
	  cout << "    Eta: " << jet->Eta << endl;
	  cout << "    BTag: " << bool(jet->BTag&1) << " | " << bool(jet->BTag&2) << endl;
          cout << "    TauTag: " << jet->TauTag << endl;
	  //	  cout << "    Area: " << jet->AreaP4().Pt() << endl;
	  cout << "    Jet Pileup ID" << endl;
	  cout << "      Beta*: " << jet->BetaStar << endl;
	  cout << "      Fractional pT in annuli (<0.1, 0.1-0.2, ..., 0.4-0.5) " << jet->FracPt[0] << " " << jet->FracPt[1] << " " << jet->FracPt[2] << " " << jet->FracPt[3] << " " << jet->FracPt[4] << endl;
	  cout << "      <dR^2>: " << jet->MeanSqDeltaR << endl;
          cout << "      NNeutrals: " << jet->NNeutrals << endl;
	  cout << "      NCharged: " << jet->NCharged << endl;
	  cout << "    Number of constituents: " << jet->Constituents.GetEntries() << endl;
*/
// 	  if (listJetTowers && branchEFlowTrack && branchEFlowTower && branchEFlowMuon) {
// 	    for (int j = 0 ; j < jet->Constituents.GetEntries() ; j++) {
// 	      TObject *obj = jet->Constituents[j];
// 	      if (obj && obj->IsA() == Tower::Class()) {
// 		Tower *tow = static_cast<Tower *> ( obj ) ;
// 		cout << "     Jet constituent Et Eta Phi Time (at calo) " << tow->ET << " " << tow->Eta << " " << tow->Phi << " " << tow->TOuter << endl;
// 	      } else {
// 		//		cout << "  not a tower - could check if it's a track instead (cf. Example3.C)" << endl;
// 	      }
// 	    }
// 	  } 
	} //jet pT 40
      } //loop over jets
      // } // verbose 

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

      tr.SetInt("njets20",(int)mt2jets_px.size());

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



}
