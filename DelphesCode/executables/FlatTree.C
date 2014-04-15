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
  //  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  //  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  //  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");

  bool verbose =false;
  bool listJetTowers = false;

  SimpleTree tr("simpleTree.root");
  tr.AddVariable("HT");
  tr.AddVariable("njets30");
  tr.AddArray("jetPt",20);

  // Loop over all events
  for(Long64_t entry = 0; entry < numberOfEntries; ++entry) {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
  
    if (verbose||entry%5000==0) cout << "Event " << entry << " / " << numberOfEntries << endl;
    
    float HT = 0;
    int njets30=0;

      // Loop over jets
      for (int i = 0 ; i < branchJet->GetEntries() ; i++) {
	Jet *jet = (Jet*) branchJet->At(i);
	if (jet->PT > 30.) {
	  HT += jet->PT;
	  if (njets30<20) tr.Set("jetPt",njets30,jet->PT);
	  njets30++;

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
	}
      }
      // } // verbose 

      tr.Set("HT",HT);
      tr.Set("njets30",float(njets30));

      tr.Fill();
  } // event

}
