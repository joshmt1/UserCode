#define basicLoop_cxx
#include "basicLoop.h"
#include <TH1.h>
#include <TH2.h>
#include <TDatime.h>
#include <TFile.h>
#include <TStyle.h>
#include <fstream>

// just leave this code here, untouched.
// it can be copied and pasted as an event loop template
void basicLoop::exampleLoop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries(); //jmt: remove Fast

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if (Cut(ientry) < 0) continue; //jmt use cut

   }
}

/*
print a cut flow table
lumi is set in basicLoop.h
*/
void basicLoop::cutflow()
{

  std::vector<int> npass;
  
  if (fChain == 0) return;
  
  const   double sigma = getCrossSection(findInputName());
  
  Long64_t nentries = fChain->GetEntries(); //jmt: remove Fast
  
  const   double   weight = lumi * sigma / double(nentries); //calculate weight
  
  Long64_t nbytes = 0, nb = 0;
  
  LoadTree(0);
  nb = fChain->GetEntry(0);   nbytes += nb;
  for (unsigned int i=0 ; i<cutResults->size(); i++) {
    npass.push_back(0);
  }
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    for (unsigned int i=0 ; i<cutResults->size(); i++) {
      if (cutRequired(i) && cutResults->at(i) )   npass.at(i) = npass.at(i) +1;
      else if (cutRequired(i) && !cutResults->at(i) ) break;
    }
    
  }

  TString samplename=  getSampleName(findInputName());
  TString filename="cutflow_"; 
  filename+=CutSchemeNames_[theCutScheme_]; filename+=".";
  filename+=samplename;
  filename+=".dat";
  ofstream file(filename.Data());
  
  for (unsigned int i=0 ; i<npass.size(); i++) {
    
    if (cutRequired(i)) {
      
      //error on n is sqrt n
      double error = sqrt(npass.at(i));
      double weighted = npass.at(i) * weight;
      double weighted_error = error*weight;
      
      char ccc[150];
      sprintf(ccc,"%20s %15d | Weighted = %f +/- %f",cutnames_.at(i).Data(),npass.at(i),weighted,weighted_error);
      cout<<ccc<<endl;

      file <<  weighted<<"\t" << weighted_error<<endl;
    }
  }
  
  file.close();
}

/*
ABCD tree maker
*/
void basicLoop::ABCDtree()
{
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntries(); //jmt: remove Fast

  TString inname=findInputName(); //uses fChain
  std::cout<<"Got an input file name as: "<<inname<<std::endl;
  
  double sigma = getCrossSection(inname);
  TString sampleName = getSampleName(inname);
  if (sigma<=0) return;
  
  if (nentries == 0) {std::cout<<"Chain has no entries!"<<std::endl; return;}
  
  double   weight = lumi * sigma / double(nentries); //calculate weight
  
  //open output file
  TString outfilename="ABCDtree."; outfilename+=sampleName; outfilename+=".root";
  TFile fout(outfilename,"RECREATE");

  
  // == make ABCD tree ==
  double myMET;
  double minDeltaPhi;
  double minDeltaRbj;
  double DeltaPhiMPTMET;

  TTree ABCDtree("ABCDtree","ABCD tree");
  ABCDtree.Branch("weight",&weight,"weight/D");
  ABCDtree.Branch("MET",&myMET,"MET/D");
  ABCDtree.Branch("minDeltaPhi",&minDeltaPhi,"minDeltaPhi/D");
  ABCDtree.Branch("minDeltaRbj",&minDeltaRbj,"minDeltaRbj/D");
  ABCDtree.Branch("DeltaPhiMPTMET",&DeltaPhiMPTMET,"DeltaPhiMPTMET/D");
  //some additional variables, not really for ABCD....
  int ntightjets;
  int ntightMCbjets;
  int nloosejets;
  int nlooseMCbjets;
  ABCDtree.Branch("ntightjets",&ntightjets,"ntightjets/I");
  ABCDtree.Branch("ntightMCbjets",&ntightMCbjets,"ntightMCbjets/I");
  ABCDtree.Branch("nloosejets",&nloosejets,"nloosejets/I");
  ABCDtree.Branch("nlooseMCbjets",&nlooseMCbjets,"nlooseMCbjets/I");

  //
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    if (Cut(ientry) < 0) continue; //jmt use cut
    myMET = MET; //could and should modify to use tcMET or caloMET
    minDeltaPhi = getMinDeltaPhiMET(3);
    minDeltaRbj = getOverallMinDeltaR_bj();
    DeltaPhiMPTMET = getDeltaPhiMPTMET();

    ntightjets=jetPt->size();
    nloosejets=loosejetPt->size();
    ntightMCbjets=0;
    nlooseMCbjets=0;
    for (unsigned int ijet=0; ijet<jetFlavor->size(); ijet++)
      if (abs(jetFlavor->at(ijet))==5) ntightMCbjets++;
    for (unsigned int ijet=0; ijet<loosejetFlavor->size(); ijet++)
      if (abs(loosejetFlavor->at(ijet))==5) nlooseMCbjets++;
    
    ABCDtree.Fill(); 
  }

  fout.Write();
  fout.Close();
  
}

/* 
main plot-making loop
*/
void basicLoop::Loop()
{
  const double pi=4*atan(1);
  //
 
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries(); //jmt: remove Fast


   TString inname=findInputName(); //uses fChain
   std::cout<<"Got an input file name as: "<<inname<<std::endl;
   
   double sigma = getCrossSection(inname);
   TString sampleName = getSampleName(inname);
   if (sigma<=0) return;

   if (nentries == 0) {std::cout<<"Chain has no entries!"<<std::endl; return;}

   double   weight = lumi * sigma / double(nentries); //calculate weight

   //open output file
   TString outfilename="plots."; 
   outfilename+=getCutDescriptionString();
   outfilename+=".";    outfilename+=sampleName; 
   outfilename+=".root";
   TFile fout(outfilename,"RECREATE");

  
   // === make some histograms ===
   const   int offset=3;
   int njets_bins=10-offset;
   TH1D Hnjets_nocuts("Hnjets_nocuts","N of jets (no cuts)",njets_bins,offset,njets_bins+offset);
   TH1D Hnjets("Hnjets","N of jets (RA2)",njets_bins,offset,njets_bins+offset);
   TH1D Hnjets_ge1b("Hnjets_ge1b","N of jets (RA2 && >=1b)",njets_bins,offset,njets_bins+offset);
   TH1D Hnjets_ge2b("Hnjets_ge2b","N of jets (RA2 && >=2b)",njets_bins,offset,njets_bins+offset);
   TH1D Hnjets_ge3b("Hnjets_ge3b","N of jets (RA2 && >=3b)",njets_bins,offset,njets_bins+offset);

   //DeltaPhi(MPT,MET)
   int nbins=50;

   TH1D HdeltaPhiMPTMET("HdeltaPhiMPTMET","DeltaPhi(MET,MPT) (RA2)",nbins,0,pi);
   TH1D HdeltaPhiMPTMET_ge2b("HdeltaPhiMPTMET_ge2b","DeltaPhi(MET,MPT) (RA2 && >=2b)",nbins,0,pi);

   TH1D HminDeltaPhiMETb_ge1b("HminDeltaPhiMETb_ge1b","minDeltaPhi(b,MET) (RA2 && >=1b)",nbins,0,pi);
   TH1D HminDeltaPhiMETj_ge1b("HminDeltaPhiMETj_ge1b","minDeltaPhi(j,MET) (RA2 && >=1b)",nbins,0,pi);
   TH2D HdeltaPhib1b2_minDeltaPhiMETb("HdeltaPhib1b2_minDeltaPhiMETb","DeltaPhi(b1,b2) v minDeltaPhi(b,MET) (RA2 && >=2b)",nbins,0,pi,nbins,0,pi);

   TH1D HminDeltaPhiMETb_ge2b("HminDeltaPhiMETb_ge2b","minDeltaPhi(b,MET) (RA2 && >=2b)",nbins,0,pi);
   TH1D HminDeltaPhiMETj_ge2b("HminDeltaPhiMETj_ge2b","minDeltaPhi(j,MET) (RA2 && >=2b)",nbins,0,pi);

   double met_max=500,met_min=0;
   TH1D H_MHT("H_MHT","MHT (RA2)",nbins, met_min , met_max);
   TH1D H_MET("H_MET","MET (RA2)",nbins, met_min , met_max);

   TH1D H_MHT_ge1b("H_MHT_ge1b","MHT (RA2 && >=1b)",nbins, met_min , met_max);
   TH1D H_MET_ge1b("H_MET_ge1b","MET (RA2 && >=1b)",nbins, met_min , met_max);

   TH1D H_MHT_ge2b("H_MHT_ge2b","MHT (RA2 && >=2b)",nbins, met_min , met_max);
   TH1D H_MET_ge2b("H_MET_ge2b","MET (RA2 && >=2b)",nbins, met_min , met_max);

   TH1D H_MHT_ge3b("H_MHT_ge3b","MHT (RA2 && >=3b)",nbins, met_min , met_max);
   TH1D H_MET_ge3b("H_MET_ge3b","MET (RA2 && >=3b)",nbins, met_min , met_max);

   TH2D HpassMET45("HpassMET45","pass MET45 versus MET (RA2)",nbins,met_min,met_max,2,0,2);
   TH2D HpassMET45_ge1b("HpassMET45_ge1b","pass MET45 versus MET (RA2 && >=1b)",nbins,met_min,met_max,2,0,2);
   TH2D HpassMET45_ge2b("HpassMET45_ge2b","pass MET45 versus MET (RA2 && >=2b)",nbins,met_min,met_max,2,0,2);

   TH2D HdeltaPhib1b2_MET("HdeltaPhib1b2_MET","DeltaPhi(b1,b2) v MET (>=2b)",nbins,met_min,met_max,nbins,0,pi);
   TH2D HminDeltaPhiMETb_MET_ge2b("HminDeltaPhiMETb_MET_ge2b","minDeltaPhi(b,MET) v MET (>=2b)",nbins,met_min,met_max,nbins,0,pi);
   TH2D HminDeltaPhiMETj_MET_ge2b("HminDeltaPhiMETj_MET_ge2b","minDeltaPhi(j,MET) v MET (>=2b)",nbins,met_min,met_max,nbins,0,pi);

   TH2D HminDeltaPhiMETb_HminDeltaPhiMETj_ge2b("HminDeltaPhiMETb_HminDeltaPhiMETj_ge2b","minDPhi(b,MET) v minDPhi(j,MET) (>=2b)",
					       nbins,0,pi,nbins,0,pi);
   TH2D HdeltaPhiMPTMET_MET_ge2b("HdeltaPhiMPTMET_MET_ge2b","DeltaPhi(MET,MPT) v MET (RA2 && >=2b)",nbins,met_min,met_max,nbins,0,pi);

   TH1D HtopDecayCategory_ge2b("HtopDecayCategory_ge2b","top decay category (RA2 && >=2b)",nTopCategories,-0.5,nTopCategories-0.5);

   TH1D HdeltaPhi_bj_ge2b("HdeltaPhi_bj_ge2b","deltaPhi(b,j) (RA2 && >=2b)",nbins,0,pi);

   float dr_min=0.5;
   float dr_max=5.5;

   TH1D HdeltaR_bj_ge2b("HdeltaR_bj_ge2b","deltaR(b,j) (RA2 && >=2b)",nbins,dr_min,dr_max);
   TH2D HdeltaR_bj_vTopCat_ge2b("HdeltaR_bj_vTopCat_ge2b","deltaR(b,j) (RA2 && >=2b)",nTopCategories,-0.5,nTopCategories-0.5,nbins,dr_min,dr_max);

   TH1D HminDeltaR_bj_ge2b("HminDeltaR_bj_ge2b","minimum deltaR(b,j) (RA2 && >=2b)",nbins,0,5);
   TH2D HminDeltaR_bj_vTopCat_ge2b("HminDeltaR_bj_vTopCat_ge2b","minimum deltaR(b,j) (RA2 && >=2b)",nTopCategories,-0.5,nTopCategories-0.5,nbins,dr_min,dr_max);

   TH2D HminDeltaR_bj_vMET_ge1b("HminDeltaR_bj_vMET_ge1b","minimum deltaR(b,j) v. MET (RA2 && >=1b)",nbins,met_min,met_max,nbins,dr_min,dr_max);
   TH2D HminDeltaR_bj_vMET_ge2b("HminDeltaR_bj_vMET_ge2b","minimum deltaR(b,j) v. MET (RA2 && >=2b)",nbins,met_min,met_max,nbins,dr_min,dr_max);

   int nbinsDR=20;
   nbins=40;
   met_min=40;
   met_max=met_min+nbins*10;
   TH2D HminDeltaR_bj_vMET_ABCD_ge1b("HminDeltaR_bj_vMET_ABCD_ge1b","minimum deltaR(b,j) v. MET (RA2 && >=1b)",nbins,met_min,met_max,nbinsDR,dr_min,dr_max);
   TH2D HminDeltaR_bj_vMET_ABCD_ge2b("HminDeltaR_bj_vMET_ABCD_ge2b","minimum deltaR(b,j) v. MET (RA2 && >=2b)",nbins,met_min,met_max,nbinsDR,dr_min,dr_max);

   double ht_min=0,ht_max=900;
   int ht_nbins=300;
   TH2D HpassHT100U("HpassHT100U","pass HT100U versus HT (RA2)",ht_nbins,ht_min,ht_max,2,0,2);
   TH2D HpassHT100U_ge1b("HpassHT100U_ge1b","pass HT100U versus HT (RA2 && >=1b)",ht_nbins,ht_min,ht_max,2,0,2);
   TH2D HpassHT100U_ge2b("HpassHT100U_ge2b","pass HT100U versus HT (RA2 && >=2b)",ht_nbins,ht_min,ht_max,2,0,2);

   //as usual, perhaps we should manage our histos with HistHolder (but let's keep it simple instead)
   Hnjets.Sumw2();
   Hnjets_nocuts.Sumw2();
   Hnjets_ge1b.Sumw2();
   Hnjets_ge2b.Sumw2();
   Hnjets_ge3b.Sumw2();

   HdeltaPhiMPTMET.Sumw2();
   HdeltaPhiMPTMET_ge2b.Sumw2();
   HdeltaPhib1b2_minDeltaPhiMETb.Sumw2();
   HminDeltaPhiMETb_ge1b.Sumw2();
   HminDeltaPhiMETj_ge1b.Sumw2();

   HminDeltaPhiMETb_ge2b.Sumw2();
   HminDeltaPhiMETj_ge2b.Sumw2();

   H_MHT.Sumw2();
   H_MET.Sumw2();
   H_MHT_ge1b.Sumw2();
   H_MET_ge1b.Sumw2();
   H_MHT_ge2b.Sumw2();
   H_MET_ge2b.Sumw2();
   H_MHT_ge3b.Sumw2();
   H_MET_ge3b.Sumw2();

   HdeltaPhiMPTMET_MET_ge2b.Sumw2();

   HdeltaPhib1b2_MET.Sumw2();
   HminDeltaPhiMETb_MET_ge2b.Sumw2();
   HminDeltaPhiMETj_MET_ge2b.Sumw2();

   HminDeltaPhiMETb_HminDeltaPhiMETj_ge2b.Sumw2();
   HtopDecayCategory_ge2b.Sumw2();
   HdeltaPhi_bj_ge2b.Sumw2();
   HdeltaR_bj_ge2b.Sumw2();
   HdeltaR_bj_vTopCat_ge2b.Sumw2();

   HminDeltaR_bj_ge2b.Sumw2();
   HminDeltaR_bj_vTopCat_ge2b.Sumw2();

   HminDeltaR_bj_vMET_ge1b.Sumw2();
   HminDeltaR_bj_vMET_ge2b.Sumw2();

   HminDeltaR_bj_vMET_ABCD_ge1b.Sumw2();
   HminDeltaR_bj_vMET_ABCD_ge2b.Sumw2();

   HpassHT100U.Sumw2();
   HpassHT100U_ge1b.Sumw2();
   HpassHT100U_ge2b.Sumw2();

   HpassMET45.Sumw2();
   HpassMET45_ge1b.Sumw2();
   HpassMET45_ge2b.Sumw2();

   //keep track of performance
   TDatime starttime; //default ctor is for current time

   //event loop
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries ;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) {assert(0);}
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      Hnjets_nocuts.Fill( jetPt->size(), weight );
      if (Cut(ientry) < 0) continue; //jmt use cut

      //calculate things
      double dp_MPTMET = getDeltaPhiMPTMET();

      Hnjets.Fill( jetPt->size(), weight );
      HdeltaPhiMPTMET.Fill( dp_MPTMET,weight);
      H_MHT.Fill(MHT ,weight);
      H_MET.Fill(MET ,weight);

      int passHTtrig = passTrigger->at(0) ? 1:0; //index 0 is HLT_HT100U
      int passMETtrig = passTrigger->at(2) ? 1:0; //index 2 is HLT_MET45
      HpassHT100U.Fill(HT,passHTtrig); //note -- NOT using weight here!
      HpassMET45.Fill(MET,passMETtrig);

      if ( nbSSVM < 1) continue; //cut on the number of b tags
      Hnjets_ge1b.Fill( jetPt->size(), weight );
      H_MHT_ge1b.Fill(MHT ,weight);
      H_MET_ge1b.Fill(MET ,weight);

      HpassHT100U_ge1b.Fill(HT,passHTtrig); //note -- NOT using weight here!
      HpassMET45_ge1b.Fill(MET,passMETtrig);

      double minDeltaPhi_b_MET= getMinDeltaPhibMET();
      double minDeltaPhi_j_MET= getMinDeltaPhiMET(3);
      HminDeltaPhiMETb_ge1b.Fill(minDeltaPhi_b_MET,weight);
      HminDeltaPhiMETj_ge1b.Fill(minDeltaPhi_j_MET,weight);

      int topcat = getTopDecayCategory();
      double minDeltaR_bj=999;
      //note that all tight jet vectors should have the same size
      for (unsigned int ib = 0; ib< jetPhi->size(); ib++) {
	if ( passBCut(ib)) { //refind the b jets
	  double deltaPhi_bj=getMinDeltaPhi_bj(ib);
	  if ( nbSSVM < 2) {	  HdeltaPhi_bj_ge2b.Fill(deltaPhi_bj,weight);}

	  double mdr=getMinDeltaR_bj(ib);
	  if ( nbSSVM < 2) {
	    HdeltaR_bj_ge2b.Fill(mdr,weight);
	    HdeltaR_bj_vTopCat_ge2b.Fill(topcat,mdr,weight);
	  }
	  if (mdr<minDeltaR_bj) minDeltaR_bj=mdr;
	}
      }

      HminDeltaR_bj_vMET_ge1b.Fill(MET, minDeltaR_bj,weight);
      HminDeltaR_bj_vMET_ABCD_ge1b.Fill(MET, minDeltaR_bj,weight);

      if ( nbSSVM < 2) continue; //cut on the number of b tags

      Hnjets_ge2b.Fill( jetPt->size(), weight );
      H_MHT_ge2b.Fill(MHT ,weight);
      H_MET_ge2b.Fill(MET ,weight);
      HdeltaPhiMPTMET_ge2b.Fill( dp_MPTMET,weight);

      HpassHT100U_ge2b.Fill(HT,passHTtrig);
      HpassMET45_ge2b.Fill(MET,passMETtrig);

      HminDeltaPhiMETb_ge2b.Fill(minDeltaPhi_b_MET,weight);
      HminDeltaPhiMETj_ge2b.Fill(minDeltaPhi_j_MET,weight);

      double deltaPhi_b1b2 = getDeltaPhib1b2();
      HdeltaPhib1b2_minDeltaPhiMETb.Fill(minDeltaPhi_b_MET,deltaPhi_b1b2,weight);

      HdeltaPhib1b2_MET.Fill(MET,deltaPhi_b1b2,weight);
      HminDeltaPhiMETb_MET_ge2b.Fill(MET,minDeltaPhi_b_MET,weight);
      HminDeltaPhiMETj_MET_ge2b.Fill(MET,minDeltaPhi_j_MET,weight);

      HminDeltaPhiMETb_HminDeltaPhiMETj_ge2b.Fill(minDeltaPhi_j_MET,minDeltaPhi_b_MET,weight);

      HdeltaPhiMPTMET_MET_ge2b.Fill(MET,dp_MPTMET,weight);

      HtopDecayCategory_ge2b.Fill(topcat,weight);

      //make a plot that has only one entry per event
      HminDeltaR_bj_ge2b.Fill(minDeltaR_bj,weight);
      HminDeltaR_bj_vTopCat_ge2b.Fill(topcat,minDeltaR_bj,weight);

      HminDeltaR_bj_vMET_ge2b.Fill(MET, minDeltaR_bj,weight);
      HminDeltaR_bj_vMET_ABCD_ge2b.Fill(MET, minDeltaR_bj,weight);

      if ( nbSSVM < 3) continue; //cut on the number of b tags
      Hnjets_ge3b.Fill( jetPt->size(), weight );
      H_MHT_ge3b.Fill(MHT ,weight);
      H_MET_ge3b.Fill(MET ,weight);
   }
   TDatime stoptime; //default ctor is for current time
   UInt_t elapsed= stoptime.Convert() - starttime.Convert();

   cout<<"events / time = "<<nentries<<" / "<<elapsed<<" = "<<double(nentries)/double(elapsed)<<" Hz"<<endl;
   fout.Write();
   fout.Close();

}


void basicLoop::nbLoop()
{
//   In a ROOT session, you can do:
//      Root > .L basicLoop.C
//      Root > basicLoop t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   TFile fout("SUSYbjets.MMSSM.root","RECREATE");
   
   // define some histograms here (one for before cuts, one for after)
   TH1D HSUSY_nb("HSUSY_nb","#b from SUSY decay (no cuts)",10,0,10);
   TH1D HSUSY_nb_RA2("HSUSY_nb_RA2","#b from SUSY decay (RA2 cuts)",10,0,10);

   Long64_t nentries = fChain->GetEntries(); //jmt: remove Fast

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      //fill before cut histograms
      HSUSY_nb.Fill(SUSY_nb);

      if (Cut(ientry) < 0) continue; //jmt use cut
      //fill after cut histograms
      HSUSY_nb_RA2.Fill(SUSY_nb);
   }

  fout.Write();
  fout.Close();

}


void basicLoop::compareRA2()
{
  Long64_t nbad=0,ngood=0;

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries(); //jmt: remove Fast

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      bool passJMT = Cut(ientry) > 0;
      bool passDon = SUSY_cutResults->at(12);

      if ( passJMT == passDon) {ngood++;}
      else {nbad++;}

   }

   std::cout<<"Good = "<<ngood<<std::endl;
   std::cout<<"Bad  = "<<nbad<<std::endl;

}
