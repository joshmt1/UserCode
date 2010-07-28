#define basicLoop_cxx
#include "basicLoop.h"
#include <TH1.h>
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

void basicLoop::cutflow()
{

  const double lumi=100;
  
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
      //dirty hack to include b tagging in the cut flow
      //	else if ( i>=13 && cutResults->at(i) )  npass.at(i) = npass.at(i) +1;
      //	else if ( i>=13 && !cutResults->at(i) )  break;
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

void basicLoop::Loop()
{
  const double pi=4*atan(1);
  //
 
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries(); //jmt: remove Fast

   //we want to weight to 100 pb^-1
   const double lumi=100;

   TString inname=findInputName(); //uses fChain
   std::cout<<"Got an input file name as: "<<inname<<std::endl;
   
   double sigma = getCrossSection(inname);
   TString sampleName = getSampleName(inname);
   if (sigma<=0) return;

   if (nentries == 0) {std::cout<<"Chain has no entries!"<<std::endl; return;}

   double   weight = lumi * sigma / double(nentries); //calculate weight

   //open output file
   TString outfilename="plots."; outfilename+=sampleName; outfilename+=".root";
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

   //as usual, perhaps we should manage our histos with HistHolder (but let's keep it simple instead)
   Hnjets.Sumw2();
   Hnjets_nocuts.Sumw2();
   Hnjets_ge1b.Sumw2();
   Hnjets_ge2b.Sumw2();
   Hnjets_ge3b.Sumw2();

   HdeltaPhiMPTMET.Sumw2();
   HdeltaPhiMPTMET_ge2b.Sumw2();

   //event loop
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) {assert(0);}
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      Hnjets_nocuts.Fill( jetPt->size(), weight );
      if (Cut(ientry) < 0) continue; //jmt use cut

      //calculate things
      double dp_MPTMET = getDeltaPhiMPTMET();

      Hnjets.Fill( jetPt->size(), weight );
      HdeltaPhiMPTMET.Fill( dp_MPTMET,weight);

      if ( nbSSVM < 1) continue; //cut on the number of b tags
      Hnjets_ge1b.Fill( jetPt->size(), weight );

      if ( nbSSVM < 2) continue; //cut on the number of b tags
      Hnjets_ge2b.Fill( jetPt->size(), weight );
      HdeltaPhiMPTMET_ge2b.Fill( dp_MPTMET,weight);

      if ( nbSSVM < 3) continue; //cut on the number of b tags
      Hnjets_ge3b.Fill( jetPt->size(), weight );
   }

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
