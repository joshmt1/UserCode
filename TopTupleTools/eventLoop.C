#define eventLoop_cxx
#include "eventLoop.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TMath.h"
#include "TMatrixT.h"
#include "TVectorT.h"

#include <fstream>

using namespace std;

//btag cuts
const Float_t kTCHEL=1.7,kTCHEM=3.3,kTCHET=10.2,kSSVL=1.05,kSSVM=1.74,kSSVT=3.05,
  kJBPL=0.988,kJBPM=1.83,kJBPT=1.95,kCSVL=0.38,kCSVM=0.75,kCSVT=0.921;

void eventLoop::Loop(TString outfile,TString weightfile, bool isSignal)
{

//   In a ROOT session, you can do:
//      Root > .L eventLoop.C
//      Root > eventLoop t
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

   TFile fin(weightfile);
   TH1D* Hweight=0;
   if (weightfile != "") {
     if (!fin.IsZombie()) {
       std::cout<<"Loading weights from: "<<weightfile<<std::endl;
       Hweight = (TH1D*) fin.Get("Hweight");
     }
     else {
       std::cout<<"ERROR loading weights from: "<<weightfile<<std::endl;
       return;
     }
   }
   
   TString asciifile=outfile;
   ofstream* fascii=0;
   if (writeASCII_) {
     asciifile.Replace( asciifile.Index(".root"), 5, ".dat");
     cout<<"[eventLoop] going to write to ASCII file: "<<asciifile<<endl;
     fascii = new ofstream(asciifile.Data());

     vector<TString> asciinames;
     //careful! this list needs to match the variables at the end! and in the same order!
     asciinames.push_back("HT");
     asciinames.push_back("MET");

     asciinames.push_back("jet_pT1");
     asciinames.push_back("jet_pT2");
     asciinames.push_back("jet_pT3");

     asciinames.push_back("ST");
     asciinames.push_back("minDeltaPhi");
     asciinames.push_back("maxDeltaPhi");

     asciinames.push_back("bjet_eta1");
     asciinames.push_back("DeltaPhi_MPTLoose_MET");

     *fascii<<asciinames.size()<<endl;
     for (vector<TString>::const_iterator iname=asciinames.begin(); iname!=asciinames.end(); ++iname) {
       *fascii<<*iname<<" ";
     }
     *fascii<<endl;
   }

   cout<<"[eventLoop] about to open file: "<<outfile<<endl;
   TFile fout(outfile,"RECREATE");

   Int_t ngoodjets;
   Float_t eventweight=1;
   Float_t jet_pT1,jet_pT2,jet_pT3;   
   Float_t jet_eta1,jet_eta2,jet_eta3;
   Float_t bjet_eta1,bjet_eta2;
   Float_t bjet_pT1,bjet_pT2,bjet_pT3;
   Float_t jet_EMf1,jet_EMf2,jet_EMf3;
   Float_t ST,HT,Meff;
   Float_t minDeltaPhi,maxDeltaPhi,DeltaPhi1,DeltaPhi2,DeltaPhi3;
   Float_t DeltaPhi_b1b2,DeltaPhi_j1b1;
   Float_t MPT,DeltaPhi_MPTLoose_MET;

   Int_t nbTCHEL=0,nbTCHEM=0,nbTCHET=0,nbSSVL=0,nbSSVM=0,nbSSVT=0,
    nbJBPL=0,nbJBPM=0,nbJBPT=0,nbCSVL=0,nbCSVM=0,nbCSVT=0;

   Int_t MC_Nb=0,MC_Ng=0,MC_Nuds=0,MC_Nc=0,MC_Nunmatched=0;

   Int_t type;

   TTree tree("tree","tree");
   tree.Branch("eventweight",&eventweight,"eventweight");
   tree.Branch("MET",&metPt,"MET");
   tree.Branch("MC_MET",&genMetPt,"MC_MET");

   tree.Branch("HT",&HT,"HT");
   tree.Branch("Meff",&Meff,"Meff");

   tree.Branch("ngoodjets",&ngoodjets,"ngoodjets/I");
   tree.Branch("jet_pT1",&jet_pT1,"jet_pT1");
   tree.Branch("jet_pT2",&jet_pT2,"jet_pT2");
   tree.Branch("jet_pT3",&jet_pT3,"jet_pT3");

   tree.Branch("bjet_pT1",&bjet_pT1,"bjet_pT1");
   tree.Branch("bjet_pT2",&bjet_pT2,"bjet_pT2");
   tree.Branch("bjet_pT3",&bjet_pT3,"bjet_pT3");

   tree.Branch("jet_eta1",&jet_eta1,"jet_eta1");
   tree.Branch("jet_eta2",&jet_eta2,"jet_eta2");
   tree.Branch("jet_eta3",&jet_eta3,"jet_eta3");

   tree.Branch("bjet_eta1",&bjet_eta1,"bjet_eta1");
   tree.Branch("bjet_eta2",&bjet_eta2,"bjet_eta2");

   tree.Branch("jet_EMf1",&jet_EMf1,"jet_EMf1");
   tree.Branch("jet_EMf2",&jet_EMf2,"jet_EMf2");
   tree.Branch("jet_EMf3",&jet_EMf3,"jet_EMf3");

   tree.Branch("minDeltaPhi",&minDeltaPhi,"minDeltaPhi");
   tree.Branch("maxDeltaPhi",&maxDeltaPhi,"maxDeltaPhi");
   tree.Branch("DeltaPhi1",&DeltaPhi1,"DeltaPhi1");
   tree.Branch("DeltaPhi2",&DeltaPhi2,"DeltaPhi2");
   tree.Branch("DeltaPhi3",&DeltaPhi3,"DeltaPhi3");

   //angles relating to b jets
   // -- * angle between two lead b jets
   tree.Branch("DeltaPhi_b1b2",&DeltaPhi_b1b2,"DeltaPhi_b1b2");
   // -- * angle between lead b jet and lead jet
   tree.Branch("DeltaPhi_j1b1",&DeltaPhi_j1b1,"DeltaPhi_j1b1");

   tree.Branch("MPT",&MPT,"MPT");
   tree.Branch("DeltaPhi_MPTLoose_MET",&DeltaPhi_MPTLoose_MET,"DeltaPhi_MPTLoose_MET");

   tree.Branch("nbTCHEL",&nbTCHEL,"nbTCHEL/I");
   tree.Branch("nbTCHEM",&nbTCHEM,"nbTCHEM/I");
   tree.Branch("nbTCHET",&nbTCHET,"nbTCHET/I");
   tree.Branch("nbSSVL",&nbSSVL,"nbSSVL/I");
   tree.Branch("nbSSVM",&nbSSVM,"nbSSVM/I");
   tree.Branch("nbSSVT",&nbSSVT,"nbSSVT/I");
   tree.Branch("nbJBPL",&nbJBPL,"nbJBPL/I");
   tree.Branch("nbJBPM",&nbJBPM,"nbJBPM/I");
   tree.Branch("nbJBPT",&nbJBPT,"nbJBPT/I");
   tree.Branch("nbCSVL",&nbCSVL,"nbCSVL/I");
   tree.Branch("nbCSVM",&nbCSVM,"nbCSVM/I");
   tree.Branch("nbCSVT",&nbCSVT,"nbCSVT/I");

   tree.Branch("MC_Ng",&MC_Ng,"MC_Ng/I");
   tree.Branch("MC_Nuds",&MC_Nuds,"MC_Nuds/I");
   tree.Branch("MC_Nc",&MC_Nc,"MC_Nc/I");
   tree.Branch("MC_Nb",&MC_Nb,"MC_Nb/I");
   tree.Branch("MC_Nunmatched",&MC_Nunmatched,"MC_Nunmatched/I");


   //event shape
   tree.Branch("ST",&ST,"ST");
   //these are all 'automatically filled' (just a copy from source ntuple)
   tree.Branch("sphericity",&topo_sphericity,"sphericity");
   tree.Branch("aplanarity",&topo_aplanarity,"aplanarity");
   tree.Branch("sphericity_e",&topo_sphericity_e,"sphericity_e");
   tree.Branch("aplanarity_e",&topo_aplanarity_e,"aplanarity_e");
   tree.Branch("ht_precomp",&topo_ht,"ht_precomp");
   //there are certainly other interesting variables in the source ntuple
   //e.g. met significance
   //can look at them later

   tree.Branch("type",&type,"type/I"); //special signal versus bkd qualifier
   type = (isSignal) ? 1 : 0;

   Long64_t nentries = fChain->GetEntries(); //i don't trust Fast
   cout<<"Chain has n="<<nentries<<endl;
   Long64_t nbytes = 0, nb = 0;
   //Loop over events
   for (Long64_t jentry=0; jentry<nentries; jentry++) {
     
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     
     if (Cut() < 0) continue;

     if (Hweight != 0) {
       int bin=Hweight->FindBin(processPtHat);
       eventweight= Hweight->GetBinContent(bin);
     }
     else {eventweight=1;}

     // ==== fill simple tree variables ====
     //MET is filled 'automatically
     HT=0;
     Meff = metPt; //we will add jet et later

     ngoodjets=nGoodJets();
     //note that these variables are for the 3 lead jets,
     //whether they are in the 'good' eta range or not!
     jet_pT1 = jetPt[0];
     jet_pT2 = jetPt[1];
     jet_pT3 = jetPt[2];
     jet_EMf1 = jetEMEnergyFraction[0];
     jet_EMf2 = jetEMEnergyFraction[1];
     jet_EMf3 = jetEMEnergyFraction[2];
     jet_eta1 = jetEta[0];
     jet_eta2 = jetEta[1];
     jet_eta3 = jetEta[2];

     DeltaPhi1 = getDeltaPhi(jetPhi[0],metPhi);
     DeltaPhi2 = getDeltaPhi(jetPhi[1],metPhi);
     DeltaPhi3 = getDeltaPhi(jetPhi[2],metPhi);

     //loop over jets
     minDeltaPhi=4,maxDeltaPhi=-1;
     MC_Nb=0; MC_Ng=0; MC_Nuds=0; MC_Nc=0; MC_Nunmatched=0;
     nbTCHEL=0;nbTCHEM=0;nbTCHET=0;nbSSVL=0;nbSSVM=0;nbSSVT=0;
     nbJBPL=0;nbJBPM=0;nbJBPT=0;nbCSVL=0;nbCSVM=0;nbCSVT=0;
     Float_t S11=0,S12=0, S22=0;
     vector<int> goodbjets;
     for (int j=0 ; j<numJet ; j++) {

       //when counting true jet flavor, i don't care if the jet is 'good' or not
       if ( fabs(jetPID[j]) == 21 ) { MC_Ng++; }
       else if ( fabs(jetPID[j]) == 5 ) {MC_Nb++;}
       else if ( fabs(jetPID[j]) == 4 ) {MC_Nc++;}
       else if ( fabs(jetPID[j]) == 3 || fabs(jetPID[j]) == 2 || fabs(jetPID[j]) == 1 ) {MC_Nuds++;}
       else if ( fabs(jetPID[j]) == 0 ) {MC_Nunmatched++;}
       else {
	 MC_Nunmatched++;
	 cout<<"jet with jetPID = "<<jetPID[j]<<endl;
       }

       if (isJetGood(j)) {
	 HT+=jetPt[j];
	 Meff+=jetEt[j];

	 float dp = getDeltaPhi(jetPhi[j],metPhi);
	 if (dp < minDeltaPhi) minDeltaPhi = dp;
	 if (dp > maxDeltaPhi) maxDeltaPhi = dp;

	 S11 += jetPx[j] * jetPx[j];
	 S12 += jetPx[j] * jetPy[j];
	 S22 += jetPy[j] * jetPy[j];

	 //count b jets of various kinds
	 bool isB=false;
	 if (jetBtagTrackCountHighEff[j] > kTCHEL) nbTCHEL++;
	 if (jetBtagTrackCountHighEff[j] > kTCHEM) nbTCHEM++;
	 if (jetBtagTrackCountHighEff[j] > kTCHET) nbTCHET++;

	 if (jetBtagSecondaryVertex[j] > kSSVL) nbSSVL++;
	 if (jetBtagSecondaryVertex[j] > kSSVM) {nbSSVM++; isB=true;}
	 if (jetBtagSecondaryVertex[j] > kSSVT) nbSSVT++;
	 
	 if (jetBtagProbability[j] > kJBPL) nbJBPL++;
	 if (jetBtagProbability[j] > kJBPM) nbJBPM++;
	 if (jetBtagProbability[j] > kJBPT) nbJBPT++;

	 if (isB) {
	   goodbjets.push_back(j);
	 }
       }
     } //end loop over jets

     //this code allows for cutting on the number of b jets!
     //it makes exclusive samples for a cut of 1 or 2
     //for 3 it makes inclusive nb>=3
     if (goodBJets_ == 2 || goodBJets_==1) {
       if (goodbjets.size() != goodBJets_) continue;       
     }
     else if (goodBJets_ >= 3) {
       if (goodbjets.size() < goodBJets_) continue;       
     }

     TMatrixT<float> Smat(2,2);
     Smat[0][0]=S11;
     Smat[0][1]=S12;
     Smat[1][0]=S12;
     Smat[1][1]=S22;
     TVectorT<float> eigenvalues;
     TMatrixT<float> eigenvectors = Smat.EigenVectors(eigenvalues);
     ST=2*eigenvalues[1] / (eigenvalues[0]+eigenvalues[1]);

     bjet_pT1=-1; bjet_pT2=-1; bjet_pT3=-1;
     for (vector<int>::const_iterator ij=goodbjets.begin(); ij!=goodbjets.end(); ++ij) {
       if      (bjet_pT1==-1) bjet_pT1=jetPt[*ij];
       else if (bjet_pT2==-1) bjet_pT2=jetPt[*ij];
       else if (bjet_pT3==-1) bjet_pT3=jetPt[*ij];
     }

     if (goodbjets.size() >= 2) {
       DeltaPhi_b1b2 = getDeltaPhi( jetPhi[goodbjets[0]], jetPhi[goodbjets[1]]);
       bjet_eta2 = jetEta[ goodbjets[1]];
     }
     else {DeltaPhi_b1b2 = -1; bjet_eta2=-1;}
     if (goodbjets.size() >= 1) {
       //here i am finding the angle between the lead b jet that has passed
       //jet cuts, and the lead jet (not necessarily passing cuts)
       //is this what i want?
       DeltaPhi_j1b1 = getDeltaPhi( jetPhi[goodbjets[0]], jetPhi[0]);
       bjet_eta1 = jetEta[ goodbjets[0]];
     }
     else {DeltaPhi_j1b1=-1; bjet_eta1=-1;}

     //loop over tracks
     Float_t MPTx=0,MPTy=0;
     for (int i=0; i<numGeneralTracks; i++) {
       //loose cuts; should be codified into a configurable parameter
       if ( generalTracksPt[i] > 5 && fabs(generalTracksEta[i]) < 5) {
	 MPTx += generalTracksPt[i] * cos(generalTracksPhi[i]);
	 MPTy += generalTracksPt[i] * sin(generalTracksPhi[i]);
       }
     }

     MPT = sqrt(MPTx*MPTx + MPTy*MPTy);
     DeltaPhi_MPTLoose_MET = getDeltaPhi( metPhi, TMath::ATan2(MPTy,MPTx));

     if (jentry<10)  cout<<"MC_Nb = "<<MC_Nb<<endl;
     tree.Fill();

     if (writeASCII_) {
       *fascii
	 <<HT<<" "
	 <<metPt<<" " //MET
	 <<jet_pT1<<" "
	 <<jet_pT2<<" "
	 <<jet_pT3<<" "
	 <<ST<<" "
	 <<minDeltaPhi<<" "
	 <<maxDeltaPhi<<" "
	 <<bjet_eta1<<" "
	 <<DeltaPhi_MPTLoose_MET<<" "
	 <<eventweight<<" "<<type<<endl;
     }
   } //loop over events

   fout.Write();
   fout.Close();
   fin.Close();

   if (writeASCII_) fascii->close();

}
