//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep 30 11:17:35 2010 by ROOT version 5.22/00d
// from TTree tree/tree
// found on file: /cu1/joshmt/BasicNtuples/V00-00-04b/LM13/BasicNtuple_1_1_3yK.root
//////////////////////////////////////////////////////////

#ifndef basicLoop_h
#define basicLoop_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// ========================================== begin
#include "/afs/cern.ch/user/j/joshmt/root/util/MiscUtil.cxx" //need to either make this file public or remove this hard-coding
//this code will have to be regenerated when changing the ntuple structure
//custom code is marked with these 'begin' and 'end' markers
// ---- this version is compatible with ntuple tag: V00-00-04, V00-00-04b, and V00-00-05 ----
// ---- ::Loop() and ::ABCDtree() are now compatible with data
#include <iostream>
#include <vector>

/*
improvement to make:
I already bifurcated the btagging from the "CutScheme" into its own thing. This is good.
I think what I really need is more than one 'cut scheme' enum.
So an enum for 'METtype' with options: MHT, MET, tcMET
and an enum for 'DeltaPhiType' with options: DeltaPhi, minDeltaPhi, etc

but i don't want to do this for now out of fear that i will screw something up.
(do it after the SUSY mtg talk)
*/

//avoid spaces and funny characters!
const char *CutSchemeNames_[]={"RA2"};
const char *METTypeNames_[]={"MHT", "MET",  "tcMET"};
const char *METRangeNames_[]={"med",  "high", "wide"}; //'no cut' is not given, because the cut can always be skipped!

const char *dpTypeNames_[]={"DeltaPhi", "minDP",  "MPT"};

//we want to weight to 50 pb^-1
const double lumi=50;
// ========================================== end

class basicLoop {
public :
  // ========================================== begin
  enum CutScheme {kRA2=0, nCutSchemes}; //after introducing these other enums, there is only one cut scheme!
  CutScheme theCutScheme_;
  enum METType {kMHT=0, kMET, ktcMET};
  METType theMETType_;
  enum METRange {kMedium=0, kHigh, kWide};
  METRange theMETRange_;
  enum dpType {kDeltaPhi=0, kminDP, kMPT};
  dpType theDPType_;
  unsigned int nBcut_;

  enum theCutFlow {cutInclusive=0,cutTrigger=1,cutPV=2,cut3Jets=3,cutJetPt1=4,cutJetPt2=5,cutJetPt3=6,cutHT=7,cutMET=8,cutMHT=9,
		   cutMuVeto=10,cutEleVeto=11,cutDeltaPhi=12,cut1B=13,cut2B=14,cut3B=15};
  std::vector<int> ignoredCut_; //allow more than 1 ignored cut!
  //if theCutFlow changes, be sure to change cutnames_ as well
  std::vector<TString> cutnames_;

  enum TopDecayCategory {kTTbarUnknown=0,kAllLeptons=1,kAllHadronic=2,kOneElectron=3,kOneMuon=4,kOneTauE=5,kOneTauMu=6,kOneTauHadronic=7,kAllTau=8,kTauPlusLepton=9, nTopCategories=10};
  // ========================================== end
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   //   vector<bool>    *SUSY_cutResults;
   vector<bool>    *cutResults;
   vector<bool>    *passTrigger;
   vector<int>     *looseJetIndex;
   vector<int>     *verylooseJetIndex;
   vector<float>   *jetPt;
   vector<float>   *jetEta;
   vector<float>   *jetPhi;
   vector<int>     *jetFlavor;
   vector<float>   *jetBTagDisc_trackCountingHighPurBJetTags;
   vector<float>   *jetBTagDisc_trackCountingHighEffBJetTags;
   vector<float>   *jetBTagDisc_simpleSecondaryVertexHighEffBJetTags;
   vector<float>   *jetBTagDisc_simpleSecondaryVertexHighPurBJetTags;
   vector<float>   *jetBTagDisc_simpleSecondaryVertexBJetTags;
   vector<float>   *loosejetPt;
   vector<float>   *loosejetEt;
   vector<float>   *loosejetEta;
   vector<float>   *loosejetPhi;
   vector<int>     *loosejetFlavor;
   vector<float>   *veryloosejetPtUncorr;
   vector<float>   *veryloosejetEtaUncorr;
   vector<float>   *veryloosejetPhiUncorr;
   vector<int>     *loosejetGenParticlePDGId;
   vector<float>   *loosejetInvisibleEnergy;
   vector<float>   *loosejetBTagDisc_trackCountingHighPurBJetTags;
   vector<float>   *loosejetBTagDisc_trackCountingHighEffBJetTags;
   vector<float>   *loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags;
   vector<float>   *loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags;
   vector<float>   *loosejetBTagDisc_simpleSecondaryVertexBJetTags;
   Int_t           nbSSVM;
   Float_t         HT;
   Float_t         MHT;
   Float_t         MHTphi;
   Float_t         DeltaPhi_JetMHT1;
   Float_t         DeltaPhi_JetMHT2;
   Float_t         DeltaPhi_JetMHT3;
   Float_t         MET;
   Float_t         METphi;
   Float_t         tcMET;
   Float_t         tcMETphi;
   vector<float>   *trackPt;
   vector<float>   *trackEta;
   vector<float>   *trackPhi;
   Int_t           SUSY_nb;
   Float_t         qScale;
   vector<int>     *topDecayCode;

   // List of branches
   //   TBranch        *b_SUSY_cutResults;   //!
   TBranch        *b_cutResults;   //!
   TBranch        *b_passTrigger;   //!
   TBranch        *b_looseJetIndex;   //!
   TBranch        *b_verylooseJetIndex;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetFlavor;   //!
   TBranch        *b_jetBTagDisc_trackCountingHighPurBJetTags;   //!
   TBranch        *b_jetBTagDisc_trackCountingHighEffBJetTags;   //!
   TBranch        *b_jetBTagDisc_simpleSecondaryVertexHighEffBJetTags;   //!
   TBranch        *b_jetBTagDisc_simpleSecondaryVertexHighPurBJetTags;   //!
   TBranch        *b_jetBTagDisc_simpleSecondaryVertexBJetTags;   //!
   TBranch        *b_loosejetPt;   //!
   TBranch        *b_loosejetEt;   //!
   TBranch        *b_loosejetEta;   //!
   TBranch        *b_loosejetPhi;   //!
   TBranch        *b_loosejetFlavor;   //!
   TBranch        *b_veryloosejetPtUncorr;   //!
   TBranch        *b_veryloosejetEtaUncorr;   //!
   TBranch        *b_veryloosejetPhiUncorr;   //!
   TBranch        *b_loosejetGenParticlePDGId;   //!
   TBranch        *b_loosejetInvisibleEnergy;   //!
   TBranch        *b_loosejetBTagDisc_trackCountingHighPurBJetTags;   //!
   TBranch        *b_loosejetBTagDisc_trackCountingHighEffBJetTags;   //!
   TBranch        *b_loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags;   //!
   TBranch        *b_loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags;   //!
   TBranch        *b_loosejetBTagDisc_simpleSecondaryVertexBJetTags;   //!
   TBranch        *b_nbSSVM;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_MHT;   //!
   TBranch        *b_MHTphi;   //!
   TBranch        *b_DeltaPhi_JetMHT1;   //!
   TBranch        *b_DeltaPhi_JetMHT2;   //!
   TBranch        *b_DeltaPhi_JetMHT3;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_METphi;   //!
   TBranch        *b_tcMET;   //!
   TBranch        *b_tcMETphi;   //!
   TBranch        *b_trackPt;   //!
   TBranch        *b_trackEta;   //!
   TBranch        *b_trackPhi;   //!
   TBranch        *b_SUSY_nb;   //!
   TBranch        *b_qScale;   //!
   TBranch        *b_topDecayCode;   //!

   basicLoop(TTree *tree=0);
   virtual ~basicLoop();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   // ========================================== begin
   virtual void     Loop(unsigned int dataindex=0);
   //   virtual void     compareRA2(); //deprecated function
   virtual void     exampleLoop();
   virtual void     nbLoop();
   virtual void     ABCDtree(unsigned int dataindex=0);

   bool isVersion04b();
   void printState();

   bool setCutScheme(CutScheme cutscheme);
   void setMETType(METType mettype);
   void setMETRange(METRange metrange);
   void setDPType(dpType dptype);
   void setIgnoredCut(const int cutIndex); //use to make N-1 plots! (and other things)
   void resetIgnoredCut() ;
   void setBCut(unsigned int nb);
   TString getCutDescriptionString();
   TString getBCutDescriptionString();
   
   void cutflow();
   bool cutRequired(unsigned int cutIndex) ;
   bool passCut(unsigned int cutIndex) ;
   TString getSampleName(const TString inname) ;
   double getCrossSection(const TString inname) ;
   TString findInputName() ;

   double getDeltaPhiMPTMET();
   double getMinDeltaPhibMET() ;
   double getMinDeltaPhiMET(unsigned int maxjets) ;
   double getMinDeltaPhitcMET(unsigned int maxjets) ;
   double getMinDeltaPhiMHT(unsigned int maxjets) ;
   double getDeltaPhib1b2();
   double getUncorrectedHT(const double threshold);
   int getTopDecayCategory();
   double getMinDeltaPhi_bj(unsigned int bindex);
   double getOverallMinDeltaR_bj();
   double getMinDeltaR_bj(unsigned int bindex);
   bool passBCut( unsigned int bindex);
   double getDeltaPhi(double phi1, double phi2);
   // ========================================== end
};

#endif

#ifdef basicLoop_cxx
basicLoop::basicLoop(TTree *tree)
//====================== begin
  :  theCutScheme_(kRA2),
     theMETType_(kMET),
     theMETRange_(kHigh),
     theDPType_(kminDP),
     nBcut_(0) 
//====================== end
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/cu1/joshmt/BasicNtuples/V00-00-04b/LM13/BasicNtuple_1_1_3yK.root");
      if (!f) {
         f = new TFile("/cu1/joshmt/BasicNtuples/V00-00-04b/LM13/BasicNtuple_1_1_3yK.root");
         f->cd("/cu1/joshmt/BasicNtuples/V00-00-04b/LM13/BasicNtuple_1_1_3yK.root:/BasicTreeMaker");
      }
      tree = (TTree*)gDirectory->Get("tree");

   }
   Init(tree);
   // ========================================== begin
   //this is now in the infotree
   //not sure how to get that unless the user passes it in as an argument
   //update -- still not using the infotree, now changing these so they can easily be used in file names
   cutnames_.push_back("Inclusive");
   cutnames_.push_back("Trigger");
   cutnames_.push_back("PV");
   cutnames_.push_back(">=3Jets");
   cutnames_.push_back("JetPt1");
   cutnames_.push_back("JetPt2");
   cutnames_.push_back("JetPt3");
   cutnames_.push_back("HT");
   cutnames_.push_back("MET");
   cutnames_.push_back("MHT");
   cutnames_.push_back("MuVeto");
   cutnames_.push_back("EleVeto");
   cutnames_.push_back("DeltaPhi");
   cutnames_.push_back(">=1b");
   cutnames_.push_back(">=2b");
   cutnames_.push_back(">=3b");
   // ========================================== end
   
}

basicLoop::~basicLoop()
{
  // ========================================== begin
  //No! stupid ROOT! this leads to a double deletion!
  //if (!fChain) return;
  // delete fChain->GetCurrentFile();
  // ========================================== end
}

Int_t basicLoop::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t basicLoop::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void basicLoop::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
  //   SUSY_cutResults = 0;
   cutResults = 0;
   passTrigger = 0;
   looseJetIndex = 0;
   verylooseJetIndex = 0;
   jetPt = 0;
   jetEta = 0;
   jetPhi = 0;
   jetFlavor = 0;
   jetBTagDisc_trackCountingHighPurBJetTags = 0;
   jetBTagDisc_trackCountingHighEffBJetTags = 0;
   jetBTagDisc_simpleSecondaryVertexHighEffBJetTags = 0;
   jetBTagDisc_simpleSecondaryVertexHighPurBJetTags = 0;
   jetBTagDisc_simpleSecondaryVertexBJetTags = 0;
   loosejetPt = 0;
   loosejetEt = 0;
   loosejetEta = 0;
   loosejetPhi = 0;
   loosejetFlavor = 0;
   veryloosejetPtUncorr = 0;
   veryloosejetEtaUncorr = 0;
   veryloosejetPhiUncorr = 0;
   loosejetGenParticlePDGId = 0;
   loosejetInvisibleEnergy = 0;
   loosejetBTagDisc_trackCountingHighPurBJetTags = 0;
   loosejetBTagDisc_trackCountingHighEffBJetTags = 0;
   loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags = 0;
   loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags = 0;
   loosejetBTagDisc_simpleSecondaryVertexBJetTags = 0;
   trackPt = 0;
   trackEta = 0;
   trackPhi = 0;
   topDecayCode = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   //   fChain->SetBranchAddress("SUSY_cutResults", &SUSY_cutResults, &b_SUSY_cutResults);
   fChain->SetBranchAddress("cutResults", &cutResults, &b_cutResults);
   fChain->SetBranchAddress("passTrigger", &passTrigger, &b_passTrigger);
   fChain->SetBranchAddress("looseJetIndex", &looseJetIndex, &b_looseJetIndex);
   if (isVersion04b())   fChain->SetBranchAddress("verylooseJetIndex", &verylooseJetIndex, &b_verylooseJetIndex); //=== begin end
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetFlavor", &jetFlavor, &b_jetFlavor);
   fChain->SetBranchAddress("jetBTagDisc_trackCountingHighPurBJetTags", &jetBTagDisc_trackCountingHighPurBJetTags, &b_jetBTagDisc_trackCountingHighPurBJetTags);
   fChain->SetBranchAddress("jetBTagDisc_trackCountingHighEffBJetTags", &jetBTagDisc_trackCountingHighEffBJetTags, &b_jetBTagDisc_trackCountingHighEffBJetTags);
   fChain->SetBranchAddress("jetBTagDisc_simpleSecondaryVertexHighEffBJetTags", &jetBTagDisc_simpleSecondaryVertexHighEffBJetTags, &b_jetBTagDisc_simpleSecondaryVertexHighEffBJetTags);
   fChain->SetBranchAddress("jetBTagDisc_simpleSecondaryVertexHighPurBJetTags", &jetBTagDisc_simpleSecondaryVertexHighPurBJetTags, &b_jetBTagDisc_simpleSecondaryVertexHighPurBJetTags);
   fChain->SetBranchAddress("jetBTagDisc_simpleSecondaryVertexBJetTags", &jetBTagDisc_simpleSecondaryVertexBJetTags, &b_jetBTagDisc_simpleSecondaryVertexBJetTags);
   fChain->SetBranchAddress("loosejetPt", &loosejetPt, &b_loosejetPt);
   fChain->SetBranchAddress("loosejetEt", &loosejetEt, &b_loosejetEt);
   fChain->SetBranchAddress("loosejetEta", &loosejetEta, &b_loosejetEta);
   fChain->SetBranchAddress("loosejetPhi", &loosejetPhi, &b_loosejetPhi);
   fChain->SetBranchAddress("loosejetFlavor", &loosejetFlavor, &b_loosejetFlavor);
   if (isVersion04b()) {//=== begin end
   fChain->SetBranchAddress("veryloosejetPtUncorr", &veryloosejetPtUncorr, &b_veryloosejetPtUncorr);
   fChain->SetBranchAddress("veryloosejetEtaUncorr", &veryloosejetEtaUncorr, &b_veryloosejetEtaUncorr);
   fChain->SetBranchAddress("veryloosejetPhiUncorr", &veryloosejetPhiUncorr, &b_veryloosejetPhiUncorr);
   }//=== begin end
   fChain->SetBranchAddress("loosejetGenParticlePDGId", &loosejetGenParticlePDGId, &b_loosejetGenParticlePDGId);
   fChain->SetBranchAddress("loosejetInvisibleEnergy", &loosejetInvisibleEnergy, &b_loosejetInvisibleEnergy);
   fChain->SetBranchAddress("loosejetBTagDisc_trackCountingHighPurBJetTags", &loosejetBTagDisc_trackCountingHighPurBJetTags, &b_loosejetBTagDisc_trackCountingHighPurBJetTags);
   fChain->SetBranchAddress("loosejetBTagDisc_trackCountingHighEffBJetTags", &loosejetBTagDisc_trackCountingHighEffBJetTags, &b_loosejetBTagDisc_trackCountingHighEffBJetTags);
   fChain->SetBranchAddress("loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags", &loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags, &b_loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags);
   fChain->SetBranchAddress("loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags", &loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags, &b_loosejetBTagDisc_simpleSecondaryVertexHighPurBJetTags);
   fChain->SetBranchAddress("loosejetBTagDisc_simpleSecondaryVertexBJetTags", &loosejetBTagDisc_simpleSecondaryVertexBJetTags, &b_loosejetBTagDisc_simpleSecondaryVertexBJetTags);
   fChain->SetBranchAddress("nbSSVM", &nbSSVM, &b_nbSSVM);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("MHT", &MHT, &b_MHT);
   fChain->SetBranchAddress("MHTphi", &MHTphi, &b_MHTphi);
   fChain->SetBranchAddress("DeltaPhi_JetMHT1", &DeltaPhi_JetMHT1, &b_DeltaPhi_JetMHT1);
   fChain->SetBranchAddress("DeltaPhi_JetMHT2", &DeltaPhi_JetMHT2, &b_DeltaPhi_JetMHT2);
   fChain->SetBranchAddress("DeltaPhi_JetMHT3", &DeltaPhi_JetMHT3, &b_DeltaPhi_JetMHT3);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("METphi", &METphi, &b_METphi);
   fChain->SetBranchAddress("tcMET", &tcMET, &b_tcMET);
   fChain->SetBranchAddress("tcMETphi", &tcMETphi, &b_tcMETphi);
   fChain->SetBranchAddress("trackPt", &trackPt, &b_trackPt);
   fChain->SetBranchAddress("trackEta", &trackEta, &b_trackEta);
   fChain->SetBranchAddress("trackPhi", &trackPhi, &b_trackPhi);
   fChain->SetBranchAddress("SUSY_nb", &SUSY_nb, &b_SUSY_nb);
   fChain->SetBranchAddress("qScale", &qScale, &b_qScale);
   fChain->SetBranchAddress("topDecayCode", &topDecayCode, &b_topDecayCode);
   Notify();
}

Bool_t basicLoop::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void basicLoop::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

// ========================================== begin

void basicLoop::setBCut(unsigned int nb) {
  nBcut_=nb;
}

void basicLoop::setIgnoredCut(const int cutIndex) {

  //could do this via exception handling....but no, i won't
  if ( cutIndex>= int(cutnames_.size())) {cout<<"Invalid cutIndex"<<endl; return;}

  ignoredCut_.push_back(cutIndex);

}

void basicLoop::resetIgnoredCut() {
  ignoredCut_.clear();
}

bool basicLoop::cutRequired(const unsigned int cutIndex) {
  //this is now implemented as an enum!
  // *        0 *        0 * Inclusive *
  // *        0 *        1 *   Trigger *
  // *        0 *        2 *        PV *
  // *        0 *        3 *  >=3 Jets *
  // *        0 *        4 * First Jet *
  // *        0 *        5 * Second Je *
  // *        0 *        6 * Third Jet *
  // *        0 *        7 *    HT Cut *
  // *        0 *        8 *   MET Cut *
  // *        0 *        9 *   MHT Cut *
  // *        0 *       10 * Muon Veto *
  // *        0 *       11 * Electron  *
  // *        0 *       12 * DeltaPhi  *
  // *        0 *       13 * >=1 B-Tag *
  // *        0 *       14 * >=2 B-Tag *
  // *        0 *       15 * >=3 B-Tag *


  //check if we are *ignoring* this cut (for N-1 studies, etc)
  for (unsigned int i = 0; i< ignoredCut_.size() ; i++) {
    if ( int(cutIndex) == ignoredCut_.at(i) ) return false;
  }

  bool cutIsRequired=false;

  //RA2
  if (theCutScheme_ == kRA2 ) {
    if      (cutIndex == cutInclusive)  cutIsRequired =  true;
    else if (cutIndex == cutTrigger)  cutIsRequired =  true;
    else if (cutIndex == cutPV)  cutIsRequired =  true;
    else if (cutIndex == cut3Jets)  cutIsRequired =  true;
    else if (cutIndex == cutJetPt1)  cutIsRequired =  false;
    else if (cutIndex == cutJetPt2)  cutIsRequired =  false;
    else if (cutIndex == cutJetPt3)  cutIsRequired =  false;
    else if (cutIndex == cutHT)  cutIsRequired =  true;
    //MET
    else if (cutIndex == cutMET)  cutIsRequired =  (theMETType_== kMET || theMETType_==ktcMET);
    //MHT
    else if (cutIndex == cutMHT)  cutIsRequired =  (theMETType_==kMHT);
    else if (cutIndex == cutMuVeto) cutIsRequired =  true;
    else if (cutIndex == cutEleVeto) cutIsRequired =  true;
    else if (cutIndex == cutDeltaPhi) cutIsRequired =  true;
    else if (cutIndex == cut1B) cutIsRequired =  nBcut_ >=1;
    else if (cutIndex == cut2B) cutIsRequired =  nBcut_ >=2;
    else if (cutIndex == cut3B) cutIsRequired =  nBcut_ >=3;
    else assert(0);
  }
  else assert(0);

  return cutIsRequired;
}

bool basicLoop::passCut(const unsigned int cutIndex) {
  //implement special exceptions here

  //MET
  //it is now an anachronism that we treat MET and MHT as different cut flow steps
  //this block of code should now evaluate *the same* for both steps
  if (cutIndex == cutMET || cutIndex==cutMHT) {
    float mymet = 0;
    if      ( theMETType_ == kMET )   mymet=MET;
    else if ( theMETType_ == ktcMET ) mymet=tcMET;
    else if ( theMETType_ == kMHT)    mymet=MHT;
    else {assert(0);}

    if (theMETRange_ == kMedium) return (mymet>=50 && mymet<150);
    //on the next line we are redoing the default that is stored in the ntuple. I think this is ok
    else if (theMETRange_ == kHigh) return (mymet >=150);
    else if (theMETRange_ == kWide) return (mymet >=50);
    else {assert(0);}
  }
  
  // -- in case we are using MET instead of MHT, we need to alter the DeltaPhi cuts --
  if (cutIndex == cutDeltaPhi && 
      (theDPType_ == kDeltaPhi && theMETType_!=kMHT)) {
    //      ( theCutScheme_ == kRA2MET ||theCutScheme_ == kRA2tcMET ) ) {
    
    float phi_of_MET = (theMETType_ == ktcMET) ? tcMETphi : METphi;
    
    int nloosejets = loosejetPhi->size();
    //need to calculate DeltaPhi between jets and MET
    //for RA2 and MHT, this is done with *loose* jets!
    double dp0=0,dp1=0,dp2=0;
    if (nloosejets>0) {dp0=getDeltaPhi( loosejetPhi->at(0), phi_of_MET);} else {return false;}
    if (nloosejets>1) {dp1=getDeltaPhi( loosejetPhi->at(1), phi_of_MET);} else {return false;}
    if (nloosejets>2) {dp2=getDeltaPhi( loosejetPhi->at(2), phi_of_MET);} else {return false;}
    //here is the implementation of the DeltaPhi cuts
    if ( dp0 >0.3 && dp1 >0.5 && dp2 >0.3 ) { return true; } else {return false;}
  }
  else if (cutIndex == cutDeltaPhi && //replace normal DeltaPhi cut with minDeltaPhi cut
	   ( theDPType_ == kminDP ) ) {
    
    if (theMETType_ == kMHT)     return ( getMinDeltaPhiMHT(3) >= 0.3 );
    else if (theMETType_ ==kMET) return ( getMinDeltaPhiMET(3) >= 0.3 );
    else if (theMETType_ ==ktcMET) return ( getMinDeltaPhitcMET(3) >= 0.3 );
    else {assert(0);}
  }
  else if (cutIndex == cutDeltaPhi && //replace normal DeltaPhi cut with DeltaPhi(MPT,MET) cut
	   ( theDPType_ == kMPT ) ) {
    if (theMETType_ == kMET)    return ( getDeltaPhiMPTMET() < 2.0 );
    else {cout<<"DeltaPhiMPTMET not implemented for that MET type!"<<endl; assert(0);}
  }

  //in case this is not an exception, return the cut result stored in the ntuple
  return cutResults->at(cutIndex);
}

Int_t basicLoop::Cut(Long64_t entry)
{
  
  for (unsigned int i=0; i< cutResults->size(); i++) {
    if (cutRequired(i) && !passCut(i) ) return -1;
  }
  
  return 1;
  
}

double basicLoop::getUncorrectedHT(const double threshold) {

  if (!isVersion04b()) return 0; //doesn't work in older versions

  /*
I have realized that if I want to do a crude trigger simulation, then I need to store every jet
with Uncorrect Pt>20 GeV, not just jets passing my 'loose' cuts

   vector<float>   *loosejetPtUncorr;
   vector<float>   *loosejetEtaUncorr;
   vector<float>   *loosejetPhiUncorr;*/
  double ht=0;
  for (unsigned int ijet=0; ijet<veryloosejetPtUncorr->size(); ijet++) {
    double thispt =fabs(veryloosejetPtUncorr->at(ijet));
    if (thispt>=threshold)   ht+= thispt;
  }
  
   return ht;
}

double basicLoop::getMinDeltaPhiMET(unsigned int maxjets) {

  unsigned int njets=  jetPhi->size();

  if (njets < maxjets) maxjets = njets;

  //get the minimum angle between the first n jets and MET
  double mindp=99;
  for (unsigned int i=0; i< maxjets; i++) {

    double dp =  acos(cos( jetPhi->at(i) - METphi));
    if (dp<mindp) mindp=dp;

  }
  return mindp;
}

double basicLoop::getMinDeltaPhitcMET(unsigned int maxjets) {

  unsigned int njets=  jetPhi->size();

  if (njets < maxjets) maxjets = njets;

  //get the minimum angle between the first n jets and MET
  double mindp=99;
  for (unsigned int i=0; i< maxjets; i++) {

    double dp =  acos(cos( jetPhi->at(i) - tcMETphi));
    if (dp<mindp) mindp=dp;

  }
  return mindp;
}

double basicLoop::getMinDeltaPhiMHT(unsigned int maxjets) {

  unsigned int njets=  jetPhi->size();

  if (njets < maxjets) maxjets = njets;

  //get the minimum angle between the first n jets and MET
  double mindp=99;
  for (unsigned int i=0; i< maxjets; i++) {

    double dp =  acos(cos( jetPhi->at(i) - MHTphi));
    if (dp<mindp) mindp=dp;

  }
  return mindp;
}

double basicLoop::getMinDeltaPhibMET() {
  //get the minimum angle between a b jet and MET
  double mindp=99;
  for (unsigned int i=0; i<jetPhi->size(); i++) {
    if (passBCut(i) ) {
      double dp =  getDeltaPhi( jetPhi->at(i), METphi);
      if (dp<mindp) mindp=dp;
    }
  }
  return mindp;
}

double basicLoop::getDeltaPhib1b2() {
  //get the angle between the lead 2 b jets
  std::vector<float> phis;

  for (unsigned int i=0; i<jetPhi->size(); i++) {
    if (passBCut(i) ) {
      
      phis.push_back( jetPhi->at(i));

      if (phis.size() == 2) break;
      
    }
  }

  return getDeltaPhi(phis.at(0),phis.at(1));
}

double basicLoop::getDeltaPhiMPTMET() {

   //find MPT
   double MPTx=0;
   double MPTy=0;
   //we have already made minimal cuts on track pt, eta
   for (unsigned int i=0; i< trackPt->size(); i++) {
     
     MPTx -= trackPt->at(i) * cos(trackPhi->at(i));
     MPTy -= trackPt->at(i) * sin(trackPhi->at(i));
   }

   double MPTphi = atan2(MPTy,MPTx);

   return getDeltaPhi(METphi, MPTphi);
}

bool basicLoop::passBCut( unsigned int bindex) {
  //bindex is taken to be the index of a *tight* jet

  //use tagger simpleSecondaryVertexHighEffBJetTags
  //which in older samples is called simpleSecondaryVertexBJetTags

  bool pass=false;
  float oldval=jetBTagDisc_simpleSecondaryVertexBJetTags->at(bindex);
  float newval=jetBTagDisc_simpleSecondaryVertexHighEffBJetTags->at(bindex);

  //if the tagger is not valid, there seems to be a value -1000 returned by CMSSW
  //so clearly that will fail the cut
  if ( (oldval >=1.74) || (newval >=1.74)) pass=true;

  return  pass;
}

double basicLoop::getDeltaPhi(double phi1, double phi2) {

  return acos(cos(phi1-phi2));
}

double basicLoop::getMinDeltaPhi_bj(unsigned int bindex) {

  // for now we don't have a choice -- we can only use tight jets
  //because i left the loose ones out of the ntuple

  double bphi = jetPhi->at(bindex);

  double minDeltaPhi=99;
  //loop over the jets
  for (unsigned int jindex=0; jindex<jetPhi->size(); jindex++) {
    if ( !passBCut(jindex) ) { //only look at non-b jets
      double dp = getDeltaPhi(bphi, jetPhi->at(jindex));
      if (dp < minDeltaPhi) minDeltaPhi = dp;
    }
  }

  return minDeltaPhi;
}

double basicLoop::getOverallMinDeltaR_bj() {

  double minDeltaR_bj=999;
  //note that all tight jet vectors should have the same size
  for (unsigned int ib = 0; ib< jetPhi->size(); ib++) {
    if ( passBCut(ib)) { //refind the b jets
      double mdr=getMinDeltaR_bj(ib); 
      if (mdr<minDeltaR_bj) minDeltaR_bj=mdr;
    }
  }

  return minDeltaR_bj;
}

double basicLoop::getMinDeltaR_bj(unsigned int bindex) {

  double beta=jetEta->at(bindex);
  double bphi=jetPhi->at(bindex);

  double minDeltaR = 999;

  //loop over the jets
  for (unsigned int jindex=0; jindex<jetPhi->size(); jindex++) {
    if ( !passBCut(jindex) ) { //only look at non-b jets
      double dr = pow(getDeltaPhi(bphi, jetPhi->at(jindex)),2) + pow(beta - jetEta->at(jindex),2);
      dr = sqrt(dr);
      if (dr < minDeltaR) minDeltaR = dr;
    }
  }

  return minDeltaR;
}

int basicLoop::getTopDecayCategory() {
  /*
    this is a kludge on top of a kludge, I suppose.
    reclassifying one set of arbitrary integers with another set of arbitrary integers.
    ugly, but i don't have a better idea.
    i will use enums to make things a little more clear
  */

  //here are the codes defined in the ntuple
  enum TopDecayCodes {kUnknown=-1, kNoB = 0, kHadronic = 1, kElectron = 2, kMuon=3,kTauHadronic=4,kTauElectron=5,kTauMuon=6,kTauMisc=7 };
  //  TopDecayCodes mycode=kUnknown;

  //the set of codes returned here is defined in the class defn at the top of this file

  int code=-1;

  unsigned int ntop=  topDecayCode->size();

  if (ntop==2) { //this is what we expect for ttbar
    int code1 = topDecayCode->at(0);
    int code2 = topDecayCode->at(1);

    //i'm pretty sure that the 'misc' category is only filled by hadronic decays
    if (code1==kTauMisc) code1=kHadronic;
    if (code2==kTauMisc) code2=kHadronic;
    
    if ( (code1==kElectron ||code1==kMuon) && (code2==kElectron ||code2==kMuon) ) code=kAllLeptons;
    else if ( code1==kHadronic && code2==kHadronic ) code=kAllHadronic;
    else if ( (code1==kElectron && code2==kHadronic) || (code2==kElectron && code1==kHadronic) ) code=kOneElectron;
    else if ( (code1==kMuon && code2==kHadronic) || (code2==kMuon && code1==kHadronic) ) code=kOneMuon;
    else if ( (code1==kTauElectron && code2==kHadronic) || (code2==kTauElectron && code1==kHadronic) ) code=kOneTauE;
    else if ( (code1==kTauMuon && code2==kHadronic) || (code2==kTauMuon && code1==kHadronic) ) code=kOneTauMu;
    else if ( (code1==kTauHadronic && code2==kHadronic) || (code2==kTauHadronic && code1==kHadronic) ) code=kOneTauHadronic;
    //this logic depends on all of the highest categories being tau
    else if ( code1>=kTauHadronic && code2>=kTauHadronic ) code=kAllTau;
    else if ( code1>=kTauHadronic && (code2==kElectron || code2==kMuon) ) code=kTauPlusLepton;
    else if ( code2>=kTauHadronic && (code1==kElectron || code1==kMuon) ) code=kTauPlusLepton;
    else {
      code=kTTbarUnknown;
      std::cout<<"TopCode = "<<code1<<"\t"<<code2<<std::endl;
    }
  }
  //  else if (ntop==1) { //maybe we expect this for single top?
  //
  //  }
  //  else {
  //    
  //  }

  return code;
}

TString basicLoop::getSampleName(const TString inname) {
 
  if (inname.Contains("/TTbarJets/") )                     return "TTbarJets";
  else if (inname.Contains("/LM0/"))                       return "LM0";
  else if (inname.Contains("/LM1/"))                       return "LM1";
  else if (inname.Contains("/LM2/"))                       return "LM2";
  else if (inname.Contains("/LM3/"))                       return "LM3";
  else if (inname.Contains("/LM4/"))                       return "LM4";
  else if (inname.Contains("/LM5/"))                       return "LM5";
  else if (inname.Contains("/LM6/"))                       return "LM6";
  else if (inname.Contains("/LM7/"))                       return "LM7";
  else if (inname.Contains("/LM8/"))                       return "LM8";
  else if (inname.Contains("/LM9/"))                       return "LM9";
  else if (inname.Contains("/LM9t175/"))                   return "LM9t175";
  else if (inname.Contains("/LM9p/"))                      return "LM9p";
  else if (inname.Contains("/LM10/"))                      return "LM10";
  else if (inname.Contains("/LM11/"))                      return "LM11";
  else if (inname.Contains("/LM12/"))                      return "LM12";
  else if (inname.Contains("/LM13/"))                      return "LM13";
  else if (inname.Contains("/MoreMSSM/"))                  return "mMSSM";
  else if (inname.Contains("/MoreMSSMv2/"))                return "mMSSMv2";
  else if (inname.Contains("/MoreMSSMv3/"))                return "mMSSMv3";
  else if (inname.Contains("/SingleTop-tChannel/"))        return "SingleTop-tChannel";
  else if (inname.Contains("/SingleTop-tWChannel/"))       return "SingleTop-tWChannel";
  else if (inname.Contains("/QCD-Pt1000toInf-madgraph/"))  return "QCD1000";
  else if (inname.Contains("/QCD-Pt100to250-madgraph/"))   return "QCD100";
  else if (inname.Contains("/QCD-Pt250to500-madgraph/"))   return "QCD250";
  else if (inname.Contains("/QCD-Pt500to1000-madgraph/"))  return "QCD500";
  else if (inname.Contains("/WJets/"))                     return "WJets";
  else if (inname.Contains("/ZJets/"))                     return "ZJets";
  else if (inname.Contains("/Zinvisible/"))                return "Zinvisible";

  else if (inname.Contains("/DATA/"))                return "data"; //need to decide if this works
  
  std::cout<<"Cannot find sample name for this sample!"<<std::endl;
  
  return "";

}

double basicLoop::getCrossSection(const TString inname) {
  //  https://twiki.cern.ch/twiki/bin/view/CMS/ProductionReProcessingSpring10
  if (inname.Contains("/TTbarJets/") )                     return 165;
  else if (inname.Contains("/LM0/"))                       return 38.93;
  else if (inname.Contains("/LM1/"))                       return 4.888;
  else if (inname.Contains("/LM2/"))                       return 0.6027;
  else if (inname.Contains("/LM3/"))                       return 3.438;
  else if (inname.Contains("/LM4/"))                       return 1.879;
  else if (inname.Contains("/LM5/"))                       return 0.4734;
  else if (inname.Contains("/LM6/"))                       return 0.3104;
  else if (inname.Contains("/LM7/"))                       return 1.209;
  else if (inname.Contains("/LM8/"))                       return 0.7300;
  else if (inname.Contains("/LM9/"))                       return 7.134;
  else if (inname.Contains("/LM9t175/"))                   return 4.241;
  else if (inname.Contains("/LM9p/"))                      return 1.653;
  else if (inname.Contains("/LM10/"))                      return 0.04778;
  else if (inname.Contains("/LM11/"))                      return 0.8236;
  else if (inname.Contains("/LM12/"))                      return 4.414;
  else if (inname.Contains("/LM13/"))                      return 6.899;
  else if (inname.Contains("/MoreMSSM/"))                  return 1.73;
  else if (inname.Contains("/MoreMSSMv2/"))                return 2.1;
  else if (inname.Contains("/MoreMSSMv3/"))                return 2.6;
  else if (inname.Contains("/SingleTop-tChannel/"))        return 20.44;
  else if (inname.Contains("/SingleTop-tWChannel/"))       return 10.6;
  else if (inname.Contains("/QCD-Pt1000toInf-madgraph/"))  return 83;
  else if (inname.Contains("/QCD-Pt100to250-madgraph/"))   return 7000000;
  else if (inname.Contains("/QCD-Pt250to500-madgraph/"))   return 171000;
  else if (inname.Contains("/QCD-Pt500to1000-madgraph/"))  return 5200;
  else if (inname.Contains("/WJets/"))                     return 28049;
  else if (inname.Contains("/ZJets/"))                     return 2907;
  else if (inname.Contains("/Zinvisible/"))                return 4500;
  
  else if (inname.Contains("/DATA/"))                return -2;

  std::cout<<"Cannot find cross section for this sample!"<<std::endl; 
  return -1;
  
}

TString basicLoop::findInputName() {
   // ===== do something tricky to figure out which sample this is ======
   TString inname="";
   if (fChain->InheritsFrom(TChain::Class())) {
     if (((TChain*) fChain)->GetListOfFiles()->GetEntries() <1) {
       std::cout<<"Chain seems to have no files!"<<std::endl; return "error";
     }
     inname = ((TChain*) fChain)->GetListOfFiles()->At(0)->GetTitle();
   }
   else {
     inname=  fChain->GetCurrentFile()->GetName();
   }
   //this is the file name of the first file in the chain

   return inname;
}

bool basicLoop::setCutScheme(CutScheme cutscheme) {

  if (cutscheme == nCutSchemes) return false;
  theCutScheme_ = cutscheme;
  return true;
}

void basicLoop::setMETType(METType mettype) {

  theMETType_ = mettype;
}

void basicLoop::setMETRange(METRange metrange) {

  theMETRange_ = metrange;
}

void basicLoop::setDPType(dpType dptype) {

  theDPType_ = dptype;

}
TString basicLoop::getCutDescriptionString() {

  TString cuts = CutSchemeNames_[theCutScheme_];
  cuts += "_";
  cuts += METTypeNames_[theMETType_];
  cuts += METRangeNames_[theMETRange_];
  cuts += "_";
  cuts+= dpTypeNames_[theDPType_];
  for (unsigned int icut=0; icut<ignoredCut_.size() ; icut++) {
    cuts+="_No";
    cuts += jmt::fortranize(cutnames_.at(ignoredCut_.at(icut)));
  }
  return cuts; 
}

TString basicLoop::getBCutDescriptionString() {

  TString thecut="ge";
  thecut += nBcut_;
  thecut+="b";
  return thecut;
}

bool basicLoop::isVersion04b() {
  
  TString fn=findInputName();
  if (fn.Contains("/V00-00-04b/")
      || fn.Contains("/V00-00-05/")) return true;
  
  return false;
  
}

void basicLoop::printState() {
  cout<<"Cut scheme set to:    "<<CutSchemeNames_[theCutScheme_]<<endl;
  cout<<"MET type set to:      "<<METTypeNames_[theMETType_]<<endl;
  cout<<"MET range set to:     "<<METRangeNames_[theMETRange_]<<endl;
  cout<<"DeltaPhi type set to: "<<dpTypeNames_[theDPType_]<<endl;
  cout<<"Requiring at least    "<<nBcut_<<" b tags"<<endl;
  for (unsigned int i = 0; i< ignoredCut_.size() ; i++) {
    cout<<"Will ignore cut:    "<<cutnames_.at(i)<<endl;
  }
  
}

// ========================================== end


#endif // #ifdef basicLoop_cxx
