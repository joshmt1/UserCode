//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 12 14:05:33 2010 by ROOT version 5.22/00d
// from TTree tree/tree
// found on file: /cu1/joshmt/BasicNtuple/V00-00-01/MoreMSSM/Spring10/BasicNtuple.root
//////////////////////////////////////////////////////////

#ifndef basicLoop_h
#define basicLoop_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream>
#include <vector>

//avoid spaces and funny characters!
const char *CutSchemeNames_[]={"RA2","RA2withBtagging"};

class basicLoop {
 public :
  enum CutScheme {kRA2=0, kRA2withB, nCutSchemes};
  CutScheme theCutScheme_;
  std::vector<TString> cutnames_;

  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
   // Declaration of leaf types
   std::vector<std::string>  *SUSY_cutNames;
   std::vector<bool>    *SUSY_cutResults;
   std::vector<bool>    *cutResults;
   std::vector<float>   *jetPt;
   std::vector<float>   *jetEta;
   std::vector<float>   *jetPhi;
   std::vector<float>   *jetBTagSSV;
   Int_t           nbSSVM;
   Float_t         HT;
   Float_t         MHT;
   Float_t         MHTphi;
   Float_t         DeltaPhi_JetMHT1;
   Float_t         DeltaPhi_JetMHT2;
   Float_t         DeltaPhi_JetMHT3;
   Float_t         MET;
   Float_t         METphi;
   std::vector<float>   *trackPt;
   std::vector<float>   *trackEta;
   std::vector<float>   *trackPhi;
   Int_t           SUSY_nb;

   // List of branches
   TBranch        *b_SUSY_cutNames;   //!
   TBranch        *b_SUSY_cutResults;   //!
   TBranch        *b_cutResults;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetBTagSSV;   //!
   TBranch        *b_nbSSVM;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_MHT;   //!
   TBranch        *b_MHTphi;   //!
   TBranch        *b_DeltaPhi_JetMHT1;   //!
   TBranch        *b_DeltaPhi_JetMHT2;   //!
   TBranch        *b_DeltaPhi_JetMHT3;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_METphi;   //!
   TBranch        *b_trackPt;   //!
   TBranch        *b_trackEta;   //!
   TBranch        *b_trackPhi;   //!
   TBranch        *b_SUSY_nb;   //!

   basicLoop(TTree *tree=0);
   virtual ~basicLoop();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual void     compareRA2();
   virtual void     exampleLoop();
   virtual void     nbLoop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   bool setCutScheme(CutScheme cutscheme);

   void cutflow();
   bool cutRequired(unsigned int cutIndex) ;
   TString getSampleName(const TString inname) ;
   double getCrossSection(const TString inname) ;
   TString findInputName() ;
   double getDeltaPhiMPTMET();
};

#endif

#ifdef basicLoop_cxx
basicLoop::basicLoop(TTree *tree) :
  theCutScheme_(kRA2)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/cu1/joshmt/BasicNtuple/V00-00-01/MoreMSSM/Spring10/BasicNtuple.root");
      if (!f) {
         f = new TFile("/cu1/joshmt/BasicNtuple/V00-00-01/MoreMSSM/Spring10/BasicNtuple.root");
         f->cd("/cu1/joshmt/BasicNtuple/V00-00-01/MoreMSSM/Spring10/BasicNtuple.root:/BasicTreeMaker");
      }
      tree = (TTree*)gDirectory->Get("tree");

   }

   std::cout<<"In basicLoop ctor. About to do Init"<<std::endl;
   Init(tree);

   //this may well be stored in the event for now too.
   //think about it
   cutnames_.push_back("Inclusive");
   cutnames_.push_back("Trigger");
   cutnames_.push_back("PV");
   cutnames_.push_back(">= 3 Jets");
   cutnames_.push_back("First Jet");
   cutnames_.push_back("2nd Jet");
   cutnames_.push_back("3rd Jet");
   cutnames_.push_back("HT Cut");
   cutnames_.push_back("MET Cut");
   cutnames_.push_back("MHT Cut");
   cutnames_.push_back("Muon veto");
   cutnames_.push_back("electron veto");
   cutnames_.push_back("Delta Phi");
   cutnames_.push_back(">=1 b");
   cutnames_.push_back(">=2 b");
   cutnames_.push_back(">=3 b");

}

basicLoop::~basicLoop()
{
  //No! stupid ROOT! this leads to a double deletion!
  //   if (!fChain) return;
  //   delete fChain->GetCurrentFile();
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
   SUSY_cutNames = 0;
   SUSY_cutResults = 0;
   cutResults = 0;
   jetPt = 0;
   jetEta = 0;
   jetPhi = 0;
   jetBTagSSV = 0;
   trackPt = 0;
   trackEta = 0;
   trackPhi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("SUSY_cutNames", &SUSY_cutNames, &b_SUSY_cutNames);
   fChain->SetBranchAddress("SUSY_cutResults", &SUSY_cutResults, &b_SUSY_cutResults);
   fChain->SetBranchAddress("cutResults", &cutResults, &b_cutResults);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetBTagSSV", &jetBTagSSV, &b_jetBTagSSV);
   fChain->SetBranchAddress("nbSSVM", &nbSSVM, &b_nbSSVM);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("MHT", &MHT, &b_MHT);
   fChain->SetBranchAddress("MHTphi", &MHTphi, &b_MHTphi);
   fChain->SetBranchAddress("DeltaPhi_JetMHT1", &DeltaPhi_JetMHT1, &b_DeltaPhi_JetMHT1);
   fChain->SetBranchAddress("DeltaPhi_JetMHT2", &DeltaPhi_JetMHT2, &b_DeltaPhi_JetMHT2);
   fChain->SetBranchAddress("DeltaPhi_JetMHT3", &DeltaPhi_JetMHT3, &b_DeltaPhi_JetMHT3);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("METphi", &METphi, &b_METphi);
   fChain->SetBranchAddress("trackPt", &trackPt, &b_trackPt);
   fChain->SetBranchAddress("trackEta", &trackEta, &b_trackEta);
   fChain->SetBranchAddress("trackPhi", &trackPhi, &b_trackPhi);
   fChain->SetBranchAddress("SUSY_nb", &SUSY_nb, &b_SUSY_nb);
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

bool basicLoop::cutRequired(unsigned int cutIndex) {
  //not quite sure what the best (most elegant) way to implement this is....
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
  
  //RA2
  if (theCutScheme_ == kRA2 || theCutScheme_==kRA2withB) {
    if      (cutIndex == 0)  return true;
    else if (cutIndex == 1)  return true;
    else if (cutIndex == 2)  return true;
    else if (cutIndex == 3)  return true;
    else if (cutIndex == 4)  return false;
    else if (cutIndex == 5)  return false;
    else if (cutIndex == 6)  return false;
    else if (cutIndex == 7)  return true;
    else if (cutIndex == 8)  return false;
    else if (cutIndex == 9)  return true;
    else if (cutIndex == 10) return true;
    else if (cutIndex == 11) return true;
    else if (cutIndex == 12) return true;
    else if (cutIndex == 13) return (theCutScheme_==kRA2withB);
    else if (cutIndex == 14) return (theCutScheme_==kRA2withB);
    else if (cutIndex == 15) return (theCutScheme_==kRA2withB);
  }

  //should not get past here!
  assert(0);
  return false;
}

Int_t basicLoop::Cut(Long64_t entry)
{
  
  for (unsigned int i=0; i< cutResults->size(); i++) {
    if (cutRequired(i) && !cutResults->at(i) ) return -1;
  }
  
  return 1;
  
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

   return acos(cos(METphi - MPTphi));
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
  else if (inname.Contains("/SingleTop-tChannel/"))        return "SingleTop-tChannel";
  else if (inname.Contains("/SingleTop-tWChannel/"))       return "SingleTop-tWChannel";
  else if (inname.Contains("/QCD-Pt1000toInf-madgraph/"))  return "QCD1000";
  else if (inname.Contains("/QCD-Pt100to250-madgraph/"))   return "QCD100";
  else if (inname.Contains("/QCD-Pt250to500-madgraph/"))   return "QCD250";
  else if (inname.Contains("/QCD-Pt500to1000-madgraph/"))  return "QCD500";
  else if (inname.Contains("/WJets/"))                     return "WJets";
  else if (inname.Contains("/ZJets/"))                     return "ZJets";
  else if (inname.Contains("/Zinvisible/"))                return "Zinvisible";
  
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
  else if (inname.Contains("/SingleTop-tChannel/"))        return 20.44;
  else if (inname.Contains("/SingleTop-tWChannel/"))       return 10.6;
  else if (inname.Contains("/QCD-Pt1000toInf-madgraph/"))  return 83;
  else if (inname.Contains("/QCD-Pt100to250-madgraph/"))   return 7000000;
  else if (inname.Contains("/QCD-Pt250to500-madgraph/"))   return 171000;
  else if (inname.Contains("/QCD-Pt500to1000-madgraph/"))  return 5200;
  else if (inname.Contains("/WJets/"))                     return 28049;
  else if (inname.Contains("/ZJets/"))                     return 2907;
  else if (inname.Contains("/Zinvisible/"))                return 4500;
  
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

  cout<<"Cut scheme set to: "<<CutSchemeNames_[theCutScheme_]<<endl;

  return true;
}


#endif // #ifdef basicLoop_cxx
