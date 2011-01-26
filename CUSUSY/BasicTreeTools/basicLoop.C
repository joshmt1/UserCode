#define basicLoop_cxx
#include "basicLoop.h"
#include <TH1.h>
#include <TH2.h>
#include <TDatime.h>
#include <TFile.h>
#include <TStyle.h>
#include <fstream>
#include <iomanip>

#include "TMath.h"

// just leave this code here, untouched.
// it can be copied and pasted as an event loop template
void basicLoop::exampleLoop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries(); //jmt: remove Fast

   Long64_t nbytes = 0, nb = 0;
   startTimer();  //keep track of performance
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      if (jentry%1000000==0) checkTimer(jentry,nentries);
      nb = GetEntry(jentry);   nbytes += nb; //use member function GetEntry instead of fChain->

      if (Cut(ientry) < 0) continue; //jmt use cut
   }
   stopTimer(nentries);


}

/*
print a cut flow table
lumi is set in basicLoop.h

This code assumes that each successive cut is a subset of the previous cut.
So it is hard to implement the mutually exclusive ==1b and >=2b categories
*/
void basicLoop::cutflow(bool writeFiles)
{
  printState();

  std::vector<int> npass;
  
  if (fChain == 0) return;
  
  const   double sigma = getCrossSection(findInputName());
  
  Long64_t nentries = fChain->GetEntries(); //jmt: remove Fast
  
  const   double   weight = isData_ ? 1 : lumi * sigma / double(nentries); //calculate weight

  Long64_t nbytes = 0, nb = 0;
  
  LoadTree(0);
  nb = GetEntry(0);   nbytes += nb; //use member function GetEntry instead of fChain->

  std::vector<TString> textfilenames; //for the writeFiles option
  std::vector<ofstream*> textfiles;   //for the writeFiles option
  //end of stuff for writeFiles option
  for (unsigned int i=0 ; i<cutTags_.size(); i++) {
    npass.push_back(0);
    if (writeFiles) {
      TString textfilename="/cu3/joshmt/cutflow.";  //FIXME hard-coded path
      textfilename+=getCutDescriptionString();
      textfilename+=".";    textfilename+=getSampleName(findInputName());
      textfilename+=".";
      textfilename+=cutTags_[i];
      textfilenames.push_back(textfilename);
      textfiles.push_back( new ofstream(textfilename.Data()));
    }
  }

  cout<<"Running..."<<endl;  
  
  startTimer();  //keep track of performance
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    if (jentry%1000000==0) checkTimer(jentry,nentries);
    nb = GetEntry(jentry);   nbytes += nb; //use member function GetEntry instead of fChain->

    //cout<<"== "<<jentry<<endl;
  
    //i hate to do this, but I'm going to put a hack here to check the run number
    //    if (runNumber != 143962) continue;

    for (unsigned int i=0 ; i<cutTags_.size(); i++) {
      //cout<<i<<endl;
      if (cutRequired(cutTags_[i]) && passCut(cutTags_[i]) )   npass.at(i) = npass.at(i) +1;
      else if (cutRequired(cutTags_[i]) && !passCut(cutTags_[i]) ) break;

      //optional code to dump events to file
      if (writeFiles)     *textfiles[i] <<runNumber<<" "<<lumiSection<<" "<<eventNumber<<endl;
      //the structure of this code means we will dump files for cuts that are not required.
      //that's a bit wasteful but i'm not going to worry about it
    }
    
  }
  cout<<endl;
  stopTimer(nentries);

  TString samplename=  getSampleName(findInputName());
  TString outfilename="cutflow."; 
  TString outfilenameU="cutflowUnweighted."; 

  outfilename+=getCutDescriptionString();
  outfilename+=".";    outfilename+=samplename; 
  outfilename+=".dat";

  outfilenameU+=getCutDescriptionString();
  outfilenameU+=".";    outfilenameU+=samplename; 
  outfilenameU+=".dat";
  ofstream file(outfilename.Data());
  ofstream fileU(outfilenameU.Data());
  
  for (unsigned int i=0 ; i<npass.size(); i++) {
    
    if (cutRequired(cutTags_[i])) {
      
      //error on n is sqrt n
      double error = sqrt(npass.at(i));
      double weighted = npass.at(i) * weight;
      double weighted_error = error*weight;
      
      char ccc[150];
      sprintf(ccc,"%20s %15d | %.2f | Weighted = %f +/- %f",cutNames_[cutTags_[i]].Data(),npass.at(i),100*double(npass.at(i))/double(npass.at(0)),weighted,weighted_error);
      cout<<ccc<<endl;

      //now including a decription string too. not perfect but better than nothing
      file <<cutNames_[cutTags_[i]].Data()<<"\t"<<setprecision(20) << weighted<<"\t" << weighted_error<<endl;
      fileU <<cutNames_[cutTags_[i]].Data()<<"\t"<<setprecision(20) << npass.at(i)<<endl;
    }
  }
  
  file.close();
  fileU.close();
  if (writeFiles) {
    for (unsigned int i=0 ; i<cutTags_.size(); i++) {
      textfiles[i]->close(); //i might be leaking the pointers here, but i don't care
    }
  }

}

/*
ABCD tree maker
*/
void basicLoop::ABCDtree(unsigned int dataindex)
{
  if (fChain == 0) return;
  printState();

  
  Long64_t nentries = fChain->GetEntries(); //jmt: remove Fast

  TString inname=findInputName(); //uses fChain
  std::cout<<"Got an input file name as: "<<inname<<std::endl;
  
  double sigma = getCrossSection(inname);
  TString sampleName = getSampleName(inname);
  if (sigma<=0 && !isData_) return;

  if (nentries == 0) {std::cout<<"Chain has no entries!"<<std::endl; return;}
  
  double   weight = isData_ ? 1 : lumi * sigma / double(nentries); //calculate weight

  //open output file
  //FIXME hardcoded for dellcmscornell here
  //TString outfilename="/cu1/joshmt/ABCDtrees/ABCDtree.";
  TString outfilename="/cu1/kreis/ABCDtrees/36_Jan20/ABCDtree.";
  outfilename+=getCutDescriptionString();
  outfilename+=".";    outfilename+=getBCutDescriptionString(); 
  outfilename+=".";    outfilename+=sampleName; 
  if (isData_) {
    outfilename+="-";
    outfilename+=dataindex;
  }
  outfilename+=".root";
  TFile fout(outfilename,"RECREATE");
  
  // == make ABCD tree ==
  double myHT;
  double myMET;
  double myMHT;
  double minDeltaPhiMET;
  double minDeltaPhiMETAll;
  double minDeltaPhiMET30All;
  double minDeltaPhiMET30_eta5All;
  double minDeltaPhiMET30_eta5_noIdAll;
  double minDeltaPhiMHT;
  double minDeltaRbj;
  double DeltaPhiMPTMET;

  TTree ABCDtree("ABCDtree","ABCD tree");
  ABCDtree.Branch("weight",&weight,"weight/D");
  ABCDtree.Branch("HT",&myHT,"HT/D");
  ABCDtree.Branch("MET",&myMET,"MET/D");
  ABCDtree.Branch("MHT",&myMHT,"MHT/D");
  ABCDtree.Branch("minDeltaPhiMET",&minDeltaPhiMET,"minDeltaPhiMET/D");
  ABCDtree.Branch("minDeltaPhiMETAll",&minDeltaPhiMETAll,"minDeltaPhiMETAll/D");
  ABCDtree.Branch("minDeltaPhiMET30All",&minDeltaPhiMET30All,"minDeltaPhiMET30All/D");
  ABCDtree.Branch("minDeltaPhiMET30_eta5All",&minDeltaPhiMET30_eta5All,"minDeltaPhiMET30_eta5All/D");
  ABCDtree.Branch("minDeltaPhiMET30_eta5_noIdAll",&minDeltaPhiMET30_eta5_noIdAll,"minDeltaPhiMET30_eta5_noIdAll/D");
  ABCDtree.Branch("minDeltaPhiMHT",&minDeltaPhiMHT,"minDeltaPhiMHT/D");
  ABCDtree.Branch("minDeltaRbj",&minDeltaRbj,"minDeltaRbj/D");
  ABCDtree.Branch("DeltaPhiMPTMET",&DeltaPhiMPTMET,"DeltaPhiMPTMET/D");
  ABCDtree.Branch("nbSSVM", &nbSSVM, "nbSSVM/I");

  startTimer();  //keep track of performance
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;     
    if (jentry%1000000==0) checkTimer(jentry,nentries);
    nb = GetEntry(jentry);   nbytes += nb; //use member function GetEntry instead of fChain->

    if (Cut(ientry) < 0) continue; //jmt use cut
    myHT = getHT();
    myMET = getMET();
    myMHT = getMHT();
    minDeltaPhiMET = getMinDeltaPhiMET(3);
    minDeltaPhiMETAll = getMinDeltaPhiMET(99);
    minDeltaPhiMET30All = getMinDeltaPhiMET30(99);
    minDeltaPhiMET30_eta5All = getMinDeltaPhiMET30_eta5(99);
    minDeltaPhiMET30_eta5_noIdAll = getMinDeltaPhiMET30_eta5_noId(99);
    
    //MHT is just a type of MET
    METType userMETType = theMETType_;
    setMETType( kMHT );
    minDeltaPhiMHT = getMinDeltaPhiMET(3);
    setMETType( userMETType);

    minDeltaRbj = getOverallMinDeltaR_bj();
    DeltaPhiMPTMET = getDeltaPhiMPTMET();

    ABCDtree.Fill(); 
  }
  stopTimer(nentries);

  fout.Write();
  fout.Close();
  
}

/*
experimental new loop -- for use with Baseline0 cut scheme

this is the code used for the normalized distributions in the AN
*/
void basicLoop::cutflowPlotter()
//no data index -- data must be passed in as one big chain!
{
 
   if (fChain == 0) return;
   printState();

   Long64_t nentries = fChain->GetEntries(); //jmt: remove Fast

   std::cout<<"Got an input file name as: "<<findInputName()<<std::endl;
   
   double sigma = getCrossSection();
   TString sampleName = getSampleName();
   if (sigma<=0 && !isData_) {std::cout<<"Cross section problem! "<<sigma<<std::endl; return;}

   if (nentries == 0) {std::cout<<"Chain has no entries!"<<std::endl; return;}

   double   weight = isData_ ? 1 : lumi * sigma / double(nentries); //calculate weight

   //open output file
   TString outfilename="cutflowPlots.";
   outfilename+=getCutDescriptionString();
   outfilename+=".";    outfilename+=sampleName; 
   outfilename+=".root";
   TFile fout(outfilename,"RECREATE");

   TH1::SetDefaultSumw2(); //trick to turn on Sumw2 for all histos

   //make histograms (and call sumw2)
   TString name;
   std::map<TString, TH1D*> H_HT;
   std::map<TString, TH1D*> H_NJets;
   std::map<TString, TH1D*> H_NElectrons;
   std::map<TString, TH1D*> H_NMuons;
   std::map<TString, TH1D*> H_MET;
   std::map<TString, TH1D*> H_minDeltaPhi;
   std::map<TString, TH1D*> H_NBJets;
   std::map<TString, TH1D*> H_topDecayCategory;

   for (unsigned int i=0 ; i<cutTags_.size(); i++) {
     if (cutRequired(cutTags_[i])  ) {
       //name will indicate that the histogram is *after* the cut in the name
       name="H_HT_";    name += cutTags_[i];
       H_HT[name]=new TH1D(name,name,300,0,1200);

       name="H_NJets_";    name += cutTags_[i];
       H_NJets[name]=new TH1D(name,name,10,0,10);

       name="H_NBJets_";    name += cutTags_[i];
       H_NBJets[name]=new TH1D(name,name,10,0,10);

       name="H_NElectrons_"; name += cutTags_[i];
       H_NElectrons[name]=new TH1D(name,name,10,0,10);

       name="H_NMuons_"; name += cutTags_[i];
       H_NMuons[name]=new TH1D(name,name,10,0,10);

       name="H_MET_"; name += cutTags_[i];
       H_MET[name]=new TH1D(name,name,100,0,400);

       name="H_minDeltaPhi_"; name += cutTags_[i];
       H_minDeltaPhi[name]=new TH1D(name,name,50,0,TMath::Pi());

       //only useful for top, obviously
       name="H_topDecayCategory_"; name += cutTags_[i];
       H_topDecayCategory[name]=new TH1D(name,name,nTopCategories-1,0.5,nTopCategories-0.5);
     }
   }
      

   //keep track of performance
   startTimer();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries ;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) {assert(0);}
      nb = GetEntry(jentry);   nbytes += nb; //use member function GetEntry instead of fChain->

      //calculate some quantities
      //unlike in the old ::Loop, here be sure to plot the same quanities used
      //by the Baseline0 scheme when making cuts
      //therefore, this is only certain to work in Baseline0
      float ht = getHT();
    
      int njets = nGoodJets();
      int nbjets = countBJets();
      int nelectrons = countEleSync1();
      int nmuons = countMuSync1();
      float met = getMET() ;
      float minDeltaPhi = getMinDeltaPhiMET(3);

      int topcat = isData_ ? -99 : getTopDecayCategory(); //no sense in looking at MC truth in the data

      //don't call the Cut method but rather do the cut flow step by step
      //just like in ::cutflow()  
      for (unsigned int i=0 ; i<cutTags_.size(); i++) {
	if (cutRequired(cutTags_[i]) && passCut(cutTags_[i]) ) {
	  //event passes that step of cut flow!
	  name="H_HT_";	  name += cutTags_[i];
	  H_HT[name]->Fill(ht,weight);

	  name="H_NJets_";	  name += cutTags_[i];
	  H_NJets[name]->Fill(njets,weight);

	  name="H_NBJets_";	  name += cutTags_[i];
	  H_NBJets[name]->Fill(nbjets,weight);

	  name="H_NElectrons_";	  name += cutTags_[i];
	  H_NElectrons[name]->Fill(nelectrons,weight);

	  name="H_NMuons_";	  name += cutTags_[i];
	  H_NMuons[name]->Fill(nmuons,weight);

	  name="H_MET_";	  name += cutTags_[i];
	  H_MET[name]->Fill(met,weight);

	  name="H_minDeltaPhi_";	  name += cutTags_[i];
	  H_minDeltaPhi[name]->Fill(minDeltaPhi,weight);

	  name="H_topDecayCategory_";	  name += cutTags_[i];
	  H_topDecayCategory[name]->Fill(topcat,weight);
	}
	else if (cutRequired(cutTags_[i]) && !passCut(cutTags_[i]) ) break;
      }
     
   }
   stopTimer(nentries);
 

   fout.Write();
   fout.Close();

}

 //will copy and paste a lot of code from ::Loop here
 //i want to start fresh!
//This is the code used for the N-1 and other stacked plots in the AN
void basicLoop::Nminus1plots()
{
  
  if (fChain == 0) return;
  
  resetIgnoredCut(); //this code is not going to allow the user to specify an ignored cut!
  setBCut(0);
  printState();
  
  Long64_t nentries = fChain->GetEntries(); //jmt: remove Fast
  
  
  TString inname=findInputName(); //uses fChain
  std::cout<<"Got an input file name as: "<<inname<<std::endl;
  
  double sigma = getCrossSection(inname);
  TString sampleName = getSampleName(inname);
  if (sigma<=0 && !isData_) return;
  
  if (nentries == 0) {std::cout<<"Chain has no entries!"<<std::endl; return;}
  
  double   weight = isData_ ? 1 : lumi * sigma / double(nentries); //calculate weight
  
  //open output file
  TString outfilename="Nminus1plots."; 
  outfilename+=getCutDescriptionString();
  outfilename+=".";    outfilename+=sampleName; 
  outfilename+=".root";
  TFile fout(outfilename,"RECREATE");
  
  //TH1::SetDefaultSumw2(); //turn on sumw2 for all histos
  
  //make some histograms (each with >=1,2 ==1 b tags)
  int nbins = 100; //can always rebin
  //HT
  int multiplier=3;
  double ht_min=0,ht_max=1200;
  TH1D H_HT("H_HT","HT",nbins*multiplier,ht_min,ht_max);
  TH1D H_HT_ge1b("H_HT_ge1b","HT (>=1b)",nbins*multiplier,ht_min,ht_max);
  TH1D H_HT_ge2b("H_HT_ge2b","HT (>=2b)",nbins*multiplier,ht_min,ht_max);
  TH1D H_HT_eq1b("H_HT_eq1b","HT (==1b)",nbins*multiplier,ht_min,ht_max);

   //n jets
   int njets_bins=10;
   TH1D Hnjets("Hnjets","N of jets (N-1)",njets_bins,0,njets_bins);
   TH1D Hnjets_ge1b("Hnjets_ge1b","N of jets (N-1 && >=1b)",njets_bins,0,njets_bins);
   TH1D Hnjets_ge2b("Hnjets_ge2b","N of jets (N-1 && >=2b)",njets_bins,0,njets_bins);
   TH1D Hnjets_eq1b("Hnjets_eq1b","N of jets (N-1 && ==1b)",njets_bins,0,njets_bins);

   //some electron variables (do this for an IMV/IEV region only)
//    double l_pt_min=0,l_pt_max=100;
//    TH1D HelePt("HelePt","electron pT (N-1)",nbins,l_pt_min,l_pt_max);
//    TH1D HeleEta("HeleEta","electron eta (N-1)",nbins,-5,5);
//    TH1D HeleIso("HeleIso","electron isolation (N-1)",nbins,0,5);

//    TH1D HelePt_ge1b("HelePt_ge1b","electron pT (N-1)",nbins,l_pt_min,l_pt_max);
//    TH1D HeleEta_ge1b("HeleEta_ge1b","electron eta (N-1)",nbins,-5,5);
//    TH1D HeleIso_ge1b("HeleIso_ge1b","electron isolation (N-1)",nbins,0,5);

//    TH1D HelePt_ge2b("HelePt_ge2b","electron pT (N-1)",nbins,l_pt_min,l_pt_max);
//    TH1D HeleEta_ge2b("HeleEta_ge2b","electron eta (N-1)",nbins,-5,5);
//    TH1D HeleIso_ge2b("HeleIso_ge2b","electron isolation (N-1)",nbins,0,5);

//    //some muon variables (pT and isolation)
//    TH1D HmuonPt("HmuonPt","muon pT (N-1)",nbins,l_pt_min,l_pt_max);
//    TH1D HmuonEta("HmuonEta","muon eta (N-1)",nbins,-5,5);
//    TH1D HmuonIso("HmuonIso","muon isolation (N-1)",nbins,0,5);

//    TH1D HmuonPt_ge1b("HmuonPt_ge1b","muon pT (N-1)",nbins,l_pt_min,l_pt_max);
//    TH1D HmuonEta_ge1b("HmuonEta_ge1b","muon eta (N-1)",nbins,-5,5);
//    TH1D HmuonIso_ge1b("HmuonIso_ge1b","muon isolation (N-1)",nbins,0,5);

//    TH1D HmuonPt_ge2b("HmuonPt_ge2b","muon pT (N-1)",nbins,l_pt_min,l_pt_max);
//    TH1D HmuonEta_ge2b("HmuonEta_ge2b","muon eta (N-1)",nbins,-5,5);
//    TH1D HmuonIso_ge2b("HmuonIso_ge2b","muon isolation (N-1)",nbins,0,5);
   TH1D HnElectrons("HnElectrons","N of Electrons (N-1)",njets_bins,0,njets_bins);
   TH1D HnElectrons_ge1b("HnElectrons_ge1b","N of Electrons (N-1 && >=1b)",njets_bins,0,njets_bins);
   TH1D HnElectrons_ge2b("HnElectrons_ge2b","N of Electrons (N-1 && >=2b)",njets_bins,0,njets_bins);
   TH1D HnElectrons_eq1b("HnElectrons_eq1b","N of Electrons (N-1 && ==1b)",njets_bins,0,njets_bins);

   TH1D HnMuons("HnMuons","N of Muons (N-1)",njets_bins,0,njets_bins);
   TH1D HnMuons_ge1b("HnMuons_ge1b","N of Muons (N-1 && >=1b)",njets_bins,0,njets_bins);
   TH1D HnMuons_ge2b("HnMuons_ge2b","N of Muons (N-1 && >=2b)",njets_bins,0,njets_bins);
   TH1D HnMuons_eq1b("HnMuons_eq1b","N of Muons (N-1 && ==1b)",njets_bins,0,njets_bins);

   //MET
   double met_max=500,met_min=0;
   TH1D H_MET("H_MET","MET",nbins, met_min , met_max);
   TH1D H_METphi("H_METphi","MET phi",nbins, -TMath::Pi() , TMath::Pi());

   TH1D H_MET_ge1b("H_MET_ge1b","MET (RA2 && >=1b)",nbins, met_min , met_max);
   TH1D H_METphi_ge1b("H_METphi_ge1b","MET phi (>=1b)",nbins, -TMath::Pi() , TMath::Pi());

   TH1D H_MET_ge2b("H_MET_ge2b","MET (RA2 && >=2b)",nbins, met_min , met_max);
   TH1D H_METphi_ge2b("H_METphi_ge2b","MET phi (>=2b)",nbins, -TMath::Pi() , TMath::Pi());

   TH1D H_MET_eq1b("H_MET_eq1b","MET (RA2 && ==1b)",nbins, met_min , met_max);
   TH1D H_METphi_eq1b("H_METphi_eq1b","MET phi (==1b)",nbins, -TMath::Pi() , TMath::Pi());

   //minDeltaPhi
   TH1D HminDeltaPhiMETj("HminDeltaPhiMETj","minDeltaPhi(j123,MET) (N-1)",nbins,0,TMath::Pi());
   TH1D HminDeltaPhiMETj_ge1b("HminDeltaPhiMETj_ge1b","minDeltaPhi(j123,MET) (N-1)",nbins,0,TMath::Pi());
   TH1D HminDeltaPhiMETj_ge2b("HminDeltaPhiMETj_ge2b","minDeltaPhi(j123,MET) (N-1)",nbins,0,TMath::Pi());
   TH1D HminDeltaPhiMETj_eq1b("HminDeltaPhiMETj_eq1b","minDeltaPhi(j123,MET) (N-1)",nbins,0,TMath::Pi());

   TH1D HminDeltaPhiMETjAll("HminDeltaPhiMETjAll","minDeltaPhi(all jets,MET) (N-1)",nbins,0,TMath::Pi());
   TH1D HminDeltaPhiMETjAll_ge1b("HminDeltaPhiMETjAll_ge1b","minDeltaPhi(all jets,MET) (N-1)",nbins,0,TMath::Pi());
   TH1D HminDeltaPhiMETjAll_ge2b("HminDeltaPhiMETjAll_ge2b","minDeltaPhi(all jets,MET) (N-1)",nbins,0,TMath::Pi());
   TH1D HminDeltaPhiMETjAll_eq1b("HminDeltaPhiMETjAll_eq1b","minDeltaPhi(all jets,MET) (N-1)",nbins,0,TMath::Pi());

   TH1D HminDeltaPhiMETjAll30("HminDeltaPhiMETjAll30","minDeltaPhi(all jets,MET) (N-1)",nbins,0,TMath::Pi());
   TH1D HminDeltaPhiMETjAll30_ge1b("HminDeltaPhiMETjAll30_ge1b","minDeltaPhi(all jets,MET) (N-1)",nbins,0,TMath::Pi());
   TH1D HminDeltaPhiMETjAll30_ge2b("HminDeltaPhiMETjAll30_ge2b","minDeltaPhi(all jets,MET) (N-1)",nbins,0,TMath::Pi());

   int vnbins=8;
   double vbins[]={0, 0.15, 0.3, 0.5, 0.7, 1, 1.5, 2, TMath::Pi()};
   TH1D HVminDeltaPhiMETj("HVminDeltaPhiMETj","minDeltaPhi(j,MET) (RA2)",vnbins,vbins);
   TH1D HVminDeltaPhiMETj_ge1b("HVminDeltaPhiMETj_ge1b","minDeltaPhi(j,MET) (RA2 && >=1b)",vnbins,vbins);
   TH1D HVminDeltaPhiMETj_ge2b("HVminDeltaPhiMETj_ge2b","minDeltaPhi(j,MET) (RA2 && >=2b)",vnbins,vbins);
   TH1D HVminDeltaPhiMETj_eq1b("HVminDeltaPhiMETj_eq1b","minDeltaPhi(j,MET) (RA2 && ==1b)",vnbins,vbins);
   int vnbins2=2;
   double vbins2[]={0,  0.3, TMath::Pi()};
   TH1D HV2minDeltaPhiMETj("HV2minDeltaPhiMETj","minDeltaPhi(j,MET) (RA2)",vnbins2,vbins2);
   TH1D HV2minDeltaPhiMETj_ge1b("HV2minDeltaPhiMETj_ge1b","minDeltaPhi(j,MET) (RA2 && >=1b)",vnbins2,vbins2);
   TH1D HV2minDeltaPhiMETj_ge2b("HV2minDeltaPhiMETj_ge2b","minDeltaPhi(j,MET) (RA2 && >=2b)",vnbins2,vbins2);
   TH1D HV2minDeltaPhiMETj_eq1b("HV2minDeltaPhiMETj_eq1b","minDeltaPhi(j,MET) (RA2 && ==1b)",vnbins2,vbins2);

   // ratio of passing minDeltaPhi
   TH1D   HminDeltaPhiRatio("HminDeltaPhiRatio","pass minDeltaPhi / fail minDeltaPhi",nbins, met_min,met_max);
   TH1D   HminDeltaPhiRatio_ge1b("HminDeltaPhiRatio_ge1b","pass minDeltaPhi / fail minDeltaPhi",nbins, met_min,met_max);
   TH1D   HminDeltaPhiRatio_ge2b("HminDeltaPhiRatio_ge2b","pass minDeltaPhi / fail minDeltaPhi",nbins, met_min,met_max);

   TH1D   HminDeltaPhiAllRatio("HminDeltaPhiAllRatio","pass minDeltaPhiAll / fail minDeltaPhiAll",nbins, met_min,met_max);
   TH1D   HminDeltaPhiAllRatio_ge1b("HminDeltaPhiAllRatio_ge1b","pass minDeltaPhiAll / fail minDeltaPhiAll",nbins, met_min,met_max);
   TH1D   HminDeltaPhiAllRatio_ge2b("HminDeltaPhiAllRatio_ge2b","pass minDeltaPhiAll / fail minDeltaPhiAll",nbins, met_min,met_max);

   //histograms just used for calculations
   TH1D   HminDeltaPhiPass("HminDeltaPhiPass","pass minDeltaPhi",nbins, met_min,met_max);
   TH1D   HminDeltaPhiPass_ge1b("HminDeltaPhiPass_ge1b","pass minDeltaPhi",nbins, met_min,met_max);
   TH1D   HminDeltaPhiPass_ge2b("HminDeltaPhiPass_ge2b","pass minDeltaPhi",nbins, met_min,met_max);
   TH1D   HminDeltaPhiFail("HminDeltaPhiFail","fail minDeltaPhi",nbins, met_min,met_max);
   TH1D   HminDeltaPhiFail_ge1b("HminDeltaPhiFail_ge1b","fail minDeltaPhi",nbins, met_min,met_max);
   TH1D   HminDeltaPhiFail_ge2b("HminDeltaPhiFail_ge2b","fail minDeltaPhi",nbins, met_min,met_max);

   TH1D   HminDeltaPhiAllPass("HminDeltaPhiAllPass","pass minDeltaPhiAll",nbins, met_min,met_max);
   TH1D   HminDeltaPhiAllPass_ge1b("HminDeltaPhiAllPass_ge1b","pass minDeltaPhiAll",nbins, met_min,met_max);
   TH1D   HminDeltaPhiAllPass_ge2b("HminDeltaPhiAllPass_ge2b","pass minDeltaPhiAll",nbins, met_min,met_max);
   TH1D   HminDeltaPhiAllFail("HminDeltaPhiAllFail","fail minDeltaPhiAll",nbins, met_min,met_max);
   TH1D   HminDeltaPhiAllFail_ge1b("HminDeltaPhiAllFail_ge1b","fail minDeltaPhiAll",nbins, met_min,met_max);
   TH1D   HminDeltaPhiAllFail_ge2b("HminDeltaPhiAllFail_ge2b","fail minDeltaPhiAll",nbins, met_min,met_max);

   //similar histograms...we want to plot r(MET) for jet undermeasurement, jet overmeasurement
   TH1D   HminDeltaPhiRatio_Over("HminDeltaPhiRatio_Over","pass minDeltaPhi / fail minDeltaPhi",nbins, met_min,met_max);
   TH1D   HminDeltaPhiRatio_Under("HminDeltaPhiRatio_Under","pass minDeltaPhi / fail minDeltaPhi",nbins, met_min,met_max);
   TH1D   HminDeltaPhiRatio_Close("HminDeltaPhiRatio_Close","pass minDeltaPhi / fail minDeltaPhi",nbins, met_min,met_max);

   TH1D   HminDeltaPhiPass_Over("HminDeltaPhiPass_Over","pass minDeltaPhi",nbins, met_min,met_max);
   TH1D   HminDeltaPhiFail_Over("HminDeltaPhiFail_Over","fail minDeltaPhi",nbins, met_min,met_max);
   TH1D   HminDeltaPhiPass_Under("HminDeltaPhiPass_Under","pass minDeltaPhi",nbins, met_min,met_max);
   TH1D   HminDeltaPhiFail_Under("HminDeltaPhiFail_Under","fail minDeltaPhi",nbins, met_min,met_max);

   TH1D   HminDeltaPhiPass_Close("HminDeltaPhiPass_Close","pass minDeltaPhi",nbins, met_min,met_max);
   TH1D   HminDeltaPhiFail_Close("HminDeltaPhiFail_Close","fail minDeltaPhi",nbins, met_min,met_max);

   //MPT , MET
   TH1D HdeltaPhiMPTMET_MET150("HdeltaPhiMPTMET_MET150","DeltaPhi(MPT,MET) MET>150 no mindp",nbins,0,TMath::Pi()+0.01);
   TH1D HdeltaPhiMPTMET_MET200("HdeltaPhiMPTMET_MET200","DeltaPhi(MPT,MET) MET>200 no mindp",nbins,0,TMath::Pi()+0.01);
   TH1D HdeltaPhiMPTMET_MET250("HdeltaPhiMPTMET_MET250","DeltaPhi(MPT,MET) MET>250 no mindp",nbins,0,TMath::Pi()+0.01);

   TH1D   HminDeltaPhiRatio_MPTMET_MET150("HminDeltaPhiRatio_MPTMET_MET150","pass minDeltaPhi / fail minDeltaPhi",nbins,0,TMath::Pi()+0.01);
   TH1D   HminDeltaPhiPass_MPTMET_MET150("HminDeltaPhiPass_MPTMET_MET150","pass minDeltaPhi",nbins,0,TMath::Pi()+0.01);
   TH1D   HminDeltaPhiFail_MPTMET_MET150("HminDeltaPhiFail_MPTMET_MET150","fail minDeltaPhi",nbins,0,TMath::Pi()+0.01 );

   TH1D   HminDeltaPhiRatio_MPTMET_MET200("HminDeltaPhiRatio_MPTMET_MET200","pass minDeltaPhi / fail minDeltaPhi",nbins,0,TMath::Pi()+0.01);
   TH1D   HminDeltaPhiPass_MPTMET_MET200("HminDeltaPhiPass_MPTMET_MET200","pass minDeltaPhi",nbins,0,TMath::Pi()+0.01);
   TH1D   HminDeltaPhiFail_MPTMET_MET200("HminDeltaPhiFail_MPTMET_MET200","fail minDeltaPhi",nbins,0,TMath::Pi()+0.01 );

   TH1D   HminDeltaPhiRatio_MPTMET_MET250("HminDeltaPhiRatio_MPTMET_MET250","pass minDeltaPhi / fail minDeltaPhi",nbins,0,TMath::Pi()+0.01);
   TH1D   HminDeltaPhiPass_MPTMET_MET250("HminDeltaPhiPass_MPTMET_MET250","pass minDeltaPhi",nbins,0,TMath::Pi()+0.01);
   TH1D   HminDeltaPhiFail_MPTMET_MET250("HminDeltaPhiFail_MPTMET_MET250","fail minDeltaPhi",nbins,0,TMath::Pi()+0.01 );

   //n b jets
   TH1D Hnbjets("Hnbjets","N of jets (RA2)",njets_bins,0,njets_bins);

   //also other angular variables (plotted after all other cuts)
   TH1D HdeltaPhiMPTMET("HdeltaPhiMPTMET","DeltaPhi(MET,MPT) (RA2)",nbins,0,TMath::Pi());
   TH1D HdeltaPhiMPTMET_ge1b("HdeltaPhiMPTMET_ge1b","DeltaPhi(MET,MPT) (RA2 && >=1b)",nbins,0,TMath::Pi());
   TH1D HdeltaPhiMPTMET_ge2b("HdeltaPhiMPTMET_ge2b","DeltaPhi(MET,MPT) (RA2 && >=2b)",nbins,0,TMath::Pi());
   TH1D HdeltaPhiMPTMET_eq1b("HdeltaPhiMPTMET_eq1b","DeltaPhi(MET,MPT) (RA2 && ==1b)",nbins,0,TMath::Pi());
   TH1D HdeltaPhib1b2_ge2b("HdeltaPhib1b2_ge2b","DeltaPhi(b1,b2) (>=2b)",nbins,0,TMath::Pi());

   double pt_min=30;
   double pt_max=700;
   TH1D Hjetpt1("Hjetpt1","pT of lead jet",nbins,pt_min,pt_max);
   TH1D Hjetpt1_ge1b("Hjetpt1_ge1b","pT of lead jet (>=1b)",nbins,pt_min,pt_max);
   TH1D Hjetpt1_ge2b("Hjetpt1_ge2b","pT of lead jet (>=2b)",nbins,pt_min,pt_max);
   TH1D Hjetpt1_eq1b("Hjetpt1_eq1b","pT of lead jet (==1b)",nbins,pt_min,pt_max);

   TH1D Hjetphi1("Hjetphi1","phi of lead jet",nbins,-TMath::Pi(),TMath::Pi());
   TH1D Hjetphi1_ge1b("Hjetphi1_ge1b","phi of lead jet (>=1b)",nbins,-TMath::Pi(),TMath::Pi());
   TH1D Hjetphi1_ge2b("Hjetphi1_ge2b","phi of lead jet (>=2b)",nbins,-TMath::Pi(),TMath::Pi());
   TH1D Hjetphi1_eq1b("Hjetphi1_eq1b","phi of lead jet (==1b)",nbins,-TMath::Pi(),TMath::Pi());

   double eta_max = 2.4;
   TH1D Hjeteta1("Hjeteta1","eta of lead jet",nbins,-eta_max,eta_max);
   TH1D Hjeteta1_ge1b("Hjeteta1_ge1b","eta of lead jet (>=1b)",nbins,-eta_max,eta_max);
   TH1D Hjeteta1_ge2b("Hjeteta1_ge2b","eta of lead jet (>=2b)",nbins,-eta_max,eta_max);
   TH1D Hjeteta1_eq1b("Hjeteta1_eq1b","eta of lead jet (==1b)",nbins,-eta_max,eta_max);

   TH1D Hbjetpt1_ge1b("Hbjetpt1_ge1b","pT of lead b jet (>=1b)",nbins,pt_min,pt_max);
   TH1D Hbjetpt1_ge2b("Hbjetpt1_ge2b","pT of lead b jet (>=2b)",nbins,pt_min,pt_max);
   TH1D Hbjetpt1_eq1b("Hbjetpt1_eq1b","pT of lead b jet (==1b)",nbins,pt_min,pt_max);

   TH1D Hbjetphi1_ge1b("Hbjetphi1_ge1b","phi of lead b jet (>=1b)",nbins,-TMath::Pi(),TMath::Pi());
   TH1D Hbjetphi1_ge2b("Hbjetphi1_ge2b","phi of lead b jet (>=2b)",nbins,-TMath::Pi(),TMath::Pi());
   TH1D Hbjetphi1_eq1b("Hbjetphi1_eq1b","phi of lead b jet (==1b)",nbins,-TMath::Pi(),TMath::Pi());

   TH1D Hbjeteta1_ge1b("Hbjeteta1_ge1b","eta of lead b jet (>=1b)",nbins,-eta_max,eta_max);
   TH1D Hbjeteta1_ge2b("Hbjeteta1_ge2b","eta of lead b jet (>=2b)",nbins,-eta_max,eta_max);
   TH1D Hbjeteta1_eq1b("Hbjeteta1_eq1b","eta of lead b jet (==1b)",nbins,-eta_max,eta_max);

   //MC truth distributions
   TH1D HgenInvisibleEnergy_SR("HgenInvisibleEnergy_SR","MC truth invis energy (SR)",nbins,met_min,met_max);
   TH1D HgenInvisibleEnergy_D("HgenInvisibleEnergy_D","MC truth invis energy (high MET control)",nbins,met_min,met_max);
   TH1D HgenInvisibleEnergy_lowMET("HgenInvisibleEnergy_lowMET","MC truth invis energy (low MET, good mindp)",nbins,met_min,met_max);
   TH1D HgenInvisibleEnergy_lowlow("HgenInvisibleEnergy_lowlow","MC truth invis energy (low MET, bad mindp)",nbins,met_min,met_max);

   TH1D HgenInvisibleMissingEnergy_SR("HgenInvisibleMissingEnergy_SR","MC truth invis energy vector (SR)",nbins,met_min,met_max);
   TH1D HgenInvisibleMissingEnergy_D("HgenInvisibleMissingEnergy_D","MC truth invis energy vector (high MET control)",nbins,met_min,met_max);
   TH1D HgenInvisibleMissingEnergy_lowMET("HgenInvisibleMissingEnergy_lowMET","MC truth invis energy vector (low MET, good mindp)",nbins,met_min,met_max);
   TH1D HgenInvisibleMissingEnergy_lowlow("HgenInvisibleMissingEnergy_lowlow","MC truth invis energy vector (low MET, bad mindp)",nbins,met_min,met_max);

   TH1D HgenMET_SR("HgenMET_SR","true MET",nbins,met_min,met_max);
   TH1D HgenMET_D("HgenMET_D","true MET",nbins,met_min,met_max);
   TH1D HgenMET_lowMET("HgenMET_lowMET","true MET",nbins,met_min,met_max);
   TH1D HgenMET_lowlow("HgenMET_lowlow","true MET",nbins,met_min,met_max);

   TH1D HmaxJetRecoError_SR("HmaxJetRecoError_SR","max(pT_reco - pT_gen)",nbins,-met_max,met_max);
   TH1D HmaxJetRecoError_D("HmaxJetRecoError_D","max(pT_reco - pT_gen)",nbins,-met_max,met_max);
   TH1D HmaxJetRecoError_lowMET("HmaxJetRecoError_lowMET","max(pT_reco - pT_gen)",nbins,-met_max,met_max);
   TH1D HmaxJetRecoError_lowlow("HmaxJetRecoError_lowlow","max(pT_reco - pT_gen)",nbins,-met_max,met_max);

   TH1D HmaxDeltaPhiMETjAll_SR("HmaxDeltaPhiMETjAll_SR","maxDeltaPhi(all jets,MET)",nbins,0,TMath::Pi()+0.01);
   TH1D HmaxDeltaPhiMETjAll_D("HmaxDeltaPhiMETjAll_D","maxDeltaPhi(all jets,MET)",nbins,0,TMath::Pi()+0.01);
   TH1D HmaxDeltaPhiMETjAll_lowMET("HmaxDeltaPhiMETjAll_lowMET","maxDeltaPhi(all jets,MET)",nbins,0,TMath::Pi()+0.01);
   TH1D HmaxDeltaPhiMETjAll_lowlow("HmaxDeltaPhiMETjAll_lowlow","maxDeltaPhi(all jets,MET)",nbins,0,TMath::Pi()+0.01);


   H_HT.Sumw2();
   H_HT_ge1b.Sumw2();
   H_HT_ge2b.Sumw2();
   H_HT_eq1b.Sumw2();
   Hnjets.Sumw2();
   Hnjets_ge1b.Sumw2();
   Hnjets_ge2b.Sumw2();
   Hnjets_eq1b.Sumw2();
   HnElectrons.Sumw2();
   HnElectrons_ge1b.Sumw2();
   HnElectrons_ge2b.Sumw2();
   HnElectrons_eq1b.Sumw2();
   HnMuons.Sumw2();
   HnMuons_ge1b.Sumw2();
   HnMuons_ge2b.Sumw2();
   HnMuons_eq1b.Sumw2();
   H_MET.Sumw2();
   H_METphi.Sumw2();
   H_MET_ge1b.Sumw2();
   H_METphi_ge1b.Sumw2();
   H_MET_ge2b.Sumw2();
   H_METphi_ge2b.Sumw2();
   H_MET_eq1b.Sumw2();
   H_METphi_eq1b.Sumw2();
   HminDeltaPhiMETj.Sumw2();
   HminDeltaPhiMETj_ge1b.Sumw2();
   HminDeltaPhiMETj_ge2b.Sumw2();
   HminDeltaPhiMETj_eq1b.Sumw2();
   HminDeltaPhiMETjAll.Sumw2();
   HminDeltaPhiMETjAll_ge1b.Sumw2();
   HminDeltaPhiMETjAll_ge2b.Sumw2();
   HminDeltaPhiMETjAll_eq1b.Sumw2();

   HminDeltaPhiMETjAll30.Sumw2();
   HminDeltaPhiMETjAll30_ge1b.Sumw2();
   HminDeltaPhiMETjAll30_ge2b.Sumw2();

   HVminDeltaPhiMETj.Sumw2();
   HVminDeltaPhiMETj_ge1b.Sumw2();
   HVminDeltaPhiMETj_ge2b.Sumw2();
   HVminDeltaPhiMETj_eq1b.Sumw2();
   HV2minDeltaPhiMETj.Sumw2();
   HV2minDeltaPhiMETj_ge1b.Sumw2();
   HV2minDeltaPhiMETj_ge2b.Sumw2();
   HV2minDeltaPhiMETj_eq1b.Sumw2();
   Hnbjets.Sumw2();
   HdeltaPhiMPTMET.Sumw2();
   HdeltaPhiMPTMET_ge1b.Sumw2();
   HdeltaPhiMPTMET_ge2b.Sumw2();
   HdeltaPhiMPTMET_eq1b.Sumw2();
   HdeltaPhib1b2_ge2b.Sumw2();
   Hjetpt1.Sumw2();
   Hjetpt1_ge1b.Sumw2();
   Hjetpt1_ge2b.Sumw2();
   Hjetpt1_eq1b.Sumw2();
   Hjetphi1.Sumw2();
   Hjetphi1_ge1b.Sumw2();
   Hjetphi1_ge2b.Sumw2();
   Hjetphi1_eq1b.Sumw2();
   Hjeteta1.Sumw2();
   Hjeteta1_ge1b.Sumw2();
   Hjeteta1_ge2b.Sumw2();
   Hjeteta1_eq1b.Sumw2();
   Hbjetpt1_ge1b.Sumw2();
   Hbjetpt1_ge2b.Sumw2();
   Hbjetpt1_eq1b.Sumw2();
   Hbjetphi1_ge1b.Sumw2();
   Hbjetphi1_ge2b.Sumw2();
   Hbjetphi1_eq1b.Sumw2();
   Hbjeteta1_ge1b.Sumw2();
   Hbjeteta1_ge2b.Sumw2();
   Hbjeteta1_eq1b.Sumw2();

   HminDeltaPhiRatio.Sumw2();
   HminDeltaPhiRatio_ge1b.Sumw2();
   HminDeltaPhiRatio_ge2b.Sumw2();

   HminDeltaPhiAllRatio.Sumw2();
   HminDeltaPhiAllRatio_ge1b.Sumw2();
   HminDeltaPhiAllRatio_ge2b.Sumw2();

   HminDeltaPhiRatio_Over.Sumw2();
   HminDeltaPhiRatio_Under.Sumw2();
   HminDeltaPhiRatio_Close.Sumw2();

   HminDeltaPhiPass_Over.Sumw2();
   HminDeltaPhiFail_Over.Sumw2();
   HminDeltaPhiPass_Under.Sumw2();
   HminDeltaPhiFail_Under.Sumw2();
   HminDeltaPhiPass_Close.Sumw2();
   HminDeltaPhiFail_Close.Sumw2();

   //histograms just used for calculations
   HminDeltaPhiPass.Sumw2();
   HminDeltaPhiPass_ge1b.Sumw2();
   HminDeltaPhiPass_ge2b.Sumw2();
   HminDeltaPhiFail.Sumw2();
   HminDeltaPhiFail_ge1b.Sumw2();
   HminDeltaPhiFail_ge2b.Sumw2();

   HminDeltaPhiAllPass.Sumw2();
   HminDeltaPhiAllPass_ge1b.Sumw2();
   HminDeltaPhiAllPass_ge2b.Sumw2();
   HminDeltaPhiAllFail.Sumw2();
   HminDeltaPhiAllFail_ge1b.Sumw2();
   HminDeltaPhiAllFail_ge2b.Sumw2();

   HgenInvisibleEnergy_SR.Sumw2();
   HgenInvisibleEnergy_D.Sumw2();
   HgenInvisibleEnergy_lowMET.Sumw2();
   HgenInvisibleEnergy_lowlow.Sumw2();

   HgenInvisibleMissingEnergy_SR.Sumw2();
   HgenInvisibleMissingEnergy_D.Sumw2();
   HgenInvisibleMissingEnergy_lowMET.Sumw2();
   HgenInvisibleMissingEnergy_lowlow.Sumw2();

   HgenMET_SR.Sumw2();
   HgenMET_D.Sumw2();
   HgenMET_lowMET.Sumw2();
   HgenMET_lowlow.Sumw2();

   HmaxJetRecoError_SR.Sumw2();
   HmaxJetRecoError_D.Sumw2();
   HmaxJetRecoError_lowMET.Sumw2();
   HmaxJetRecoError_lowlow.Sumw2();

   HmaxDeltaPhiMETjAll_SR.Sumw2();
   HmaxDeltaPhiMETjAll_D.Sumw2();
   HmaxDeltaPhiMETjAll_lowMET.Sumw2();
   HmaxDeltaPhiMETjAll_lowlow.Sumw2();

   HdeltaPhiMPTMET_MET150.Sumw2();
   HdeltaPhiMPTMET_MET200.Sumw2();
   HdeltaPhiMPTMET_MET250.Sumw2();
   
   HminDeltaPhiRatio_MPTMET_MET150.Sumw2();
   HminDeltaPhiPass_MPTMET_MET150.Sumw2();
   HminDeltaPhiFail_MPTMET_MET150.Sumw2();
   
   HminDeltaPhiRatio_MPTMET_MET200.Sumw2();
   HminDeltaPhiPass_MPTMET_MET200.Sumw2();
   HminDeltaPhiFail_MPTMET_MET200.Sumw2();
   
   HminDeltaPhiRatio_MPTMET_MET250.Sumw2();
   HminDeltaPhiPass_MPTMET_MET250.Sumw2();
   HminDeltaPhiFail_MPTMET_MET250.Sumw2();
   


   startTimer();   //keep track of performance
   //event loop
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries ;jentry++) {
     if (jentry%1000000==0) checkTimer(jentry,nentries);
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) {assert(0);}
      nb = GetEntry(jentry);   nbytes += nb; //use member function GetEntry instead of fChain->

      //calculate things
      HT = getHT();
      int njets = nGoodJets(); //update to use baseline0
      double minDeltaPhi_j_MET = getMinDeltaPhiMET(3);
      double minDeltaPhi_j_MET_All = getMinDeltaPhiMET(99);
      double minDeltaPhi_j_MET_All30 = getMinDeltaPhiMET30(99);
      double dp_MPTMET = getDeltaPhiMPTMET();
      double deltaPhi_b1b2 =-99; //to be calculated later

      double leadJetPt= jetPtOfN(1);
      double leadJetPhi= jetPhiOfN(1);
      double leadJetEta= jetEtaOfN(1);

      int nelectrons=countEleSync1();
      int nmuons=countMuSync1();

      //HT
      setIgnoredCut("cutHT");
      if (Cut(ientry) >= 0) {
	H_HT.Fill(HT,weight);
	if (nbSSVM >=1)  H_HT_ge1b.Fill(HT,weight);
	if (nbSSVM >=2)  H_HT_ge2b.Fill(HT,weight);
	if (nbSSVM ==1)  H_HT_eq1b.Fill(HT,weight);
      }
      resetIgnoredCut();

      //n jets
      setIgnoredCut("cut3Jets");
      if (Cut(ientry) >= 0) {
	Hnjets.Fill(njets,weight);
	if (nbSSVM >=1)  Hnjets_ge1b.Fill(njets,weight);
	if (nbSSVM >=2)  Hnjets_ge2b.Fill(njets,weight);
	if (nbSSVM ==1)  Hnjets_eq1b.Fill(njets,weight);
      }
      resetIgnoredCut();

      //electron veto
      setIgnoredCut("cutEleVeto");
      if (Cut(ientry) >= 0) {
	HnElectrons.Fill( nelectrons,weight);
	if (nbSSVM >=1)  HnElectrons_ge1b.Fill(nelectrons,weight);
	if (nbSSVM >=2)  HnElectrons_ge2b.Fill(nelectrons,weight);
	if (nbSSVM ==1)  HnElectrons_eq1b.Fill(nelectrons,weight);
      }
      resetIgnoredCut();

      //mu veto
      setIgnoredCut("cutMuVeto");
      if (Cut(ientry) >= 0) {
	HnMuons.Fill( nmuons,weight);
	if (nbSSVM >=1)  HnMuons_ge1b.Fill(nmuons,weight);
	if (nbSSVM >=2)  HnMuons_ge2b.Fill(nmuons,weight);
	if (nbSSVM ==1)  HnMuons_eq1b.Fill(nmuons,weight);
      }
      resetIgnoredCut();

      //MET
      setIgnoredCut("cutMET");
      if (Cut(ientry) >= 0) {
	H_MET.Fill(getMET(),weight);
	if (nbSSVM >=1)  H_MET_ge1b.Fill(getMET(),weight);
	if (nbSSVM >=2)  H_MET_ge2b.Fill(getMET(),weight);
	if (nbSSVM ==1)  H_MET_eq1b.Fill(getMET(),weight);
      }
      resetIgnoredCut();

      //minDeltaPhi
      setIgnoredCut("cutDeltaPhi");
      if (Cut(ientry) >= 0) {
	HminDeltaPhiMETj.Fill( minDeltaPhi_j_MET,weight);
	HminDeltaPhiMETjAll.Fill( minDeltaPhi_j_MET_All,weight);
	HVminDeltaPhiMETj.Fill( minDeltaPhi_j_MET,weight);
	HV2minDeltaPhiMETj.Fill( minDeltaPhi_j_MET,weight);

	HminDeltaPhiMETjAll30.Fill( minDeltaPhi_j_MET_All30,weight);

	//MET > 150
	HdeltaPhiMPTMET_MET150.Fill(dp_MPTMET,weight);
	if (minDeltaPhi_j_MET > 0.3) HminDeltaPhiPass_MPTMET_MET150.Fill(dp_MPTMET,weight);
	else HminDeltaPhiFail_MPTMET_MET150.Fill(dp_MPTMET,weight);

	if (getMET() >200) {
	  HdeltaPhiMPTMET_MET200.Fill(dp_MPTMET,weight);
	  if (minDeltaPhi_j_MET > 0.3) HminDeltaPhiPass_MPTMET_MET200.Fill(dp_MPTMET,weight);
	  else HminDeltaPhiFail_MPTMET_MET200.Fill(dp_MPTMET,weight);
	}
	if (getMET() >250) {
	  HdeltaPhiMPTMET_MET250.Fill(dp_MPTMET,weight);
	  if (minDeltaPhi_j_MET > 0.3) HminDeltaPhiPass_MPTMET_MET250.Fill(dp_MPTMET,weight);
	  else HminDeltaPhiFail_MPTMET_MET250.Fill(dp_MPTMET,weight);
	}

	if (nbSSVM >=1)  HminDeltaPhiMETj_ge1b.Fill(minDeltaPhi_j_MET,weight);
	if (nbSSVM >=2)  HminDeltaPhiMETj_ge2b.Fill(minDeltaPhi_j_MET,weight);
	if (nbSSVM ==1)  HminDeltaPhiMETj_eq1b.Fill(minDeltaPhi_j_MET,weight);

	if (nbSSVM >=1)  HminDeltaPhiMETjAll_ge1b.Fill(minDeltaPhi_j_MET_All,weight);
	if (nbSSVM >=2)  HminDeltaPhiMETjAll_ge2b.Fill(minDeltaPhi_j_MET_All,weight);
	if (nbSSVM ==1)  HminDeltaPhiMETjAll_eq1b.Fill(minDeltaPhi_j_MET_All,weight);

	if (nbSSVM >=1)  HminDeltaPhiMETjAll30_ge1b.Fill(minDeltaPhi_j_MET_All30,weight);
	if (nbSSVM >=2)  HminDeltaPhiMETjAll30_ge2b.Fill(minDeltaPhi_j_MET_All30,weight);

	if (nbSSVM >=1)  HVminDeltaPhiMETj_ge1b.Fill(minDeltaPhi_j_MET,weight);
	if (nbSSVM >=2)  HVminDeltaPhiMETj_ge2b.Fill(minDeltaPhi_j_MET,weight);
	if (nbSSVM ==1)  HVminDeltaPhiMETj_eq1b.Fill(minDeltaPhi_j_MET,weight);

	if (nbSSVM >=1)  HV2minDeltaPhiMETj_ge1b.Fill(minDeltaPhi_j_MET,weight);
	if (nbSSVM >=2)  HV2minDeltaPhiMETj_ge2b.Fill(minDeltaPhi_j_MET,weight);
	if (nbSSVM ==1)  HV2minDeltaPhiMETj_eq1b.Fill(minDeltaPhi_j_MET,weight);
      }
      resetIgnoredCut();

      //without minDeltaPhi and without MET
      //this is for plotting the ratio of (pass minDeltaPhi) / (fail minDeltaPhi), etc
      setIgnoredCut("cutDeltaPhi");
      setIgnoredCut("cutMET");
      if (Cut(ientry) >= 0) {

	bool passMinDeltaPhi = minDeltaPhi_j_MET > 0.3;	//this could be passCut("cutDeltaPhi")
	//but then the code wouldn't be invariant on changes of the def'n of cutDeltaPhi...it just depends on what we want
	bool passMinDeltaPhiAll = minDeltaPhi_j_MET_All > 0.3;
	//	bool passMinDeltaPhiAll30 = minDeltaPhi_j_MET_All30 > 0.3; //not using this at the moment
	double largestJetRecoError = getLargestJetPtRecoError(3); //look at only the first 3 jets

	//same question here: use passCut("cutMET") or just apply a hard-coded limit?
	//opt for hard-coded
	bool passMET = getMET() > 150;

	if (passMinDeltaPhi) {
	  HminDeltaPhiPass.Fill(getMET(),weight);
	  if (nbSSVM >=1)  HminDeltaPhiPass_ge1b.Fill(getMET(),weight);
	  if (nbSSVM >=2)  HminDeltaPhiPass_ge2b.Fill(getMET(),weight);

	  if      (largestJetRecoError>50)  HminDeltaPhiPass_Over.Fill(getMET(), weight);
	  else if (largestJetRecoError<-50) HminDeltaPhiPass_Under.Fill(getMET(), weight);
	  else                              HminDeltaPhiPass_Close.Fill(getMET(), weight);
 	}
	else {
	  HminDeltaPhiFail.Fill(getMET(),weight);
	  if (nbSSVM >=1)  HminDeltaPhiFail_ge1b.Fill(getMET(),weight);
	  if (nbSSVM >=2)  HminDeltaPhiFail_ge2b.Fill(getMET(),weight);

	  if      (largestJetRecoError>50)  HminDeltaPhiFail_Over.Fill(getMET(), weight);
	  else if (largestJetRecoError<-50) HminDeltaPhiFail_Under.Fill(getMET(), weight);
	  else                              HminDeltaPhiFail_Close.Fill(getMET(), weight);
	}

	if (passMinDeltaPhiAll) {
	  HminDeltaPhiAllPass.Fill(getMET(),weight);
	  if (nbSSVM >=1)  HminDeltaPhiAllPass_ge1b.Fill(getMET(),weight);
	  if (nbSSVM >=2)  HminDeltaPhiAllPass_ge2b.Fill(getMET(),weight);
 	}
	else {
	  HminDeltaPhiAllFail.Fill(getMET(),weight);
	  if (nbSSVM >=1)  HminDeltaPhiAllFail_ge1b.Fill(getMET(),weight);
	  if (nbSSVM >=2)  HminDeltaPhiAllFail_ge2b.Fill(getMET(),weight);
	}

	//calculate scalar and vector sums of jet invisible energy (MC truth)
	double jetInvisibleEnergyHT = getJetInvisibleEnergyHT(); //scalar sum
	double jetInvisibleEnergyMHT = getJetInvisibleEnergyMHT(); //vector sum
	double maxdp = getMaxDeltaPhiMET30(99);
	if (passMinDeltaPhi && passMET) {//signal region; using All30 variant
	  HgenInvisibleEnergy_SR.Fill(jetInvisibleEnergyHT ,weight);
	  HgenInvisibleMissingEnergy_SR.Fill(jetInvisibleEnergyMHT ,weight);
	  HgenMET_SR.Fill(getGenMET() ,weight);
	  HmaxJetRecoError_SR.Fill(largestJetRecoError, weight);
	  HmaxDeltaPhiMETjAll_SR.Fill(maxdp, weight);
	}
	else if (!passMinDeltaPhi && passMET) { //region 'D'
	  HgenInvisibleEnergy_D.Fill(jetInvisibleEnergyHT ,weight);
	  HgenInvisibleMissingEnergy_D.Fill(jetInvisibleEnergyMHT ,weight);
	  HgenMET_D.Fill(getGenMET() ,weight);
	  HmaxJetRecoError_D.Fill(largestJetRecoError, weight);
	  HmaxDeltaPhiMETjAll_D.Fill(maxdp, weight);
	}
	else if (passMinDeltaPhi && !passMET) { // "low MET" region
	  HgenInvisibleEnergy_lowMET.Fill(jetInvisibleEnergyHT ,weight);
	  HgenInvisibleMissingEnergy_lowMET.Fill(jetInvisibleEnergyMHT ,weight);
	  HgenMET_lowMET.Fill(getGenMET() ,weight);
	  HmaxJetRecoError_lowMET.Fill(largestJetRecoError, weight);
	  HmaxDeltaPhiMETjAll_lowMET.Fill(maxdp, weight);
	}
	else if (!passMinDeltaPhi && !passMET) { // "lowlow" region
	  HgenInvisibleEnergy_lowlow.Fill(jetInvisibleEnergyHT ,weight);
	  HgenInvisibleMissingEnergy_lowlow.Fill(jetInvisibleEnergyMHT ,weight);
	  HgenMET_lowlow.Fill(getGenMET() ,weight);
	  HmaxJetRecoError_lowlow.Fill(largestJetRecoError, weight);
	  HmaxDeltaPhiMETjAll_lowlow.Fill(maxdp, weight);
	}
	//in fact the ABCD regions might exclude the very low MET events, but i'm not excluding them here
	else { cout<<"Found an event that did not fall into one of the ABCD regions!"<<endl; assert(0);}
      }
      resetIgnoredCut();
      
      //finally, apply all cuts (no b tagging)
      if (Cut(ientry) >= 0) {

	Hnbjets.Fill(nbSSVM, weight);
	HdeltaPhiMPTMET.Fill(dp_MPTMET,weight);
	H_METphi.Fill(getMETphi(),weight);
	Hjetpt1.Fill(leadJetPt,weight);
	Hjetphi1.Fill(leadJetPhi,weight);
	Hjeteta1.Fill(leadJetEta,weight);

	if (nbSSVM >=1)  {
	  H_METphi_ge1b.Fill(getMETphi(),weight);
	  HdeltaPhiMPTMET_ge1b.Fill(dp_MPTMET,weight);
	  Hjetpt1_ge1b.Fill(leadJetPt,weight);
	  Hjetphi1_ge1b.Fill(leadJetPhi,weight);
	  Hjeteta1_ge1b.Fill(leadJetEta,weight);

	  double bjetpt1=0,bjetphi1=-99,bjeteta1=-6;
	  int nbjetsfound=0;
	  //now rewriting this for kBaseline0
	  for (unsigned int ib = 0; ib< loosejetPhi->size(); ib++) {
	    if ( isGoodJet30(ib) && passSSVM(ib)) { //refind the b jets
	      nbjetsfound++;
	      if (nbjetsfound==1) { //if this is the *lead* b jet
		bjetpt1=getLooseJetPt(ib);
		bjetphi1=loosejetPhi->at(ib);
		bjeteta1=loosejetEta->at(ib);
		break; //only interested in the lead b jet at the moment
	      }
	    }
	  }
	  Hbjetpt1_ge1b.Fill(bjetpt1,weight);
	  Hbjetphi1_ge1b.Fill(bjetphi1,weight);
	  Hbjeteta1_ge1b.Fill(bjeteta1,weight);
	  
	  
	  if (nbSSVM ==1)  {
	    H_METphi_eq1b.Fill(getMETphi(),weight);
	    HdeltaPhiMPTMET_eq1b.Fill(dp_MPTMET,weight);
	    Hjetpt1_eq1b.Fill(leadJetPt,weight);
	    Hjetphi1_eq1b.Fill(leadJetPhi,weight);
	    Hjeteta1_eq1b.Fill(leadJetEta,weight);
	    
	    Hbjetpt1_eq1b.Fill(bjetpt1,weight);
	    Hbjetphi1_eq1b.Fill(bjetphi1,weight);
	    Hbjeteta1_eq1b.Fill(bjeteta1,weight);
	    
	  }
	  if (nbSSVM >=2)  {
	    H_METphi_ge2b.Fill(getMETphi(),weight);
	    HdeltaPhiMPTMET_ge2b.Fill(dp_MPTMET,weight);
	    Hjetpt1_ge2b.Fill(leadJetPt,weight);
	    Hjetphi1_ge2b.Fill(leadJetPhi,weight);
	    Hjeteta1_ge2b.Fill(leadJetEta,weight);
	    
	    Hbjetpt1_ge2b.Fill(bjetpt1,weight);
	    Hbjetphi1_ge2b.Fill(bjetphi1,weight);
	    Hbjeteta1_ge2b.Fill(bjeteta1,weight);
	    
	    deltaPhi_b1b2 = getDeltaPhib1b2(); //now updated for kBaseline0
	    HdeltaPhib1b2_ge2b.Fill(deltaPhi_b1b2,weight);
	  }
	}
      }
   }
   stopTimer(nentries);

   HminDeltaPhiRatio.Divide(&HminDeltaPhiPass,&HminDeltaPhiFail);
   HminDeltaPhiRatio_ge1b.Divide(&HminDeltaPhiPass_ge1b,&HminDeltaPhiFail_ge1b);
   HminDeltaPhiRatio_ge2b.Divide(&HminDeltaPhiPass_ge2b,&HminDeltaPhiFail_ge2b);

   HminDeltaPhiAllRatio.Divide(&HminDeltaPhiAllPass,&HminDeltaPhiAllFail);
   HminDeltaPhiAllRatio_ge1b.Divide(&HminDeltaPhiAllPass_ge1b,&HminDeltaPhiAllFail_ge1b);
   HminDeltaPhiAllRatio_ge2b.Divide(&HminDeltaPhiAllPass_ge2b,&HminDeltaPhiAllFail_ge2b);

   HminDeltaPhiRatio_Over.Divide(&HminDeltaPhiPass_Over,&HminDeltaPhiFail_Over);
   HminDeltaPhiRatio_Under.Divide(&HminDeltaPhiPass_Under,&HminDeltaPhiFail_Under);
   HminDeltaPhiRatio_Close.Divide(&HminDeltaPhiPass_Close,&HminDeltaPhiFail_Close);

   HminDeltaPhiRatio_MPTMET_MET150.Divide(&HminDeltaPhiPass_MPTMET_MET150,&HminDeltaPhiFail_MPTMET_MET150);
   HminDeltaPhiRatio_MPTMET_MET200.Divide(&HminDeltaPhiPass_MPTMET_MET200,&HminDeltaPhiFail_MPTMET_MET200);
   HminDeltaPhiRatio_MPTMET_MET250.Divide(&HminDeltaPhiPass_MPTMET_MET250,&HminDeltaPhiFail_MPTMET_MET250);

   fout.Write();
   fout.Close();

}


/* 
main plot-making loop

this code is quite old and may not be plotting the correct quantities in the Sync1 or Baseline0 schemes
[now largely fixed for Baseline0]

I'm not using this loop for making plots for the AN
*/
void basicLoop::Loop(unsigned int dataindex)
{
 
   if (fChain == 0) return;
   printState();

   Long64_t nentries = fChain->GetEntries(); //jmt: remove Fast


   TString inname=findInputName(); //uses fChain
   std::cout<<"Got an input file name as: "<<inname<<std::endl;
   
   double sigma = getCrossSection(inname);
   TString sampleName = getSampleName(inname);
   if (sigma<=0 && !isData_) return;

   if (nentries == 0) {std::cout<<"Chain has no entries!"<<std::endl; return;}

   double   weight = isData_ ? 1 : lumi * sigma / double(nentries); //calculate weight

   //open output file
   TString outfilename="plots."; 
   outfilename+=getCutDescriptionString();
   outfilename+=".";    outfilename+=sampleName; 
   if (isData_) {
     outfilename+="-";
     outfilename+=dataindex;
   }
   outfilename+=".root";
   TFile fout(outfilename,"RECREATE");

   //   TH1::SetDefaultSumw2(); //don't use the trick because of a root bug having to do with TH2

   // === make some histograms ===
   const   int offset=0;
   int njets_bins=10-offset;
   TH1D Hnjets("Hnjets","N of jets (RA2)",njets_bins,offset,njets_bins+offset);
   TH1D Hnjets_ge1b("Hnjets_ge1b","N of jets (RA2 && >=1b)",njets_bins,offset,njets_bins+offset);
   TH1D Hnjets_ge2b("Hnjets_ge2b","N of jets (RA2 && >=2b)",njets_bins,offset,njets_bins+offset);
   TH1D Hnjets_ge3b("Hnjets_ge3b","N of jets (RA2 && >=3b)",njets_bins,offset,njets_bins+offset);

   TH1D Hnbjets("Hnbjets","Number of SSVM jets",4,0,4);

   //DeltaPhi(MPT,MET)
   int nbins=50;

   TH1D HdeltaPhiMPTMET("HdeltaPhiMPTMET","DeltaPhi(MET,MPT) (RA2)",nbins,0,TMath::Pi());
   TH1D HdeltaPhiMPTMET_ge1b("HdeltaPhiMPTMET_ge1b","DeltaPhi(MET,MPT) (RA2 && >=1b)",nbins,0,TMath::Pi());
   TH1D HdeltaPhiMPTMET_ge2b("HdeltaPhiMPTMET_ge2b","DeltaPhi(MET,MPT) (RA2 && >=2b)",nbins,0,TMath::Pi());

   TH1D HminDeltaPhiMETj("HminDeltaPhiMETj","minDeltaPhi(j,MET) (RA2)",nbins,0,TMath::Pi());
   TH1D HminDeltaPhiMETb_ge1b("HminDeltaPhiMETb_ge1b","minDeltaPhi(b,MET) (RA2 && >=1b)",nbins,0,TMath::Pi());
   TH1D HminDeltaPhiMETj_ge1b("HminDeltaPhiMETj_ge1b","minDeltaPhi(j,MET) (RA2 && >=1b)",nbins,0,TMath::Pi());
   TH2D HdeltaPhib1b2_minDeltaPhiMETb("HdeltaPhib1b2_minDeltaPhiMETb","DeltaPhi(b1,b2) v minDeltaPhi(b,MET) (RA2 && >=2b)",nbins,0,TMath::Pi(),nbins,0,TMath::Pi());
   TH1D HdeltaPhib1b2("HdeltaPhib1b2","DeltaPhi(b1,b2) (>=2b)",nbins,0,TMath::Pi());

   double pt_min=30;
   double pt_max=500;
   TH1D Hjetpt1("Hjetpt1","pT of lead jet",nbins,pt_min,pt_max);
   TH1D Hjetpt1_ge1b("Hjetpt1_ge1b","pT of lead jet (>=1b)",nbins,pt_min,pt_max);
   TH1D Hjetpt1_ge2b("Hjetpt1_ge2b","pT of lead jet (>=2b)",nbins,pt_min,pt_max);

   TH1D Hjetphi1("Hjetphi1","phi of lead jet",nbins,-TMath::Pi(),TMath::Pi());
   TH1D Hjetphi1_ge1b("Hjetphi1_ge1b","phi of lead jet (>=1b)",nbins,-TMath::Pi(),TMath::Pi());
   TH1D Hjetphi1_ge2b("Hjetphi1_ge2b","phi of lead jet (>=2b)",nbins,-TMath::Pi(),TMath::Pi());

   double eta_max = 2.4;
   TH1D Hjeteta1("Hjeteta1","eta of lead jet",nbins,-eta_max,eta_max);
   TH1D Hjeteta1_ge1b("Hjeteta1_ge1b","eta of lead jet (>=1b)",nbins,-eta_max,eta_max);
   TH1D Hjeteta1_ge2b("Hjeteta1_ge2b","eta of lead jet (>=2b)",nbins,-eta_max,eta_max);

   TH1D Hbjetpt1_ge1b("Hbjetpt1_ge1b","pT of lead b jet (>=1b)",nbins,pt_min,pt_max);
   TH1D Hbjetpt1_ge2b("Hbjetpt1_ge2b","pT of lead b jet (>=2b)",nbins,pt_min,pt_max);

   TH1D Hbjetphi1_ge1b("Hbjetphi1_ge1b","phi of lead b jet (>=1b)",nbins,-TMath::Pi(),TMath::Pi());
   TH1D Hbjetphi1_ge2b("Hbjetphi1_ge2b","phi of lead b jet (>=2b)",nbins,-TMath::Pi(),TMath::Pi());

   TH1D Hbjeteta1_ge1b("Hbjeteta1_ge1b","eta of lead b jet (>=1b)",nbins,-eta_max,eta_max);
   TH1D Hbjeteta1_ge2b("Hbjeteta1_ge2b","eta of lead b jet (>=2b)",nbins,-eta_max,eta_max);


   int vnbins=8;
   double vbins[]={0, 0.15, 0.3, 0.5, 0.7, 1, 1.5, 2, TMath::Pi()};
   TH1D HVminDeltaPhiMETj("HVminDeltaPhiMETj","minDeltaPhi(j,MET) (RA2)",vnbins,vbins);
   TH1D HVminDeltaPhiMETj_ge1b("HVminDeltaPhiMETj_ge1b","minDeltaPhi(j,MET) (RA2 && >=1b)",vnbins,vbins);
   TH1D HVminDeltaPhiMETj_ge2b("HVminDeltaPhiMETj_ge2b","minDeltaPhi(j,MET) (RA2 && >=2b)",vnbins,vbins);

   TH1D HVminDeltaPhiMHTj("HVminDeltaPhiMHTj","minDeltaPhi(j,MHT) (RA2)",vnbins,vbins);
   TH1D HVminDeltaPhiMHTj_ge1b("HVminDeltaPhiMHTj_ge1b","minDeltaPhi(j,MHT) (RA2 && >=1b)",vnbins,vbins);
   TH1D HVminDeltaPhiMHTj_ge2b("HVminDeltaPhiMHTj_ge2b","minDeltaPhi(j,MHT) (RA2 && >=2b)",vnbins,vbins);

   int vnbins2=2;
   double vbins2[]={0,  0.3, TMath::Pi()};
   TH1D HV2minDeltaPhiMETj("HV2minDeltaPhiMETj","minDeltaPhi(j,MET) (RA2)",vnbins2,vbins2);
   TH1D HV2minDeltaPhiMETj_ge1b("HV2minDeltaPhiMETj_ge1b","minDeltaPhi(j,MET) (RA2 && >=1b)",vnbins2,vbins2);
   TH1D HV2minDeltaPhiMETj_ge2b("HV2minDeltaPhiMETj_ge2b","minDeltaPhi(j,MET) (RA2 && >=2b)",vnbins2,vbins2);

   TH1D HminDeltaPhiMETb_ge2b("HminDeltaPhiMETb_ge2b","minDeltaPhi(b,MET) (RA2 && >=2b)",nbins,0,TMath::Pi());
   TH1D HminDeltaPhiMETj_ge2b("HminDeltaPhiMETj_ge2b","minDeltaPhi(j,MET) (RA2 && >=2b)",nbins,0,TMath::Pi());

   double met_max=500,met_min=0;
   TH1D H_MHT("H_MHT","MHT",nbins, met_min , met_max);
   TH1D H_MET("H_MET","MET",nbins, met_min , met_max);
   TH1D H_METphi("H_METphi","MET phi",nbins, -TMath::Pi() , TMath::Pi());

   TH1D H_MHT_ge1b("H_MHT_ge1b","MHT (RA2 && >=1b)",nbins, met_min , met_max);
   TH1D H_MET_ge1b("H_MET_ge1b","MET (RA2 && >=1b)",nbins, met_min , met_max);
   TH1D H_METphi_ge1b("H_METphi_ge1b","MET phi (>=1b)",nbins, -TMath::Pi() , TMath::Pi());

   TH1D H_MHT_ge2b("H_MHT_ge2b","MHT (RA2 && >=2b)",nbins, met_min , met_max);
   TH1D H_MET_ge2b("H_MET_ge2b","MET (RA2 && >=2b)",nbins, met_min , met_max);
   TH1D H_METphi_ge2b("H_METphi_ge2b","MET phi (>=2b)",nbins, -TMath::Pi() , TMath::Pi());

   TH1D H_MHT_ge3b("H_MHT_ge3b","MHT (RA2 && >=3b)",nbins, met_min , met_max);
   TH1D H_MET_ge3b("H_MET_ge3b","MET (RA2 && >=3b)",nbins, met_min , met_max);
   TH1D H_METphi_ge3b("H_METphi_ge3b","MET phi (>=3b)",nbins, -TMath::Pi() , TMath::Pi());

   TH2D HdeltaPhib1b2_MET("HdeltaPhib1b2_MET","DeltaPhi(b1,b2) v MET (>=2b)",nbins,met_min,met_max,nbins,0,TMath::Pi());
   TH2D HminDeltaPhiMETb_MET_ge2b("HminDeltaPhiMETb_MET_ge2b","minDeltaPhi(b,MET) v MET (>=2b)",nbins,met_min,met_max,nbins,0,TMath::Pi());
   TH2D HminDeltaPhiMETj_MET_ge2b("HminDeltaPhiMETj_MET_ge2b","minDeltaPhi(j,MET) v MET (>=2b)",nbins,met_min,met_max,nbins,0,TMath::Pi());

   TH2D HminDeltaPhiMETb_HminDeltaPhiMETj_ge2b("HminDeltaPhiMETb_HminDeltaPhiMETj_ge2b","minDPhi(b,MET) v minDPhi(j,MET) (>=2b)",
					       nbins,0,TMath::Pi(),nbins,0,TMath::Pi());
   TH2D HdeltaPhiMPTMET_MET_ge2b("HdeltaPhiMPTMET_MET_ge2b","DeltaPhi(MET,MPT) v MET (RA2 && >=2b)",nbins,met_min,met_max,nbins,0,TMath::Pi());

   TH1D HtopDecayCategory_ge2b("HtopDecayCategory_ge2b","top decay category (RA2 && >=2b)",nTopCategories-1,0.5,nTopCategories-0.5);

   TH1D HdeltaPhi_bj_ge2b("HdeltaPhi_bj_ge2b","deltaPhi(b,j) (RA2 && >=2b)",nbins,0,TMath::Pi());

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

   TH1D H_HT("H_HT","HT",ht_nbins,ht_min,ht_max);
   TH1D H_HT_ge1b("H_HT_ge1b","HT (>=1b)",ht_nbins,ht_min,ht_max);
   TH1D H_HT_ge2b("H_HT_ge2b","HT (>=2b)",ht_nbins,ht_min,ht_max);

   //copy and paste the block of code with the TH1/2D definitions into afile and do
   // awk '/TH/ {print $2}' afile | awk -F\( '// {astr=sprintf("%s.Sumw2();", $1); print astr;}'
   Hnjets.Sumw2();
   Hnjets_ge1b.Sumw2();
   Hnjets_ge2b.Sumw2();
   Hnjets_ge3b.Sumw2();
   Hnbjets.Sumw2();
   HdeltaPhiMPTMET.Sumw2();
   HdeltaPhiMPTMET_ge1b.Sumw2();
   HdeltaPhiMPTMET_ge2b.Sumw2();
   HminDeltaPhiMETj.Sumw2();
   HminDeltaPhiMETb_ge1b.Sumw2();
   HminDeltaPhiMETj_ge1b.Sumw2();
   HdeltaPhib1b2_minDeltaPhiMETb.Sumw2();
   HdeltaPhib1b2.Sumw2();
   Hjetpt1.Sumw2();
   Hjetpt1_ge1b.Sumw2();
   Hjetpt1_ge2b.Sumw2();
   Hjetphi1.Sumw2();
   Hjetphi1_ge1b.Sumw2();
   Hjetphi1_ge2b.Sumw2();
   Hjeteta1.Sumw2();
   Hjeteta1_ge1b.Sumw2();
   Hjeteta1_ge2b.Sumw2();
   Hbjetpt1_ge1b.Sumw2();
   Hbjetpt1_ge2b.Sumw2();
   Hbjetphi1_ge1b.Sumw2();
   Hbjetphi1_ge2b.Sumw2();
   Hbjeteta1_ge1b.Sumw2();
   Hbjeteta1_ge2b.Sumw2();
   HVminDeltaPhiMETj.Sumw2();
   HVminDeltaPhiMETj_ge1b.Sumw2();
   HVminDeltaPhiMETj_ge2b.Sumw2();
   HVminDeltaPhiMHTj.Sumw2();
   HVminDeltaPhiMHTj_ge1b.Sumw2();
   HVminDeltaPhiMHTj_ge2b.Sumw2();
   HV2minDeltaPhiMETj.Sumw2();
   HV2minDeltaPhiMETj_ge1b.Sumw2();
   HV2minDeltaPhiMETj_ge2b.Sumw2();
   HminDeltaPhiMETb_ge2b.Sumw2();
   HminDeltaPhiMETj_ge2b.Sumw2();
   H_MHT.Sumw2();
   H_MET.Sumw2();
   H_METphi.Sumw2();
   H_MHT_ge1b.Sumw2();
   H_MET_ge1b.Sumw2();
   H_METphi_ge1b.Sumw2();
   H_MHT_ge2b.Sumw2();
   H_MET_ge2b.Sumw2();
   H_METphi_ge2b.Sumw2();
   H_MHT_ge3b.Sumw2();
   H_MET_ge3b.Sumw2();
   H_METphi_ge3b.Sumw2();
   HdeltaPhib1b2_MET.Sumw2();
   HminDeltaPhiMETb_MET_ge2b.Sumw2();
   HminDeltaPhiMETj_MET_ge2b.Sumw2();
   HminDeltaPhiMETb_HminDeltaPhiMETj_ge2b.Sumw2();
   HdeltaPhiMPTMET_MET_ge2b.Sumw2();
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
   H_HT.Sumw2();
   H_HT_ge1b.Sumw2();
   H_HT_ge2b.Sumw2();


   //keep track of performance
   TDatime starttime; //default ctor is for current time

   //event loop
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries ;jentry++) {
     if (jentry%1000000==0) cout<<int(100*double(jentry)/double(nentries))<<endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) {assert(0);}
      nb = GetEntry(jentry);   nbytes += nb; //use member function GetEntry instead of fChain->

      int njets = nGoodJets(); //update to use baseline0

      MHT = getMHT(); //here we overwrite an ntuple variable!
      HT = getHT(); //here we overwrite an ntuple variable!

      if (Cut(ientry) < 0) continue; //jmt use cut

      Hnbjets.Fill(nbSSVM,weight);

      //calculate things
      double dp_MPTMET = getDeltaPhiMPTMET();

      Hnjets.Fill( njets, weight );
      HdeltaPhiMPTMET.Fill( dp_MPTMET,weight);
      H_MHT.Fill(MHT ,weight);
      H_MET.Fill( getMET() ,weight);
      H_METphi.Fill( getMETphi() ,weight);

      double leadJetPt= jetPtOfN(1);
      double leadJetPhi= jetPhiOfN(1);
      double leadJetEta= jetEtaOfN(1);
      Hjetpt1.Fill(leadJetPt,weight);
      Hjetphi1.Fill(leadJetPhi,weight);
      Hjeteta1.Fill(leadJetEta,weight);

      H_HT.Fill(HT,weight);

      double minDeltaPhi_j_MET= getMinDeltaPhiMET(3);

      HminDeltaPhiMETj.Fill(minDeltaPhi_j_MET,weight);
      HVminDeltaPhiMETj.Fill(minDeltaPhi_j_MET,weight);
      HV2minDeltaPhiMETj.Fill(minDeltaPhi_j_MET,weight);

      if ( nbSSVM < 1) continue; //cut on the number of b tags
      Hnjets_ge1b.Fill( njets, weight );
      H_MHT_ge1b.Fill(MHT ,weight);
      H_MET_ge1b.Fill(getMET() ,weight);
      H_METphi_ge1b.Fill( getMETphi() ,weight);

      Hjetpt1_ge1b.Fill( leadJetPt,weight);
      Hjetphi1_ge1b.Fill(leadJetPhi,weight);
      Hjeteta1_ge1b.Fill(leadJetEta,weight);

      HdeltaPhiMPTMET_ge1b.Fill( dp_MPTMET,weight);

      H_HT_ge1b.Fill(HT,weight);

      double minDeltaPhi_b_MET= getMinDeltaPhibMET(); //this has not been updated for full kBaseline0 compatibility

      HminDeltaPhiMETb_ge1b.Fill(minDeltaPhi_b_MET,weight);
      HminDeltaPhiMETj_ge1b.Fill(minDeltaPhi_j_MET,weight);
      HVminDeltaPhiMETj_ge1b.Fill(minDeltaPhi_j_MET,weight);
      HV2minDeltaPhiMETj_ge1b.Fill(minDeltaPhi_j_MET,weight);

      int topcat = isData_ ? -99 : getTopDecayCategory(); //no sense in looking at MC truth in the data
      double minDeltaR_bj=999;
      int nbjetsfound=0;
      double bjetpt1=0,bjetphi1=-99,bjeteta1=-6;
      //now rewriting this for kBaseline0
      for (unsigned int ib = 0; ib< loosejetPhi->size(); ib++) {
	if ( isGoodJet30(ib) && passSSVM(ib)) { //refind the b jets
	  nbjetsfound++;
	  if (nbjetsfound==1) { //if this is the *lead* b jet
	    bjetpt1=getLooseJetPt(ib);
	    bjetphi1=loosejetPhi->at(ib);
	    bjeteta1=loosejetEta->at(ib);
	  }
	  //	  double deltaPhi_bj=getMinDeltaPhi_bj(ib);
	  //	  if ( nbSSVM < 2) {	  HdeltaPhi_bj_ge2b.Fill(deltaPhi_bj,weight);}

	  //	  double mdr=getMinDeltaR_bj(ib);
	  //	  if ( nbSSVM < 2) {
	  //	    HdeltaR_bj_ge2b.Fill(mdr,weight);
	  //	    HdeltaR_bj_vTopCat_ge2b.Fill(topcat,mdr,weight);
	  //	  }
	  //	  if (mdr<minDeltaR_bj) minDeltaR_bj=mdr;
	}
      }
      Hbjetpt1_ge1b.Fill(bjetpt1,weight);

      //      HminDeltaR_bj_vMET_ge1b.Fill(getMET(), minDeltaR_bj,weight);
      //      HminDeltaR_bj_vMET_ABCD_ge1b.Fill(getMET(), minDeltaR_bj,weight);

      if ( nbSSVM < 2) continue; //cut on the number of b tags
      Hbjetpt1_ge2b.Fill(bjetpt1,weight);
      H_METphi_ge2b.Fill( getMETphi() ,weight);

      Hnjets_ge2b.Fill( njets, weight );
      H_MHT_ge2b.Fill(MHT ,weight);
      H_MET_ge2b.Fill(getMET() ,weight);
      HdeltaPhiMPTMET_ge2b.Fill( dp_MPTMET,weight);

      Hjetpt1_ge2b.Fill(leadJetPt,weight);
      Hjetphi1_ge2b.Fill(leadJetPhi,weight);
      Hjeteta1_ge2b.Fill(leadJetEta,weight);

      H_HT_ge2b.Fill(HT,weight);

      HminDeltaPhiMETb_ge2b.Fill(minDeltaPhi_b_MET,weight);
      HminDeltaPhiMETj_ge2b.Fill(minDeltaPhi_j_MET,weight);
      HVminDeltaPhiMETj_ge2b.Fill(minDeltaPhi_j_MET,weight);
      HV2minDeltaPhiMETj_ge2b.Fill(minDeltaPhi_j_MET,weight);

      double deltaPhi_b1b2 = getDeltaPhib1b2(); //now! updated for kBaseline0
      HdeltaPhib1b2.Fill(deltaPhi_b1b2,weight);
      HdeltaPhib1b2_minDeltaPhiMETb.Fill(minDeltaPhi_b_MET,deltaPhi_b1b2,weight);
      HdeltaPhib1b2_MET.Fill(getMET(),deltaPhi_b1b2,weight);

      HminDeltaPhiMETb_MET_ge2b.Fill(getMET(),minDeltaPhi_b_MET,weight);
      HminDeltaPhiMETj_MET_ge2b.Fill(getMET(),minDeltaPhi_j_MET,weight);

      HminDeltaPhiMETb_HminDeltaPhiMETj_ge2b.Fill(minDeltaPhi_j_MET,minDeltaPhi_b_MET,weight);

      HdeltaPhiMPTMET_MET_ge2b.Fill(getMET(),dp_MPTMET,weight);

      HtopDecayCategory_ge2b.Fill(topcat,weight);

      //make a plot that has only one entry per event
      HminDeltaR_bj_ge2b.Fill(minDeltaR_bj,weight);
      HminDeltaR_bj_vTopCat_ge2b.Fill(topcat,minDeltaR_bj,weight);

      HminDeltaR_bj_vMET_ge2b.Fill(getMET(), minDeltaR_bj,weight);
      HminDeltaR_bj_vMET_ABCD_ge2b.Fill(getMET(), minDeltaR_bj,weight);

      if ( nbSSVM < 3) continue; //cut on the number of b tags
      Hnjets_ge3b.Fill( njets, weight );
      H_MHT_ge3b.Fill(MHT ,weight);
      H_MET_ge3b.Fill(getMET() ,weight);
      H_METphi_ge3b.Fill( getMETphi() ,weight);
   }
   TDatime stoptime; //default ctor is for current time
   UInt_t elapsed= stoptime.Convert() - starttime.Convert();

   cout<<"events / time = "<<nentries<<" / "<<elapsed<<" = "<<double(nentries)/double(elapsed)<<" Hz"<<endl;
   fout.Write();
   fout.Close();

}


void basicLoop::nbLoop()
{

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
      nb = GetEntry(jentry);   nbytes += nb; //use member function GetEntry instead of fChain->

      //fill before cut histograms
      HSUSY_nb.Fill(SUSY_nb);

      if (Cut(ientry) < 0) continue; //jmt use cut
      //fill after cut histograms
      HSUSY_nb_RA2.Fill(SUSY_nb);
   }

  fout.Write();
  fout.Close();

}

void basicLoop::screendump()
{

  //for now hard-coded with a couple events
  /* some data events
  specifyEvent(143962, 2, 732462);
  specifyEvent(143962, 2, 810194);
  */
  // specifyEvent(148862, 75, 120899194); //event that fails my cuts but not Don

  /* LM0 events */
  /*
    specifyEvent(1, 6, 817);
    specifyEvent(1, 186, 25463);
  */

  //Fall10 ttbar events
  specifyEvent(1,1,122488);

  //LM9 events
  //  specifyEvent(1,11,4313);
  //  specifyEvent(1,12,4816);

  //  specifyEvent( 1, 325, 131361);
  //  specifyEvent( 1, 1, 45);

  ULong64_t nfound=0;

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries(); //jmt: remove Fast
   const TString sp=" ";
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = GetEntry(jentry);   nbytes += nb; //use member function GetEntry instead of fChain->

      if (true) {
	if (eventIsSpecified() ) {
	  nfound++;
	  cout<<"--- record for event: ("<<jentry <<") run,ls,ev = "<<runNumber<<sp<<lumiSection<<sp<<eventNumber<<endl;
	  //Show() doesn't cut it-- it just shows the pointer addresses!
	  
	  if (false) {
	    for (unsigned int i=0 ; i<cutTags_.size(); i++) {
	      if (cutRequired(cutTags_[i])) {
		TString passstr = passCut(cutTags_[i]) ? "pass" : "fail";
		cout<< cutNames_[cutTags_[i]] <<sp<<passstr<<endl;
	      }
	    }
	  }
	  
	  if (false) {
	    cout<<" jet info (pT, Eta, Jet ID, isGood, SSVHE) n good jets = "<<nGoodJets()<<endl;
	    for (unsigned int ijet=0; ijet<loosejetPt->size(); ijet++) {
	      TString jetisgood = isGoodJet(ijet) ? "Good" : "notGood";
	      cout<<"\tjet "<<ijet<<": "<<getLooseJetPt(ijet) <<sp<<loosejetEta->at(ijet)<<sp<<loosejetPassLooseID->at(ijet)<<sp<<jetisgood
		  <<sp<< loosejetBTagDisc_simpleSecondaryVertexHighEffBJetTags->at(ijet)<<endl;
	    }
	    cout<<"\tEvent HT = "<<getHT()<<endl;
	  }

	  if (true) {
	    //	    cout<<" jet info (pT, Eta, Jet ID, isGood, SSVHE) n good jets = "<<nGoodJets()<<endl;
	    for (unsigned int ijet=0; ijet<loosejetPt->size(); ijet++) {
	      TString jetisgood = isGoodJet(ijet) ? "Good" : "notGood";
	      cout<<"jet "<<ijet<<endl<<"pt = "<<loosejetPt->at(ijet) <<", eta = "<<loosejetEta->at(ijet)<<endl
		  <<"variable jet unc = "<<loosejetJECUncPlus->at(ijet)<<endl
		  <<"total JEC uncertainty factor = "<<1+sqrt( loosejetJECUncPlus->at(ijet) * loosejetJECUncPlus->at(ijet) +0.053 * 0.053)<<endl
		  <<"factor * pT = "<< (1+sqrt( loosejetJECUncPlus->at(ijet) * loosejetJECUncPlus->at(ijet) +0.053 * 0.053) )* loosejetPt->at(ijet) <<" compare to: "<<getLooseJetPt(ijet)<<endl;
	    }
	  }

	}
	if (nfound == specifiedEvents_.size()) break; //save time at the end
      }

      //in case we just want to dump every event number
      if (false) {
	if (runNumber==143962) {
	  cout<<"run ls event "<<runNumber<<" "<<lumiSection<<" "<<eventNumber<<endl;
	}
      }
      
   }
}

void basicLoop::triggerPlot()
{

  /*
    in some ways it is stupid to do this independently from ::Loop(). oh well.....
  */

   if (fChain == 0) return;
   printState();

   //open output file
   TString outfilename="triggerCurve."; 
   outfilename+=getCutDescriptionString();
   outfilename+=".";    outfilename+=getSampleName(findInputName());; 
   outfilename+=".root";
   TFile fout(outfilename,"RECREATE");

   double ht_min=0,ht_max=900;
   int ht_nbins=300;
   TH1D H_HT_AfterCuts("H_HT_AfterCuts","ht after cuts",ht_nbins,ht_min,ht_max);
   TH1D H_HT_AfterTrigger("H_HT_AfterTrigger","ht after cuts+trigger",ht_nbins,ht_min,ht_max);
   TH1D H_HT_AfterHTU("H_HT_AfterHTU","ht after cuts+fake trigger",ht_nbins,ht_min,ht_max);

   TH1D H_HTU_BeforeTrigger("H_HTU_BeforeTrigger","uncorrected offline ht before trigger cut",ht_nbins,ht_min,ht_max);
   TH1D H_HTU_AfterTrigger("H_HTU_AfterTrigger","uncorrected offline ht after trigger cut",ht_nbins,ht_min,ht_max);

   H_HT_AfterCuts.Sumw2();
   H_HT_AfterTrigger.Sumw2();

   H_HTU_BeforeTrigger.Sumw2();
   H_HTU_AfterTrigger.Sumw2();

   Long64_t nentries = fChain->GetEntries(); //jmt: remove Fast

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries ; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = GetEntry(jentry);   nbytes += nb; //use member function GetEntry instead of fChain->

      if (Cut(ientry) < 0) continue; //jmt use cut

      double ht = getHT_Sync1(); //recalculate HT
      double htu = getUncorrectedHT(30); //does the HT100U in MC use 20 or 30???

      H_HT_AfterCuts.Fill(ht);
      H_HTU_BeforeTrigger.Fill(htu);
      if (  passHLT("HLT_HT100U") ) {
	H_HT_AfterTrigger.Fill(ht);
	H_HTU_AfterTrigger.Fill(htu);
      }
      if ( htu > 100 ) {
	H_HT_AfterHTU.Fill(ht);
      }

      //debug
      if (true && jentry<25) {
	cout<<passHLT("HLT_HT100U")<<"\t"<<htu<<"\t"<<ht<<endl;
      }

   }

   TH1D Heff("Heff","pass trig fraction",ht_nbins,ht_min,ht_max);
   TH1D Heff_FakeHTUTrigger("Heff_FakeHTUTrigger","pass fraction using HTU (not trigger)",ht_nbins,ht_min,ht_max);
   TH1D H_HTU_turnon("H_HTU_turnon","pass trig fraction",ht_nbins,ht_min,ht_max);
   H_HTU_turnon.SetXTitle("Uncorrected offline HT");
   Heff.Sumw2();
   Heff_FakeHTUTrigger.Sumw2();
   H_HTU_turnon.Sumw2();
   Heff.Divide(&H_HT_AfterTrigger,&H_HT_AfterCuts,1,1,"B");
   Heff_FakeHTUTrigger.Divide(&H_HT_AfterTrigger,&H_HT_AfterHTU,1,1,"B");
   H_HTU_turnon.Divide(&H_HTU_AfterTrigger,&H_HTU_BeforeTrigger,1,1,"B");

   fout.Write();
   fout.Close();

}



void basicLoop::triggerPlotData()
{

  /*
    note...this code isn't going to work until I get HLT_HT100U_v3 into the ntuple
    this means rerunning over the MultiJet datasets
  */

   if (fChain == 0) return;
   printState();

   //open output file
   TString outfilename="triggerData."; 
   outfilename+=getCutDescriptionString();
   outfilename+=".";    outfilename+=getSampleName(findInputName());; 
   outfilename+=".root";
   TFile fout(outfilename,"RECREATE");

   double ht_min=0,ht_max=900;
   int ht_nbins=300;
   TH1D H_HT_AfterHT100U("H_HT_AfterHT100U","ht after prescaled trigger",ht_nbins,ht_min,ht_max);
   TH1D H_HT_AfterHT150U("H_HT_AfterHT150U","ht after physics trigger",ht_nbins,ht_min,ht_max);

   /*
   TH1D H_HT_AfterHTU("H_HT_AfterHTU","ht after cuts+fake trigger",ht_nbins,ht_min,ht_max);

   TH1D H_HTU_BeforeTrigger("H_HTU_BeforeTrigger","uncorrected offline ht before trigger cut",ht_nbins,ht_min,ht_max);
   TH1D H_HTU_AfterTrigger("H_HTU_AfterTrigger","uncorrected offline ht after trigger cut",ht_nbins,ht_min,ht_max);
   H_HTU_BeforeTrigger.Sumw2();
   H_HTU_AfterTrigger.Sumw2();
   */

   H_HT_AfterHT100U.Sumw2();
   H_HT_AfterHT150U.Sumw2();

   Long64_t nentries = fChain->GetEntries(); //jmt: remove Fast

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries ; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = GetEntry(jentry);   nbytes += nb; //use member function GetEntry instead of fChain->

      //look only at very late runs
      if (runNumber < 148783) continue;

      //now verify that these triggers existed in this menu
      if (getHLTPrescale("HLT_HT100U_v3") <1 || getHLTPrescale("HLT_HT150U_v3")<1) {
	cout<<"Problem with one of the triggers!"<<endl;
	cout<<getHLTPrescale("HLT_HT100U_v3")<<endl;
	cout<<getHLTPrescale("HLT_HT150U_v3")<<endl;
	return;
      }

      if (Cut(ientry) < 0) continue; //jmt use cut

      if ( !passHLT("HLT_HT100U_v3") ) continue;

      double ht = getHT_Sync1(); //recalculate HT

      H_HT_AfterHT100U.Fill(ht);

      if (  passHLT("HLT_HT150U_v3") ) {
	H_HT_AfterHT150U.Fill(ht);
      }

   }

   TH1D Heff("Heff","pass trig fraction",ht_nbins,ht_min,ht_max);
   Heff.Sumw2();

   Heff.Divide(&H_HT_AfterHT150U,&H_HT_AfterHT100U,1,1,"B");

   fout.Write();
   fout.Close();

}

//note that this function is now crippled because I don't fill lastTriggerPass_
void basicLoop::triggerTest()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries(); //jmt: remove Fast

   std::map<TString, unsigned int> npassHLT;

   std::map<TString, unsigned int> runNumber_min;
   std::map<TString, unsigned int> runNumber_max;
   std::map<TString, bool> printed;
   std::vector<TString> triggers;
   triggers.push_back("HLT_HT100U");
   triggers.push_back("HLT_HT120U");
   triggers.push_back("HLT_HT140U");
   triggers.push_back("HLT_HT150U_v3");
   for (unsigned int itrg = 0; itrg < triggers.size() ; itrg++) {
     TString thistrg = triggers.at(itrg);
     runNumber_min[thistrg] = 9999999;
     runNumber_max[thistrg] = 0;
     printed[thistrg]=false;
   }
   
   //open output file
   TString outfilename="trigger."; 
   outfilename+=getCutDescriptionString();
   outfilename+=".";    outfilename+=getSampleName(findInputName());; 
   /*  i think i don't need this as long as i make a big tchain of all data files
   if (isData_) {
     outfilename+="-";
     outfilename+=dataindex;
   }
   */
   outfilename+=".root";
   TFile fout(outfilename,"RECREATE");


   int minrun = 125000;
   int maxrun = 150000;
   TH1D Hrunnum("Hrunnum","run number",maxrun-minrun,minrun,maxrun);
   TH1D Hrunnum_passHLT("Hrunnum_passHLT","run number",maxrun-minrun,minrun,maxrun);

   Hrunnum.Sumw2();
   Hrunnum_passHLT.Sumw2();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = GetEntry(jentry);   nbytes += nb; //use member function GetEntry instead of fChain->

      Hrunnum.Fill(runNumber);

      //note that I have killed the Cut() function here and just call the trigger cut explicitly
      if (passHLT() ) {
	Hrunnum_passHLT.Fill(runNumber);
	if (npassHLT.find(lastTriggerPass_) == npassHLT.end()) {
	  npassHLT[lastTriggerPass_]=1;
	  cout<<"An event passed: "<<lastTriggerPass_<<endl;
	}
	else {  npassHLT[lastTriggerPass_] =  npassHLT[lastTriggerPass_] +1;}
      }

      //after this function was written, i implemented a better way to do this
      //leave the old code alone, but add some new code here:
      for (unsigned int itrg = 0; itrg < triggers.size() ; itrg++) {
	TString thistrg=triggers.at(itrg);
	int prescale= getHLTPrescale(thistrg);
	if (runNumber ==148032 && !printed[thistrg]) {
	  cout<<runNumber<<"\t"<<thistrg<<"\t"<<prescale<<endl;
	  printed[thistrg]=true;
	}
	if (prescale == 1) {
	  if ( runNumber < runNumber_min[thistrg] ) runNumber_min[thistrg] = runNumber;
	  if ( runNumber > runNumber_max[thistrg] ) runNumber_max[thistrg] = runNumber;
	}
      }

   }
   
   TH1D Hrunnum_ratio("Hrunnum_ratio","pass trig fraction",maxrun-minrun,minrun,maxrun);
   Hrunnum_ratio.Sumw2();
   Hrunnum_ratio.Divide(&Hrunnum_passHLT,&Hrunnum);

   for ( std::map<TString, unsigned int>::const_iterator ihlt = npassHLT.begin() ; ihlt!=npassHLT.end(); ++ihlt) {
     cout<< ihlt->first <<"\t"<<ihlt->second<<endl;
   }

   fout.Write();
   fout.Close();
   cout<<"\n======================="<<endl;
   for (unsigned int itrg = 0; itrg < triggers.size() ; itrg++) {
     TString thistrg=triggers.at(itrg);
     cout<<thistrg<<"\t"<<runNumber_min[thistrg]<<"\t"<<runNumber_max[thistrg]<<endl;
   }

}

