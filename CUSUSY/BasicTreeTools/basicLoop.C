#define basicLoop_cxx
#include "basicLoop.h"
#include <TH1.h>
#include <TH2.h>
#include <TDatime.h>
#include <TFile.h>
#include <TStyle.h>
#include <fstream>

#include "TMath.h"

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
      nb = GetEntry(jentry);   nbytes += nb; //use member function GetEntry instead of fChain->

      if (Cut(ientry) < 0) continue; //jmt use cut
   }
}

/*
print a cut flow table
lumi is set in basicLoop.h
*/
void basicLoop::cutflow()
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

  for (unsigned int i=0 ; i<cutTags_.size(); i++) {
    npass.push_back(0);
  }
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = GetEntry(jentry);   nbytes += nb; //use member function GetEntry instead of fChain->

    //cout<<"== "<<jentry<<endl;
  
    //i hate to do this, but I'm going to put a hack here to check the run number
    //    if (runNumber != 143962) continue;

    for (unsigned int i=0 ; i<cutTags_.size(); i++) {
      //cout<<i<<endl;
      if (cutRequired(cutTags_[i]) && passCut(cutTags_[i]) )   npass.at(i) = npass.at(i) +1;
      else if (cutRequired(cutTags_[i]) && !passCut(cutTags_[i]) ) break;

      //optional code to dump events to file
      /*
      if (cutTags_[i] == "cutJetPt1") {
	cout<<"run lumi event = "<<runNumber<<" "<<lumiSection<<" "<<eventNumber<<endl;
      }
      */
    }
    
  }

  TString samplename=  getSampleName(findInputName());
  TString outfilename="cutflow."; 
  outfilename+=getCutDescriptionString();
  outfilename+=".";    outfilename+=samplename; 
  outfilename+=".dat";
  ofstream file(outfilename.Data());
  
  for (unsigned int i=0 ; i<npass.size(); i++) {
    
    if (cutRequired(cutTags_[i])) {
      
      //error on n is sqrt n
      double error = sqrt(npass.at(i));
      double weighted = npass.at(i) * weight;
      double weighted_error = error*weight;
      
      char ccc[150];
      sprintf(ccc,"%20s %15d | %.2f | Weighted = %f +/- %f",cutNames_[cutTags_[i]].Data(),npass.at(i),100*double(npass.at(i))/double(npass.at(0)),weighted,weighted_error);
      cout<<ccc<<endl;

      file <<  weighted<<"\t" << weighted_error<<endl;
    }
  }
  
  file.close();
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
  TString outfilename="/cu1/joshmt/ABCDtrees/ABCDtree.";
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
  double myMET;
  double myMHT;
  double minDeltaPhiMET;
  double minDeltaPhiMHT;
  double minDeltaRbj;
  double DeltaPhiMPTMET;

  TTree ABCDtree("ABCDtree","ABCD tree");
  ABCDtree.Branch("weight",&weight,"weight/D");
  ABCDtree.Branch("MET",&myMET,"MET/D");
  ABCDtree.Branch("MHT",&myMHT,"MHT/D");
  ABCDtree.Branch("minDeltaPhiMET",&minDeltaPhiMET,"minDeltaPhiMET/D");
  ABCDtree.Branch("minDeltaPhiMHT",&minDeltaPhiMHT,"minDeltaPhiMHT/D");
  ABCDtree.Branch("minDeltaRbj",&minDeltaRbj,"minDeltaRbj/D");
  ABCDtree.Branch("DeltaPhiMPTMET",&DeltaPhiMPTMET,"DeltaPhiMPTMET/D");
  ABCDtree.Branch("nbSSVM", &nbSSVM, "nbSSVM/I");

  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = GetEntry(jentry);   nbytes += nb; //use member function GetEntry instead of fChain->

    if (Cut(ientry) < 0) continue; //jmt use cut
    myMET = getMET();
    myMHT = MHT;
    minDeltaPhiMET = getMinDeltaPhiMET(3);
    minDeltaPhiMHT = getMinDeltaPhiMHT(3);
    minDeltaRbj = getOverallMinDeltaR_bj();
    DeltaPhiMPTMET = getDeltaPhiMPTMET();

    ABCDtree.Fill(); 
  }

  fout.Write();
  fout.Close();
  
}

/*
experimental new loop -- for use with Baseline0 cut scheme
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
     }
   }
      

   //keep track of performance
   TDatime starttime; //default ctor is for current time
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
	}
	else if (cutRequired(cutTags_[i]) && !passCut(cutTags_[i]) ) break;
      }
     
   }
   TDatime stoptime; //default ctor is for current time
   UInt_t elapsed= stoptime.Convert() - starttime.Convert();

   cout<<"events / time = "<<nentries<<" / "<<elapsed<<" = "<<double(nentries)/double(elapsed)<<" Hz"<<endl;
 

   fout.Write();
   fout.Close();

}

/* 
main plot-making loop

this code is quite old and may not be plotting the correct quantities in the Sync1 or Baseline0 schemes
*/
void basicLoop::Loop(unsigned int dataindex)
{
  const double pi=4*atan(1);

  //
 
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
   TH1D HdeltaPhiMPTMET_ge1b("HdeltaPhiMPTMET_ge1b","DeltaPhi(MET,MPT) (RA2 && >=1b)",nbins,0,pi);
   TH1D HdeltaPhiMPTMET_ge2b("HdeltaPhiMPTMET_ge2b","DeltaPhi(MET,MPT) (RA2 && >=2b)",nbins,0,pi);

   TH1D HminDeltaPhiMETj("HminDeltaPhiMETj","minDeltaPhi(j,MET) (RA2)",nbins,0,pi);
   TH1D HminDeltaPhiMETb_ge1b("HminDeltaPhiMETb_ge1b","minDeltaPhi(b,MET) (RA2 && >=1b)",nbins,0,pi);
   TH1D HminDeltaPhiMETj_ge1b("HminDeltaPhiMETj_ge1b","minDeltaPhi(j,MET) (RA2 && >=1b)",nbins,0,pi);
   TH2D HdeltaPhib1b2_minDeltaPhiMETb("HdeltaPhib1b2_minDeltaPhiMETb","DeltaPhi(b1,b2) v minDeltaPhi(b,MET) (RA2 && >=2b)",nbins,0,pi,nbins,0,pi);

   double pt_min=50;
   double pt_max=400;
   TH1D Hjetpt1("Hjetpt1","pT of lead jet",nbins,pt_min,pt_max);
   TH1D Hjetpt1_ge1b("Hjetpt1_ge1b","pT of lead jet (>=1b)",nbins,pt_min,pt_max);
   TH1D Hjetpt1_ge2b("Hjetpt1_ge2b","pT of lead jet (>=2b)",nbins,pt_min,pt_max);

   TH1D Hbjetpt1_ge1b("Hbjetpt1_ge1b","pT of lead b jet (>=1b)",nbins,pt_min,pt_max);
   TH1D Hbjetpt1_ge2b("Hbjetpt1_ge2b","pT of lead b jet (>=2b)",nbins,pt_min,pt_max);


   int vnbins=8;
   double vbins[]={0, 0.15, 0.3, 0.5, 0.7, 1, 1.5, 2, pi};
   TH1D HVminDeltaPhiMETj("HVminDeltaPhiMETj","minDeltaPhi(j,MET) (RA2)",vnbins,vbins);
   TH1D HVminDeltaPhiMETj_ge1b("HVminDeltaPhiMETj_ge1b","minDeltaPhi(j,MET) (RA2 && >=1b)",vnbins,vbins);
   TH1D HVminDeltaPhiMETj_ge2b("HVminDeltaPhiMETj_ge2b","minDeltaPhi(j,MET) (RA2 && >=2b)",vnbins,vbins);

   TH1D HVminDeltaPhiMHTj("HVminDeltaPhiMHTj","minDeltaPhi(j,MHT) (RA2)",vnbins,vbins);
   TH1D HVminDeltaPhiMHTj_ge1b("HVminDeltaPhiMHTj_ge1b","minDeltaPhi(j,MHT) (RA2 && >=1b)",vnbins,vbins);
   TH1D HVminDeltaPhiMHTj_ge2b("HVminDeltaPhiMHTj_ge2b","minDeltaPhi(j,MHT) (RA2 && >=2b)",vnbins,vbins);

   int vnbins2=2;
   double vbins2[]={0,  0.3, pi};
   TH1D HV2minDeltaPhiMETj("HV2minDeltaPhiMETj","minDeltaPhi(j,MET) (RA2)",vnbins2,vbins2);
   TH1D HV2minDeltaPhiMETj_ge1b("HV2minDeltaPhiMETj_ge1b","minDeltaPhi(j,MET) (RA2 && >=1b)",vnbins2,vbins2);
   TH1D HV2minDeltaPhiMETj_ge2b("HV2minDeltaPhiMETj_ge2b","minDeltaPhi(j,MET) (RA2 && >=2b)",vnbins2,vbins2);

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

//    TH2D HpassMET45("HpassMET45","pass MET45 versus MET (RA2)",nbins,met_min,met_max,2,0,2);
//    TH2D HpassMET45_ge1b("HpassMET45_ge1b","pass MET45 versus MET (RA2 && >=1b)",nbins,met_min,met_max,2,0,2);
//    TH2D HpassMET45_ge2b("HpassMET45_ge2b","pass MET45 versus MET (RA2 && >=2b)",nbins,met_min,met_max,2,0,2);

   TH2D HdeltaPhib1b2_MET("HdeltaPhib1b2_MET","DeltaPhi(b1,b2) v MET (>=2b)",nbins,met_min,met_max,nbins,0,pi);
   TH2D HminDeltaPhiMETb_MET_ge2b("HminDeltaPhiMETb_MET_ge2b","minDeltaPhi(b,MET) v MET (>=2b)",nbins,met_min,met_max,nbins,0,pi);
   TH2D HminDeltaPhiMETj_MET_ge2b("HminDeltaPhiMETj_MET_ge2b","minDeltaPhi(j,MET) v MET (>=2b)",nbins,met_min,met_max,nbins,0,pi);

   TH2D HminDeltaPhiMETb_HminDeltaPhiMETj_ge2b("HminDeltaPhiMETb_HminDeltaPhiMETj_ge2b","minDPhi(b,MET) v minDPhi(j,MET) (>=2b)",
					       nbins,0,pi,nbins,0,pi);
   TH2D HdeltaPhiMPTMET_MET_ge2b("HdeltaPhiMPTMET_MET_ge2b","DeltaPhi(MET,MPT) v MET (RA2 && >=2b)",nbins,met_min,met_max,nbins,0,pi);

   TH1D HtopDecayCategory_ge2b("HtopDecayCategory_ge2b","top decay category (RA2 && >=2b)",nTopCategories-1,0.5,nTopCategories-0.5);

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
//    TH2D HpassHT100U("HpassHT100U","pass HT100U versus HT (RA2)",ht_nbins,ht_min,ht_max,2,0,2);
//    TH2D HpassHT100U_ge1b("HpassHT100U_ge1b","pass HT100U versus HT (RA2 && >=1b)",ht_nbins,ht_min,ht_max,2,0,2);
//    TH2D HpassHT100U_ge2b("HpassHT100U_ge2b","pass HT100U versus HT (RA2 && >=2b)",ht_nbins,ht_min,ht_max,2,0,2);

   TH1D H_HT("H_HT","HT",ht_nbins,ht_min,ht_max);
   TH1D H_HT_ge1b("H_HT_ge1b","HT (>=1b)",ht_nbins,ht_min,ht_max);
   TH1D H_HT_ge2b("H_HT_ge2b","HT (>=2b)",ht_nbins,ht_min,ht_max);

   //as usual, perhaps we should manage our histos with HistHolder (but let's keep it simple instead)
   Hnjets.Sumw2();
   Hnjets_nocuts.Sumw2();
   Hnjets_ge1b.Sumw2();
   Hnjets_ge2b.Sumw2();
   Hnjets_ge3b.Sumw2();

   HdeltaPhiMPTMET.Sumw2();
   HdeltaPhiMPTMET_ge1b.Sumw2();
   HdeltaPhiMPTMET_ge2b.Sumw2();
   HdeltaPhib1b2_minDeltaPhiMETb.Sumw2();
   HminDeltaPhiMETb_ge1b.Sumw2();
   HminDeltaPhiMETj_ge1b.Sumw2();
   HminDeltaPhiMETj.Sumw2();

   HminDeltaPhiMETb_ge2b.Sumw2();
   HminDeltaPhiMETj_ge2b.Sumw2();

   HVminDeltaPhiMETj.Sumw2();
   HVminDeltaPhiMETj_ge1b.Sumw2();
   HVminDeltaPhiMETj_ge2b.Sumw2();

   HVminDeltaPhiMHTj.Sumw2();
   HVminDeltaPhiMHTj_ge1b.Sumw2();
   HVminDeltaPhiMHTj_ge2b.Sumw2();

   HV2minDeltaPhiMETj.Sumw2();
   HV2minDeltaPhiMETj_ge1b.Sumw2();
   HV2minDeltaPhiMETj_ge2b.Sumw2();

   Hjetpt1.Sumw2();
   Hjetpt1_ge1b.Sumw2();
   Hjetpt1_ge2b.Sumw2();

   Hbjetpt1_ge1b.Sumw2();
   Hbjetpt1_ge2b.Sumw2();

   H_MHT.Sumw2();
   H_MET.Sumw2();
   H_MHT_ge1b.Sumw2();
   H_MET_ge1b.Sumw2();
   H_MHT_ge2b.Sumw2();
   H_MET_ge2b.Sumw2();
   H_MHT_ge3b.Sumw2();
   H_MET_ge3b.Sumw2();

   H_HT.Sumw2();
   H_HT_ge1b.Sumw2();
   H_HT_ge2b.Sumw2();

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

//    HpassHT100U.Sumw2();
//    HpassHT100U_ge1b.Sumw2();
//    HpassHT100U_ge2b.Sumw2();

//    HpassMET45.Sumw2();
//    HpassMET45_ge1b.Sumw2();
//    HpassMET45_ge2b.Sumw2();

   //keep track of performance
   TDatime starttime; //default ctor is for current time

   //event loop
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries ;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) {assert(0);}
      nb = GetEntry(jentry);   nbytes += nb; //use member function GetEntry instead of fChain->

      Hnjets_nocuts.Fill( tightJetIndex->size(), weight );
      if (Cut(ientry) < 0) continue; //jmt use cut

      //calculate things
      double dp_MPTMET = getDeltaPhiMPTMET();

      Hnjets.Fill( tightJetIndex->size(), weight );
      HdeltaPhiMPTMET.Fill( dp_MPTMET,weight);
      H_MHT.Fill(MHT ,weight);
      H_MET.Fill( getMET() ,weight);

      Hjetpt1.Fill(jetPt.at(0),weight); //FIXME

      H_HT.Fill(HT,weight);

      //FIXME -- this will not work in general
//       int passHTtrig = passTrigger->at(0) ? 1:0; //index 0 is HLT_HT100U
//       int passMETtrig = passTrigger->at(2) ? 1:0; //index 2 is HLT_MET45
//       HpassHT100U.Fill(HT,passHTtrig); //note -- NOT using weight here!
//       HpassMET45.Fill(MET,passMETtrig);

      double minDeltaPhi_j_MET= getMinDeltaPhiMET(3);
      double minDeltaPhi_j_MHT= getMinDeltaPhiMHT(3);

      HminDeltaPhiMETj.Fill(minDeltaPhi_j_MET,weight);
      HVminDeltaPhiMETj.Fill(minDeltaPhi_j_MET,weight);
      HVminDeltaPhiMHTj.Fill(minDeltaPhi_j_MHT,weight);
      HV2minDeltaPhiMETj.Fill(minDeltaPhi_j_MET,weight);

      if ( nbSSVM < 1) continue; //cut on the number of b tags
      Hnjets_ge1b.Fill( tightJetIndex->size(), weight );
      H_MHT_ge1b.Fill(MHT ,weight);
      H_MET_ge1b.Fill(getMET() ,weight);

      Hjetpt1_ge1b.Fill(jetPt.at(0),weight);

      HdeltaPhiMPTMET_ge1b.Fill( dp_MPTMET,weight);

//       HpassHT100U_ge1b.Fill(HT,passHTtrig); //note -- NOT using weight here!
//       HpassMET45_ge1b.Fill(MET,passMETtrig);
      H_HT_ge1b.Fill(HT,weight);

      double minDeltaPhi_b_MET= getMinDeltaPhibMET();

      HminDeltaPhiMETb_ge1b.Fill(minDeltaPhi_b_MET,weight);
      HminDeltaPhiMETj_ge1b.Fill(minDeltaPhi_j_MET,weight);
      HVminDeltaPhiMETj_ge1b.Fill(minDeltaPhi_j_MET,weight);
      HVminDeltaPhiMHTj_ge1b.Fill(minDeltaPhi_j_MHT,weight);
      HV2minDeltaPhiMETj_ge1b.Fill(minDeltaPhi_j_MET,weight);

      int topcat = isData_ ? -99 : getTopDecayCategory(); //no sense in looking at MC truth in the data
      double minDeltaR_bj=999;
      int nbjetsfound=0;
      double bjetpt1=0;
      //note that all tight jet vectors should have the same size
      for (unsigned int ib = 0; ib< jetPhi.size(); ib++) {
	if ( passBCut(ib)) { //refind the b jets
	  nbjetsfound++;
	  if (nbjetsfound==1) bjetpt1=jetPt.at(ib); //if this is the *lead* b jet

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
      Hbjetpt1_ge1b.Fill(bjetpt1,weight);

      HminDeltaR_bj_vMET_ge1b.Fill(getMET(), minDeltaR_bj,weight);
      HminDeltaR_bj_vMET_ABCD_ge1b.Fill(getMET(), minDeltaR_bj,weight);

      if ( nbSSVM < 2) continue; //cut on the number of b tags
      Hbjetpt1_ge2b.Fill(bjetpt1,weight);

      Hnjets_ge2b.Fill( tightJetIndex->size(), weight );
      H_MHT_ge2b.Fill(MHT ,weight);
      H_MET_ge2b.Fill(getMET() ,weight);
      HdeltaPhiMPTMET_ge2b.Fill( dp_MPTMET,weight);

      Hjetpt1_ge2b.Fill(jetPt.at(0),weight);

//       HpassHT100U_ge2b.Fill(HT,passHTtrig);
//       HpassMET45_ge2b.Fill(MET,passMETtrig);
      H_HT_ge2b.Fill(HT,weight);

      HminDeltaPhiMETb_ge2b.Fill(minDeltaPhi_b_MET,weight);
      HminDeltaPhiMETj_ge2b.Fill(minDeltaPhi_j_MET,weight);
      HVminDeltaPhiMETj_ge2b.Fill(minDeltaPhi_j_MET,weight);
      HVminDeltaPhiMHTj_ge2b.Fill(minDeltaPhi_j_MHT,weight);
      HV2minDeltaPhiMETj_ge2b.Fill(minDeltaPhi_j_MET,weight);

      double deltaPhi_b1b2 = getDeltaPhib1b2();
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
      Hnjets_ge3b.Fill( tightJetIndex->size(), weight );
      H_MHT_ge3b.Fill(MHT ,weight);
      H_MET_ge3b.Fill(getMET() ,weight);
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

  /* LM0 events */
  /*
    specifyEvent(1, 6, 817);
    specifyEvent(1, 186, 25463);
  */

  ULong64_t nfound=0;

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries(); //jmt: remove Fast
   const TString sp=" ";
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = GetEntry(jentry);   nbytes += nb; //use member function GetEntry instead of fChain->

      if (false) {
	if (eventIsSpecified() ) {
	  nfound++;
	  cout<<"--- record for event: ("<<jentry <<") run,ls,ev = "<<runNumber<<sp<<lumiSection<<sp<<eventNumber<<endl;
	  //Show() doesn't cut it-- it just shows the pointer addresses!
	  
	  if (true) {
	    for (unsigned int i=0 ; i<cutTags_.size(); i++) {
	      if (cutRequired(cutTags_[i])) {
		TString passstr = passCut(cutTags_[i]) ? "pass" : "fail";
		cout<< cutNames_[cutTags_[i]] <<sp<<passstr<<endl;
	      }
	    }
	  }
	  
	  if (true) {
	    cout<<" jet info (pT, Eta, hadFrac, isGood) n good jets = "<<nGoodJets_Sync1()<<endl;
	    for (unsigned int ijet=0; ijet<loosejetPt->size(); ijet++) {
	      TString jetisgood = isGoodJet_Sync1(ijet) ? "Good" : "notGood";
	      cout<<"\tjet "<<ijet<<": "<<loosejetPt->at(ijet)<<sp<<loosejetEta->at(ijet)<<sp<<loosejetEnergyFracHadronic->at(ijet)<<sp<<jetisgood<<endl;
	    }
	  }
	}
	if (nfound == specifiedEvents_.size()) break; //save time at the end
      }

      //in case we just want to dump every event number
      if (true) {
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

