
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TCut.h"
#include "TMath.h"

/*
should replace this code with simple call to TTree::CopyTree("cut string");
*/

#include "TDirectory.h"

//#include "TH2.h"
//#include "TH3.h"

#include <iostream>

//
// Note: If you are using this on a SLC6 machine, you need to do
//       this in your interactive root session
//
//          gSystem->AddIncludePath(" -D__USE_XOPEN2K8 ")
//
//       before doing this
//
//          .L doSkimSlim.C+
//

void DelphesSkim( const char* infile_name ) {

      TFile* infile = new TFile( infile_name ) ;
      if ( ! (infile->IsOpen()) ) return ;

      TTree* inReducedTree = (TTree*) infile->Get("simpleTree") ;

      printf("\n\n Number of entries: %llu\n\n", inReducedTree->GetEntries() ) ;

      TString filename(infile_name);

      TString skimoption=TString( gSystem->Getenv("SKIMOPTION"));
      bool skim_on_mll = skimoption.Contains("MLL");
      bool skim_on_ht = skimoption.Contains("TIGHTHT");

      TCut selectionCuts="(1)";
      if (skim_on_mll) {
	cout<<"skimming on Mll"<<endl;
	selectionCuts = "mll>0 || mll_loose >0 || mll_veryloose>0";
      }
      else      if (skim_on_ht) {
	cout<<"skimming on HT"<<endl;
	selectionCuts = selectionCuts && "HT>1000";
      }
      else {cout<<"Need to set one selection cut or the other via envvar"<<endl; assert(0);}
      TString      selectionString = selectionCuts.GetTitle();
      cout<<"cut string: "<<selectionString<<endl;

      //--- Open output file
      TString outfile_name( infile_name ) ;

      outfile_name.ReplaceAll( ".root", "-skim.root" ) ;

      if ( outfile_name.CompareTo( infile_name ) == 0 ) {
         printf("\n\n *** Input and output file names same.  Input doesn't contain .root in name?\n") ;
         printf("    input: %s \n", infile_name ) ;
         printf("    output: %s \n", outfile_name.Data() ) ;
         return ;
      }
      printf("\n\n Output file: %s\n\n", outfile_name.Data() ) ;
      char command[10000] ;
      sprintf( command, "ls %s >& /dev/null", outfile_name.Data() ) ;
      int returnstat = gSystem->Exec( command ) ;
      if ( returnstat == 0 ) {
         char mvfile[10000] ;
         sprintf( mvfile, "%s-old", outfile_name.Data() ) ;
         printf("\n\n *** Output file already exists.  Moving it to %s\n\n", mvfile ) ;
         sprintf( command, "mv %s %s", outfile_name.Data(), mvfile ) ;
         gSystem->Exec( command ) ;
      }
      TFile* outfile = new TFile( outfile_name, "recreate" ) ;
      
      //      inReducedTree -> SetBranchStatus("*",1) ; // enable all branches.

      //--- Vars needed to decide whether or not to save the event.
      /*
      float mll,mll_loose,mll_veryloose,HT ;
      inReducedTree -> SetBranchAddress("HT", &HT) ;
      inReducedTree -> SetBranchAddress("mll", &mll) ;
      inReducedTree -> SetBranchAddress("mll_loose", &mll_loose) ;
      inReducedTree -> SetBranchAddress("mll_veryloose", &mll_veryloose) ;
      */      

      outfile->cd();
      TTree* outReducedTree=inReducedTree->CopyTree(selectionString);   
   //      TTree* outReducedTree = inReducedTree->CloneTree(0) ;

      //	 bool pass_mll = (mll>0) || (mll_loose>0) || (mll_veryloose>0);
      //	 bool pass_ht = HT>=1000;

      printf("\n\n\n Done.\n\n\n") ;
      printf("%.2f percent selected by skim\n",100*double(outReducedTree->GetEntries()) / double(inReducedTree->GetEntries()));
      printf("\n\n Output file:  %s\n\n\n", outfile_name.Data() ) ;

      outReducedTree->Write();
      //      delete outfile ;
      outfile->Close();
  }

