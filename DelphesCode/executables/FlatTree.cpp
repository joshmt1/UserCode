#include <iostream>
#include <utility>
#include <vector>

#include <cassert>

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"

#include "TString.h"

#include "TStopwatch.h"
//#include "TH2.h"
//#include "THStack.h"
//#include "TLegend.h"
//#include "TPaveText.h"
#include "TClonesArray.h"
//#include "TLorentzVector.h"

#include "external/mt2analysis/SimpleTree.hh"
#include "external/mt2analysis/mt2_bisect.h"
#include "external/mt2analysis/FindHemispheres.hh"
#include "external/mt2analysis/Hemisphere.hh"
#include "external/mt2analysis/Utilities.hh"
#include "external/mt2analysis/CrossSections.hh"
#include "external/mt2analysis/McTruthInfo.hh"
#include "external/mt2analysis/MllComputer.hh"
#include "external/mt2analysis/JetLeptonCleaner.hh"
#include "external/mt2analysis/DelWeight.h"

#include "classes/DelphesClasses.h"

#include "ExRootAnalysis/ExRootTreeReader.h"
//#include "ExRootAnalysis/ExRootTreeWriter.h"
//#include "ExRootAnalysis/ExRootTreeBranch.h"
//#include "ExRootAnalysis/ExRootResult.h"
//#include "ExRootAnalysis/ExRootUtilities.h"

using namespace std;

//------------------------------------------------------------------------------

// Here you can put your analysis macro

#include "FlatTree.C"

//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char *appName = "FlatTree";

  if(argc < 2 )
  {
    cout << " Usage: " << appName << " input_file [-O output_file] [-N x y]" << endl;
    cout << " input_file - input file in ROOT format ('Delphes' tree)," << endl;
    cout << " output_file - optional name of output file" << endl;
    cout << " x,y - optional job splitting parameters. This is job x of total y" << endl;
    return 1;
  }

  gROOT->SetBatch();

  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);

  TString inputFile(argv[1]);
  TString outputFile = "simpleTree.root";

  int x=1,y=1;
  bool docleaning=true;
  bool usepujetid=false;
  for (int iarg = 2; iarg<argc; ) {

    TString arg = argv[iarg];
    if (arg=="-O") {
      outputFile = argv[iarg+1];
      iarg+=2;
    }
    else if (arg=="-N") {
      x = TString(argv[iarg+1]).Atoi();
      y = TString(argv[iarg+2]).Atoi();
      iarg+=3;
    }
    else if (arg=="-nocleaning") {
      docleaning=false;
      iarg++;
    }
    else if (arg=="-pujetid") {
      usepujetid=true;
      iarg++;
    }
    else {
      cout << " Usage: " << appName << " input_file [-O output_file] [-N x y]" << endl;
      assert(0);
    }

  }

//------------------------------------------------------------------------------

// Here you call your macro's main function 

  FlatTree(inputFile,outputFile,x,y,docleaning,usepujetid);

//------------------------------------------------------------------------------

}

