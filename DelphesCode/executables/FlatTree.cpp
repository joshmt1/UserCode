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

  if(argc != 2 && argc !=3)
  {
    cout << " Usage: " << appName << " input_file [output_file]" << endl;
    cout << " input_file - input file in ROOT format ('Delphes' tree)," << endl;
    cout << " output_file - optional name of output file" << endl;
    return 1;
  }

  gROOT->SetBatch();

  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);

  TString inputFile(argv[1]);
  TString outputFile = "simpleTree.root";
  if (argc==3) outputFile = TString(argv[2]);

//------------------------------------------------------------------------------

// Here you call your macro's main function 

  FlatTree(inputFile,outputFile);

//------------------------------------------------------------------------------

}


