#include "TString.h"
#include "TFile.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TCollection.h"

#include <iostream>

using namespace std;

void mergeQCD( TString cutstring="CutsJMTNoJetVeto" )
{

  float weights[] = {66.0342,3.703,0.1292,0.005181}; //these are weights for 100 pb^-1

  float total=0;
  for (int i=0; i<4;i++)     total+=weights[i];

  TString path="CUJEM_r1.3/";

  TString qcd[]={"","","",""};
  qcd[0] = "plots_QCD_Pt100to250-madgraph_";
  qcd[1] = "plots_QCD_Pt250to500-madgraph_";
  qcd[2] = "plots_QCD_Pt500to1000-madgraph_";
  qcd[3] = "plots_QCD_Pt1000toInf-madgraph_";
  
  TFile* fqcd[4];

  for (int i=0;i<4;i++) {
    qcd[i] += cutstring;
    qcd[i] += ".root";
    qcd[i].Prepend(path);
    fqcd[i] = new TFile(qcd[i]);
  }

  TString outname="plots_QCD_";
  outname+=cutstring;
  outname+=".root";
  TFile fout(outname,"RECREATE");

  //we've got all of the files open, so let's loop over the keys in a file and add them up!
  //TList* keys= fqcd[0]->GetListOfKeys();

  TIter nextkey(fqcd[0]->GetListOfKeys());
  TKey *key;

  while (key = (TKey*) nextkey() ) {
    cout<<"-- new key --"<<endl;
    TH1F* th1=0;
    TH2F* th2=0;
    TH1F* h1=0;
    TH2F* h2=0;
    if ( key->ReadObj()->InheritsFrom("TH2F") ) {
      cout<<"2D: Key="<<key->ReadObj()->GetName()<<endl;
      th2=(TH2F*) key->ReadObj();
    }
    else if ( key->ReadObj()->InheritsFrom("TH1F") ) {
      cout<<"1D: Key="<<key->ReadObj()->GetName()<<endl;
      th1=(TH1F*) key->ReadObj();
    }
    else {
      cout<<"Problem determining type of key="<<key->ReadObj()->GetName()<<endl;
      continue;
    }

    fout.cd();
    cout<<"-- about to clone --"<<endl;
    if (th1!=0) { 
      h1=(TH1F*) th1->Clone();
      h1->Sumw2();
      h1->Scale(weights[0]);
    }
    else if (th2!=0) {
      h2=(TH2F*) th2->Clone();
      h2->Sumw2();
      h2->Scale(weights[0]);
    }

    cout<<" -- about to add subsequent histos --"<<endl;
    for (int i=1; i<4; i++) {

      //this seems to work! we need to add the weights in here now!!

      if (h1!=0) {
	th1 = (TH1F*) fqcd[i]->Get(th1->GetName());
	//cout<<th1<<endl;
	h1->Add(th1,weights[i]);
      }
      else if (h2!=0) {
	th2 = (TH2F*) fqcd[i]->Get(th2->GetName());
	//cout<<th2<<endl;
	h2->Add(th2,weights[i]);
      }

    }

    fout.cd();
    if (h1!=0)      h1->Write();
    else if (h2!=0) h2->Write();
  }

  //    TH2F* h2=dynamic_cast<TH2F*> keys->At(i);
  cout<<"Finished write file "<<outname<<endl;
  fout.Close();

}
