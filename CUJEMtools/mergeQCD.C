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
  //input path
  const TString path="CUJEM_r1.4/"; //adjust this as needed

  vector<float> weights;
  vector<TString> qcd;

  const TString kindofsample="pythia";
  if (kindofsample=="madgraph") {
    cout<<"Using madgraph samples!"<<endl;
    //weights -- adjust as needed
    weights.push_back(66.0342); //these are weights for 100 pb^-1
    weights.push_back(3.703);
    weights.push_back(0.1292);
    weights.push_back(0.005181);  
    
    //file names -- adjust as needed
    qcd.push_back( "plots_QCD_Pt100to250-madgraph_");
    qcd.push_back( "plots_QCD_Pt250to500-madgraph_");
    qcd.push_back("plots_QCD_Pt500to1000-madgraph_");
    qcd.push_back("plots_QCD_Pt1000toInf-madgraph_");
  }
  else if (kindofsample=="pythia") {
    cout<<"Using pyhtia samples!"<<endl;

    weights.push_back(28.8377); //these are weights for 100 pb^-1
    weights.push_back(0.813167135);

    qcd.push_back("plots_QCD_Pt80_");
    qcd.push_back("plots_QCD_Pt170_");

  }
  //  float total=0;
  //  for (unsigned int i=0; i<weights.size();i++)     total+=weights.at(i);

  vector<TFile*> fqcd;
  for (unsigned int i=0;i<weights.size();i++) {
    TString filename=qcd.at(i);
    filename+= cutstring;
    filename+=".root";
    filename.Prepend(path);
    TFile* ftmp=new TFile(filename);
    if (ftmp->IsZombie()) cout<<"Problem loading file with name="<<filename<<endl;
    fqcd.push_back( ftmp);
  }
  TString outname="plots_QCD_";
  outname+=cutstring;
  outname+=".root";
  TFile fout(outname,"RECREATE");

  //we've got all of the files open, so let's loop over the keys in a file and add them up!

  TIter nextkey(fqcd[0]->GetListOfKeys());
  TKey *key;

  while (key = (TKey*) nextkey() ) {
    //    cout<<"-- new key --"<<endl;
    TH1F* th1=0;
    TH2F* th2=0;
    TH1F* h1=0;
    TH2F* h2=0;
    if ( key->ReadObj()->InheritsFrom("TH2F") ) {
      //      cout<<"2D: Key="<<key->ReadObj()->GetName()<<endl;
      th2=(TH2F*) key->ReadObj();
    }
    else if ( key->ReadObj()->InheritsFrom("TH1F") ) {
      //      cout<<"1D: Key="<<key->ReadObj()->GetName()<<endl;
      th1=(TH1F*) key->ReadObj();
    }
    else {
      cout<<"Problem determining type of key="<<key->ReadObj()->GetName()<<endl;
      continue;
    }

    fout.cd();
    //    cout<<"-- about to clone --"<<endl;
    if (th1!=0) { 
      h1=(TH1F*) th1->Clone();
      h1->Sumw2();
      h1->Scale(weights.at(0));
    }
    else if (th2!=0) {
      h2=(TH2F*) th2->Clone();
      h2->Sumw2();
      h2->Scale(weights.at(0));
    }

    //    cout<<" -- about to add subsequent histos --"<<endl;
    for (unsigned int i=1; i<fqcd.size(); i++) {

      if (h1!=0) {
	th1 = (TH1F*) fqcd.at(i)->Get(th1->GetName());
	//cout<<th1<<endl;
	h1->Add(th1,weights.at(i));
      }
      else if (h2!=0) {
	th2 = (TH2F*) fqcd.at(i)->Get(th2->GetName());
	//cout<<th2<<endl;
	h2->Add(th2,weights.at(i));
      }

    }

    fout.cd();
    if (h1!=0)      h1->Write();
    else if (h2!=0) h2->Write();
  }

  cout<<"Finished writing file "<<outname<<endl;
  fout.Close();

}
