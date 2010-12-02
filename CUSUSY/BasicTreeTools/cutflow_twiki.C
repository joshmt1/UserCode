#include "TString.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>

/*
usage:
root -b -l -q cutflow_twiki.C++
*/
//#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TAxis.h"
TCanvas * Ccutflow;

//global settings that change what the output of the code is

const TString filestub_ ="Baseline0_PF_pfMEThigh_PFLep_minDP_NoDeltaPhi_NoTrigger";
const int mode_ = 1; //mode 1 is print cut flow table ; mode 2 is print S/sqrt(S+B) table ; mode 3 is S/sqrt(B)
const bool latexMode_ = true; //otherwise TWiki


//utility function for making output more readable
TString format_nevents(double n,double e) {
  char out[100];
  if (e>=1  ||  e < 0.00001) {
    sprintf(out,"%.0f +/- %.0f",n,e);
  }
  else {
    int nfig = ceil(fabs(log10(e)));
    TString form="%.";
    form+=nfig; form+="f +/- %.";
    form+=nfig; form+="f";
    sprintf(out,form.Data(),n,e);
  }
  return TString(out);
}

//i can't think of an elegant solution, other than to just run this each time
void combineSamples(TString outputname, std::vector<TString> inputList)
{
  /*
idea here is to open the input files, combine them, and write out a combined file.
in this way, the main code can read the combined sample without any special treatment
This is exactly how QCD used to be treated in the main function, but now I want to generalize it and write the file too
  */

  unsigned int nToCombine = inputList.size();
  std::vector<TString> cutnames;
  std::map<TString, std::vector<float> > n;
  std::map<TString, std::vector<float> > err;

  //first load the event counts and errors for the input samples
  for (unsigned int i=0 ; i<nToCombine; i++) {
    TString filename = "cutflow.";
    filename+=filestub_;
    filename+="."; filename+=inputList[i];
    filename += ".dat";

    cout<<"Reading "<<filename<<endl;
    float nevt,nevterr;
    char cutdesc[100];
    ifstream file(filename.Data());
    if (!file.good()) {cout<<"bad file! "<<filename<<endl; assert(0);}
    while (file>>cutdesc>>nevt>>nevterr ) {
      if (i==0)  cutnames.push_back(cutdesc);
      n[TString(inputList[i])].push_back(nevt);
      err[TString(inputList[i])].push_back(nevterr);
    }
    file.close();
  }

  //open a file for output
  TString outfilename="cutflow."; 
  outfilename +=filestub_;
  outfilename +=".";
  outfilename +=outputname;
  outfilename +=".dat";
  remove(outfilename.Data()); //if it already exists, delete it
  ofstream ofile(outfilename.Data());

  unsigned int ncuts = cutnames.size();
  for (unsigned int i=0; i<ncuts; i++) {
    float total=0;
    float total_err=0;

    for (unsigned int j=0 ; j<nToCombine; j++) {
      total += n[inputList[j]].at(i);
      total_err += pow(err[inputList[j]].at(i),2);
    }

    total_err = sqrt(total_err);
    //now we've got to write the output file
    ofile<<cutnames.at(i)<<"\t"<<total<<"\t"<<total_err<<endl;
  }
  ofile.close();

}

void cutflow_twiki()
{
  assert(mode_==1 || mode_==2 || mode_ ==3);

  TString modeDescription="";
  if (mode_==2) modeDescription="S/sqrt(S+B)";
  else if (mode_==3) modeDescription="S/sqrt(B)";

  int icolor=1;

  //the cut names are no longer hard-coded here, but are rather read in from the cutflow files
  //so old cutflow files won't work any more!
  std::vector<TString> cutnames;

  //is this really the best way to do this?
  int nqcd = 4;
  char *qcd_list[]={"QCD100","QCD250","QCD500","QCD1000"};
  int nbackground = 5;
  char *background_list[]={"TTbarJets","SingleTop","Zinvisible","WJets","ZJets"};
  int nsignal = 5;//16; //oops, where did LM3 go?
  //  char *signal_list[]={"LM0", "LM1", "LM2", "LM4", "LM5", "LM6","LM7", "LM8","LM9","LM9p", "LM9t175", "LM10", "LM11", "LM12","LM13","mMSSM"};
  char *signal_list[]={"LM0", "LM1", "LM9","LM13","mMSSMv3"};

  //in principle we can also rewrite the code using this to combine QCD
  std::vector<TString> singletopnames;
  singletopnames.push_back("SingleTop-tChannel");
  singletopnames.push_back("SingleTop-tWChannel");
  combineSamples("SingleTop",singletopnames);

  //these are only relevant for mode==1
  const int drawsignalindex = 3; //this is the signal to use for the plot
  const TString mysig = signal_list[drawsignalindex];

  //as long as i compile, I can use the most basic stl containers
  // ... what was I thinking? what i really want is a map of these, indexed by the names.
  //can ROOT handle it? //seems ok....

  std::map<TString, std::vector<float> > qcd;
  std::map<TString, std::vector<float> > qcderr;

  std::map<TString, std::vector<float> > background;
  std::map<TString, std::vector<float> > backgrounderr;
  //probably I could have used a TH1 as a storage container from the beginning
  //but keep the old STL structure in order to avoid bugs
  std::map<TString, TH1D*> Hcutflow;

  std::map<TString, std::vector<float> > signal;
  std::map<TString, std::vector<float> > signalerr;

  //lordy, this could really be more generic!
  char cutdesc[100];
  for (int i=0 ; i<nqcd; i++) {
    TString filename = "cutflow.";
    filename+=filestub_;
    filename+="."; filename+=qcd_list[i];
    filename += ".dat";

    cout<<"Reading "<<filename<<endl;
    float nevt,nevterr;
    ifstream file(filename.Data());
    if (!file.good()) {cout<<"bad file! "<<filename<<endl; return;}
    while (file>>cutdesc>>nevt>>nevterr ) {
      if (i==0)  cutnames.push_back(cutdesc);
      qcd[TString(qcd_list[i])].push_back(nevt);
      qcderr[TString(qcd_list[i])].push_back(nevterr);
    }
    file.close();
  }

  if (mode_==1) {
    Hcutflow[TString("QCD")] = new TH1D("HQCD","QCD",cutnames.size(),0,cutnames.size());
    Hcutflow["QCD"]->SetLineColor(icolor);   Hcutflow["QCD"]->SetMarkerColor(icolor);
    Hcutflow["QCD"]->SetFillColor(icolor);
    icolor++;
  }

  if (mode_==1) {
    Hcutflow[mysig] = new TH1D("Hsignal",mysig,cutnames.size(),0,cutnames.size());
    Hcutflow[mysig]->SetLineColor(icolor);   Hcutflow[mysig]->SetMarkerColor(icolor);
    Hcutflow[mysig]->SetFillColor(icolor);
    icolor++;
  }

  for (int i=0 ; i<nbackground; i++) {
    TString filename = "cutflow.";
    filename+=filestub_;
    filename+="."; filename+=background_list[i];
    filename += ".dat";

    if (mode_==1) {
      Hcutflow[TString(background_list[i])] = new TH1D("H"+TString(background_list[i]),TString(background_list[i]),cutnames.size(),0,cutnames.size());
      cout<<"Created histogram name = "<<Hcutflow[TString(background_list[i])]->GetName()<<endl;
      Hcutflow[TString(background_list[i])]->SetLineColor(icolor); 
      Hcutflow[TString(background_list[i])]->SetMarkerColor(icolor);
      Hcutflow[TString(background_list[i])]->SetFillColor(icolor);
      icolor++;
    }
    
    cout<<"Reading "<<filename<<endl;
    float nevt,nevterr;
    ifstream file(filename.Data());
    if (!file.good()) {cout<<"bad file! "<<filename<<endl; return;}
    int ibin=1;
    while (file>>cutdesc>>nevt>>nevterr ) {
      background[TString(background_list[i])].push_back(nevt);
      backgrounderr[TString(background_list[i])].push_back(nevterr);
      if (mode_==1) {
	Hcutflow[TString(background_list[i])]->SetBinContent(ibin,nevt);
	Hcutflow[TString(background_list[i])]->SetBinError(ibin,nevterr);
	ibin++;
      }
    }
    file.close();
    //    cout<<  background[TString(background_list[i])].size()<<"\t"<<backgrounderr[TString(background_list[i])].size()<<endl;
  }
  cout<<"--"<<endl;
  cout<<  background.size()<<"\t"<<  backgrounderr.size()<<endl;

  for (int i=0 ; i<nsignal; i++) {
    TString filename = "cutflow.";
    filename+=filestub_;
    filename+="."; filename+=signal_list[i];
    filename += ".dat";

    if (mode_==2 ||mode_ ==3) {
      Hcutflow[TString(signal_list[i])] = new TH1D("H"+TString(signal_list[i]),TString(signal_list[i]),cutnames.size(),0,cutnames.size());
      cout<<"Created histogram name = "<<Hcutflow[TString(signal_list[i])]->GetName()<<endl;
      Hcutflow[TString(signal_list[i])]->SetLineColor(icolor); 
      Hcutflow[TString(signal_list[i])]->SetLineWidth(2); 
      Hcutflow[TString(signal_list[i])]->SetMarkerColor(icolor);
      //      Hcutflow[TString(signal_list[i])]->SetFillColor(icolor); //don't want this for modes 2 and 3
      icolor++;
    }
    
    cout<<"Reading "<<filename<<endl;
    float nevt,nevterr;
    ifstream file(filename.Data());
    if (!file.good()) {cout<<"bad file! "<<filename<<endl; return;}
    while (file>>cutdesc>>nevt>>nevterr ) {
      signal[TString(signal_list[i])].push_back(nevt);
      signalerr[TString(signal_list[i])].push_back(nevterr);
    }
    file.close();
    //    cout<<  signal[TString(signal_list[i])].size()<<"\t"<<signalerr[TString(signal_list[i])].size()<<endl;
  }
  cout<<"--"<<endl;
  cout<<  signal.size()<<"\t"<< signalerr.size()<<endl;
  cout<<"--"<<endl;

  //then we need to translate to a one cut at a time table
  int ncuts = signal[TString(signal_list[0])].size(); //should be the same in all cases!
  cout<<"ncuts = "<<ncuts<<endl;
  const TString pm = latexMode_ ? " \\pm " : " +/- ";
  const  TString col = latexMode_ ? " & " : " | ";

  if (!latexMode_)  cout<<col<<"n" <<col;
  cout<<"Cut"<<col;
  if (mode_==1) {
    cout<<"QCD"<<col;
    for (int ibackground=0 ; ibackground<nbackground; ibackground++)     cout<< background_list[ibackground]<<col;
  }
  for (int isignal=0 ; isignal<nsignal; isignal++)   {
    cout<<signal_list[isignal];
    if (!(latexMode_ && isignal==nsignal-1))  cout<<col;
    else cout<<" \\\\";
  }
  cout<<endl;


  for (int i=0; i<ncuts; i++) {
    if (!latexMode_)    cout<<col<< i<<col;
    cout<<cutnames.at(i)<<col;

    //cout<<setprecision(1);

    float qcd_total=0;
    float qcd_total_err=0;
    //it would be more elegant to iterate over the map...but to hell with elegance
    for (int iqcd=0 ; iqcd<nqcd; iqcd++) {
      //      cout<<     qcd[qcd_list[iqcd]].at(i) << " | ";
      qcd_total += qcd[qcd_list[iqcd]].at(i);
      qcd_total_err += pow(qcderr[qcd_list[iqcd]].at(i),2);
    }
    if (mode_ ==1) {
      Hcutflow["QCD"]->SetBinContent(i+1,qcd_total);
      Hcutflow["QCD"]->SetBinError(i+1,sqrt(qcd_total_err));
    }

    float background_total=qcd_total;
    float background_total_err=qcd_total_err; //before taking the square root

    qcd_total_err = sqrt(qcd_total_err);
    if (mode_==1) cout<< format_nevents(qcd_total , qcd_total_err)<<col;

    for (int ibackground=0 ; ibackground<nbackground; ibackground++) {
      if (mode_==1)
	cout<< format_nevents(background[background_list[ibackground]].at(i) ,backgrounderr[background_list[ibackground]].at(i)) << col;
      background_total += background[background_list[ibackground]].at(i);
      background_total_err += pow(backgrounderr[background_list[ibackground]].at(i),2);
    }

    //now that we've summed the squares, take the square root
    background_total_err = sqrt(background_total_err);

    for (int isignal=0 ; isignal<nsignal; isignal++) {
      if (mode_==1)   {
	cout<<format_nevents(signal[signal_list[isignal]].at(i) , signalerr[signal_list[isignal]].at(i)) ;
	if (!(latexMode_ && isignal==nsignal-1)) cout<< col; //get rid of trailing & in latex mode
	else cout<<" \\\\";
      }
      else if (mode_==2) {
	if (signal[signal_list[isignal]].at(i)+ background_total >0)
	  cout<<signal[signal_list[isignal]].at(i)/sqrt(signal[signal_list[isignal]].at(i)+ background_total)  << col;
	else  cout<<" - "<<endl;
      }
      else if (mode_==3) {
	if ( background_total >0)
	  cout<<signal[signal_list[isignal]].at(i)/sqrt(background_total)  << col;
	else  cout<<" - "<<endl;
      }
      if (mode_ ==1 && isignal == drawsignalindex) {
	  Hcutflow[mysig]->SetBinContent(i+1,signal[signal_list[isignal]].at(i));
	  Hcutflow[mysig]->SetBinError(i+1,signalerr[signal_list[isignal]].at(i));
      }
      else if (mode_==2) {
	if (signal[signal_list[isignal]].at(i)+ background_total >0) {
	  //not so keen on the fact that I'm duplicating the calculation code here...
	  Hcutflow[signal_list[isignal]]->SetBinContent(i+1, signal[signal_list[isignal]].at(i)/sqrt(signal[signal_list[isignal]].at(i)+ background_total));
	}
      }
      else if (mode_==3) {
	if ( background_total >0) 
	  Hcutflow[signal_list[isignal]]->SetBinContent(i+1,signal[signal_list[isignal]].at(i)/sqrt(background_total) );
      }
      
    }
    
    cout<<endl;
  }
  //this is working! need to: combine QCD into one number; 
  //do similar for other backgrounds;  -- done
  //add signal;  --done
  //compute s/root(s+b)
  //etc

  Ccutflow = new TCanvas("Ccutflow");
  //  if (mode_==1)  Ccutflow->SetLogy();
  THStack mystack("mystack","Cut Flow Steps"); //only used for mode 1
  float legx1=0.6,legy1=0.6,legx2=0.9,legy2=0.9;
  if (mode_==2 || mode_==3) {
    legx1=0.1; legx2=0.3;
  }
  TLegend leg(legx1,legy1,legx2,legy2);
  leg.SetFillColor(0);

  TString opt="";
  double highest=0;
  if (mode_==1) {
    mystack.Add( Hcutflow["QCD"]); //add qcd first
    leg.AddEntry(Hcutflow["QCD"],"QCD");
    for (int ib=0; ib<nbackground; ib++)  {
      mystack.Add( Hcutflow[ background_list[ib]]); //then backgrounds
      leg.AddEntry(Hcutflow[ background_list[ib]],background_list[ib]);
    }
    mystack.Add(Hcutflow[mysig]); //finally signal
    leg.AddEntry(Hcutflow[mysig],mysig);
  }
  else {
    for ( std::map<TString, TH1D*>::const_iterator ib=Hcutflow.begin(); ib!=Hcutflow.end(); ++ib) {
      ib->second->Draw(opt);
      opt="SAME";
      double height=ib->second->GetMaximum();
      if (height>highest) highest=height;
      leg.AddEntry(ib->second,ib->first);
    }
  }

  gStyle->SetOptStat(0);
  TAxis* ax=0;
  if (mode_==1) {
    mystack.SetMaximum(750); //hard-coded!
    mystack.Draw("HIST");
    TH1* hax=mystack.GetHistogram();
    hax->SetYTitle("Events (50 pb^{-1})");
    ax = hax->GetXaxis();
  }
  else { //this is tricky
    ax= Hcutflow.begin()->second->GetXaxis();
    Hcutflow.begin()->second->SetMaximum(highest*1.05);
    Hcutflow.begin()->second->SetTitle("");
    Hcutflow.begin()->second->SetYTitle(modeDescription);
  }
  cout<<ax<<endl;
  
  for (unsigned int icut=0; icut<cutnames.size(); icut++) {
    ax->SetBinLabel(icut+1,cutnames.at(icut));
  }
  leg.Draw();

  TString outname = filestub_;
  outname+=".";
  if (mode_==1) {
    outname+=mysig;
    outname+=".";
  }
  outname+=mode_;
  outname+=".eps";
  Ccutflow->SaveAs(outname);

  /*

  | *Step* | *Description* | *MoreMSSM* | *LM3* | *LM9* | *LM13* | *TTbarJets* | *SingleTop <br />(t,tW)* | *WJets* | *ZJets* | *Zinvisible* | *VqqJets* | *QCD* | *S(MMSSM)<br />/B* | *S(LM3)<br />/B* | *S(LM9)<br />/B* | *S(LM13)<br />/B* | *S(MMSSM)<br />/&radic;(B* | *S(LM3)<br />/&radic;(B)* | *S(LM9)<br />/&radic;(B)* | *S(LM13)<br />/&radic;(B)* | *S(MMSSM)<br />/&radic;(S+B)* | *S(LM3)<br />/&radic;(S+B)* | *S(LM9)<br />/&radic;(S+B)* | *S(LM13)<br />/&radic;(S+B)* |
| *Step 0* | Inclusive | 173 | 343.8 | 713.4 | 690 | 16500 | 3100 ± 2 | 2417000 | 280000 | 450000 | 3580 | 717628300 | 1E-4 | 1E-4 | 2E-4 | 2E-4 | 0.10 | 0.19 | 0.40 | 0.39 | 0.10 | 0.19 | 0.40 | 0.39 |

  */

}
