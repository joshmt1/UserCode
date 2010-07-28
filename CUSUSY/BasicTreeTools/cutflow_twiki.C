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

void cutflow_twiki()
{
  //first we need to load each text file (one sample at a time)
  const TString filestub ="cutflow_RA2withBtagging";

  const int mode = 2; //mode 1 is print cut flow table ; mode 2 is print significance table
  assert(mode==1 || mode==2);

  //in principle these are coded in basicLoop.C/h; do the easy thing for now
  std::vector<TString> cutnames;
  cutnames.push_back("Inclusive");
  cutnames.push_back("Trigger");
  cutnames.push_back("PV");
  cutnames.push_back(">= 3 Jets");
  cutnames.push_back("HT Cut");
  cutnames.push_back("MHT Cut");
  cutnames.push_back("Muon veto");
  cutnames.push_back("electron veto");
  cutnames.push_back("Delta Phi");
  cutnames.push_back(">=1 b");
  cutnames.push_back(">=2 b");
  cutnames.push_back(">=3 b");


  //is this really the best way to do this?
  int nqcd = 4;
  char *qcd_list[]={"QCD100","QCD250","QCD500","QCD1000"};
  int nbackground = 6;
  char *background_list[]={"TTbarJets","SingleTop-tChannel","SingleTop-tWChannel","Zinvisible","WJets","ZJets"};
  int nsignal = 16; //oops, where did LM3 go?
  char *signal_list[]={"LM0", "LM1", "LM2", "LM4", "LM5", "LM6","LM7", "LM8","LM9","LM9p", "LM9t175", "LM10", "LM11", "LM12","LM13","mMSSM"};

  //as long as i compile, I can use the most basic stl containers
  //SHIT ... what was I thinking? what i really want is a map of these, indexed by the names.
  //can ROOT handle it? //seems ok....

  std::map<TString, std::vector<float> > qcd;
  std::map<TString, std::vector<float> > qcderr;

  std::map<TString, std::vector<float> > background;
  std::map<TString, std::vector<float> > backgrounderr;

  std::map<TString, std::vector<float> > signal;
  std::map<TString, std::vector<float> > signalerr;

  //lordy, this could really be more generic!
  for (int i=0 ; i<nqcd; i++) {
    TString filename = filestub;
    filename+="."; filename+=qcd_list[i];
    filename += ".dat";

    cout<<"Reading "<<filename<<endl;
    float nevt,nevterr;
    ifstream file(filename.Data());
    if (!file.good()) {cout<<"bad file! "<<filename<<endl; return;}
    while (file>>nevt>>nevterr ) {
      qcd[TString(qcd_list[i])].push_back(nevt);
      qcderr[TString(qcd_list[i])].push_back(nevterr);
    }
    file.close();
  }

  for (int i=0 ; i<nbackground; i++) {
    TString filename = filestub;
    filename+="."; filename+=background_list[i];
    filename += ".dat";

    cout<<"Reading "<<filename<<endl;
    float nevt,nevterr;
    ifstream file(filename.Data());
    if (!file.good()) {cout<<"bad file! "<<filename<<endl; return;}
    while (file>>nevt>>nevterr ) {
      background[TString(background_list[i])].push_back(nevt);
      backgrounderr[TString(background_list[i])].push_back(nevterr);
    }
    file.close();
    //    cout<<  background[TString(background_list[i])].size()<<"\t"<<backgrounderr[TString(background_list[i])].size()<<endl;
  }
  cout<<"--"<<endl;
  cout<<  background.size()<<"\t"<<  backgrounderr.size()<<endl;

  for (int i=0 ; i<nsignal; i++) {
    TString filename = filestub;
    filename+="."; filename+=signal_list[i];
    filename += ".dat";

    cout<<"Reading "<<filename<<endl;
    float nevt,nevterr;
    ifstream file(filename.Data());
    if (!file.good()) {cout<<"bad file! "<<filename<<endl; return;}
    while (file>>nevt>>nevterr ) {
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
  //need to get the cutnames into this script....
  const  TString pm = " +/- ";
  const  TString col = " | ";

  cout<<col<<"n" <<col<<"Cut"<<col;
  if (mode==1) {
    cout<<"QCD"<<col;
    for (int ibackground=0 ; ibackground<nbackground; ibackground++)     cout<< background_list[ibackground]<<col;
  }
  for (int isignal=0 ; isignal<nsignal; isignal++)     cout<<signal_list[isignal]<<col;
  cout<<endl;

  for (int i=0; i<ncuts; i++) {
    cout<<col<< i<<col<<cutnames.at(i)<<col;

    //cout<<setprecision(1);

    float qcd_total=0;
    float qcd_total_err=0;
    //it would be more elegant to iterate over the map...but to hell with elegance
    for (int iqcd=0 ; iqcd<nqcd; iqcd++) {
      //      cout<<     qcd[qcd_list[iqcd]].at(i) << " | ";
      qcd_total += qcd[qcd_list[iqcd]].at(i);
      qcd_total_err += pow(qcderr[qcd_list[iqcd]].at(i),2);
    }
    float background_total=qcd_total;
    float background_total_err=qcd_total_err; //before taking the square root

    qcd_total_err = sqrt(qcd_total_err);
    if (mode==1) cout<< qcd_total << pm << qcd_total_err<<col;

    for (int ibackground=0 ; ibackground<nbackground; ibackground++) {
      if (mode==1)
	cout<< background[background_list[ibackground]].at(i) <<pm<<backgrounderr[background_list[ibackground]].at(i) << col;
      background_total += background[background_list[ibackground]].at(i);
      background_total_err += pow(backgrounderr[background_list[ibackground]].at(i),2);
    }

    //now that we've summed the squares, take the square root
    background_total_err = sqrt(background_total_err);

    for (int isignal=0 ; isignal<nsignal; isignal++) {
      if (mode==1)      cout<<signal[signal_list[isignal]].at(i) <<pm<< signalerr[signal_list[isignal]].at(i) << col;
      else if (mode==2) {
	if (signal[signal_list[isignal]].at(i)+ background_total >0)
	  cout<<signal[signal_list[isignal]].at(i)/sqrt(signal[signal_list[isignal]].at(i)+ background_total)  << col;
	else  cout<<" - "<<endl;
      }
    }
    
    cout<<endl;
  }
  //this is working! need to: combine QCD into one number; 
  //do similar for other backgrounds;  -- done
  //add signal;  --done
  //compute s/root(s+b)
  //etc

  /*

  | *Step* | *Description* | *MoreMSSM* | *LM3* | *LM9* | *LM13* | *TTbarJets* | *SingleTop <br />(t,tW)* | *WJets* | *ZJets* | *Zinvisible* | *VqqJets* | *QCD* | *S(MMSSM)<br />/B* | *S(LM3)<br />/B* | *S(LM9)<br />/B* | *S(LM13)<br />/B* | *S(MMSSM)<br />/&radic;(B* | *S(LM3)<br />/&radic;(B)* | *S(LM9)<br />/&radic;(B)* | *S(LM13)<br />/&radic;(B)* | *S(MMSSM)<br />/&radic;(S+B)* | *S(LM3)<br />/&radic;(S+B)* | *S(LM9)<br />/&radic;(S+B)* | *S(LM13)<br />/&radic;(S+B)* |
| *Step 0* | Inclusive | 173 | 343.8 | 713.4 | 690 | 16500 | 3100 ± 2 | 2417000 | 280000 | 450000 | 3580 | 717628300 | 1E-4 | 1E-4 | 2E-4 | 2E-4 | 0.10 | 0.19 | 0.40 | 0.39 | 0.10 | 0.19 | 0.40 | 0.39 |

  */

}
