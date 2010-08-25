#include <iostream>

#include "TString.h"
#include "TH1.h"
#include "TTree.h"

#include "PlotUtil.h"

using namespace std;

PlotUtil::PlotUtil(HistHolder* hh) :
  hh_(hh),
  debug_(false),
  lastVarname_("MET")
{
}

PlotUtil::~PlotUtil() {}

Double_t PlotUtil::ErrorOnIntegral(const TH1D* h, const Int_t lowbin, Int_t highbin) {
  
  if ( highbin == 0)  highbin = h->GetNbinsX();

  Double_t err=0;

  Double_t thisbin=0;
  for ( int i = lowbin; i<=highbin ; i++) {
    thisbin = h->GetBinError(i);
    err += thisbin*thisbin;
  }

  return sqrt( err );

}

void PlotUtil::createHistos(TString varname,int nbins, float hmin, float hmax) {
  if (varname=="") varname = lastVarname_;

  for ( map<TString,TTree*>::const_iterator isamp=samples_.begin() ; isamp!=samples_.end() ; ++isamp) {
    TString fullname;
    fullname.Form("H%s_%s",varname.Data(),isamp->first.Data());
    hh_->make(fullname.Data(),varname.Data(),nbins,hmin,hmax);
  }
  lastVarname_=varname;
}

void PlotUtil::fillHistos(TString varname,TString cut, TString drawcommand) {
  if (varname=="") varname = lastVarname_;

  if (drawcommand=="") drawcommand=varname;

  for ( map<TString,TTree*>::const_iterator i=samples_.begin() ; i!=samples_.end() ; ++i) {
    i->second->Project(getFullName(varname,i->first),drawcommand,cut);
  }
  lastVarname_=varname;
}

TString PlotUtil::getFullName(TString varname,TString sampleid) {

  TString fullname;
  fullname.Form("H%s_%s",varname.Data(),sampleid.Data());
  return fullname;
}

void PlotUtil::addQCD(TString varname) {
  if (varname=="") varname = lastVarname_;

  const  TString qcdhist="H"+varname+"_qcd";

  TH1D* hbase=hh_->find(qcdhist);
  if (hbase==0) {
    cout<<"Could not find the _qcd histogram!"<<endl;
    return;
  }

  for ( map<TString,TTree*>::const_iterator i=samples_.begin() ; i!=samples_.end() ; ++i) {
    TString fullname=getFullName(varname,i->first);
    if (  fullname!= qcdhist && fullname.Contains("_qcd") ) {
      if (debug_)   cout<<"Adding to base qcd histo: "<<fullname<<endl;
      TH1D* htoadd=hh_->find(fullname);
      hbase->Add(htoadd);
    }
  }
  lastVarname_=varname;
}

void PlotUtil::setLineColors(TString varname) {
  if (varname=="") varname = lastVarname_;

  for ( map<TString,UInt_t>::const_iterator i=sampleColors_.begin() ; i!=sampleColors_.end() ; ++i) {
    TString fullname=getFullName(varname,i->first);
    if (debug_)  {
      cout<<fullname<<endl;
      cout<<hh_->find(fullname)<<endl;
    }
    hh_->find(fullname)->SetLineColor(i->second);
  }
  
  lastVarname_=varname;
}

void PlotUtil::drawPlots(TString varname) {
  if (varname=="") varname = lastVarname_;

  //want to exclude the secondary qcd histos
  const TString qcdhist="H"+varname+"_qcd";

  //find the "tallest" histogram
  TH1D* firsttodraw=0;
  double max=0;
  for ( map<TString,TTree*>::const_iterator i=samples_.begin() ; i!=samples_.end() ; ++i) {
    TString fullname=getFullName(varname,i->first);

    if ( (!fullname.Contains("_qcd")) || (fullname ==qcdhist) ) {
      double mymax=   hh_->find(fullname)->GetMaximum();
      if (mymax>max) {
	max=mymax;
	firsttodraw = hh_->find(fullname);
      }
    }
  }
  
  firsttodraw->Draw();

  //draw the rest
  for ( map<TString,TTree*>::const_iterator i=samples_.begin() ; i!=samples_.end() ; ++i) {
    TString fullname=getFullName(varname,i->first);
    if ( (!fullname.Contains("_qcd")) || (fullname ==qcdhist) ) {
      TH1D* hist= hh_->find(fullname);
      if (hist!= firsttodraw) hist->Draw("SAME");
    }
  }

  lastVarname_=varname;
}

void PlotUtil::fillLegend(TLegend * leg, TString varname) {
  if (varname=="") varname = lastVarname_;

  if (debug_) cout<<"[fillLegend]"<<endl;
  setLineColors(varname);
  if (debug_) cout<<"[done with setLineColors()]"<<endl;

  //want to exclude the secondary qcd histos
  const  TString qcdhist="H"+varname+"_qcd";

  for ( map<TString,TTree*>::const_iterator i=samples_.begin() ; i!=samples_.end() ; ++i) {
    TString fullname=getFullName(varname,i->first);
    if ( (!fullname.Contains("_qcd")) || (fullname ==qcdhist) ) {
      leg->AddEntry(hh_->find(fullname), sampleNames_[i->first]);
    }
  }
  lastVarname_=varname;
}

void PlotUtil::addSample(TString sampleid,TTree* sampletree,TString sampleName,UInt_t kcolor) {

  samples_[sampleid]=sampletree;

  if (sampleName=="") sampleName=sampleid;
  sampleNames_[sampleid]=sampleName;

  sampleColors_[sampleid]=kcolor;

}

TString PlotUtil::fortranize(TString cut) {

  cut.ReplaceAll("==","eq");
  cut.ReplaceAll(">=","gte");
  cut.ReplaceAll("<=","lte");

  cut.ReplaceAll(">","gt");
  cut.ReplaceAll("<","lt");

  return cut;
}
