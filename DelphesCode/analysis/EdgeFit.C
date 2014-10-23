/*
--- to compile ---

.L TSelectorMultiDraw.C+
.L CrossSectionTable.cxx+
.L ConfigurationDescriptions.cxx+

.L RooEdge.cxx+
.L RooEdgeFlavSym.cxx+
.L RooTopPairProductionSpline.cxx+
.L RooCruijff.cxx+
.L RooBifurCB.cxx+
.L RooFallingSpectrum.cxx+

.L EdgeFit.C+

*/
#include "drawDelphesBase.C"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooAbsData.h"
//#include "RooGenericPdf.h"
//#include "RooProdPdf.h"

#include "TFitResultPtr.h"
#include "TFitResult.h"

#include "RooSimultaneous.h"
#include "RooCategory.h"

#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooVoigtian.h"
#include "RooGaussian.h"

#include "RooEdge.h"
#include "RooEdgeFlavSym.h"
#include "RooTopPairProductionSpline.h"
#include "RooCBShape.h"
#include "RooCruijff.h"
#include "RooBifurCB.h"
#include "RooFallingSpectrum.h"
#include "RooLognormal.h"

#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TCut.h"
#include "TRandom.h"
#include "TRegexp.h"
#include "TObjArray.h"

#include "TSystemFile.h"
#include "TSystemDirectory.h"

#include <utility>

Long_t extractSeedFromFilename(TString filename) {
  
  //now find the toy index (seed)
  TRegexp re2("_[0-9]+.root");
  int startplace=filename.Index(re2);
  TString theseedstr=filename(startplace+1,filename.Length()-startplace-6);
  Long_t seed=theseedstr.Atoi();
  
  return seed;
}

void fitSignal(TString somecuts="(1)",TString cuts_description="") {

  //  initSamples("signal");
  //TChain* signal = getTree("susyhit_Scenario1_v02");

  TChain * signal = new TChain("simpleTree");
  signal->Add("/cu3/joshmt/Upgrade/PhaseII_Configuration4v2_140PileUp/v10/reducedTree.def.naturalModel1.root");

  TCut edge_signal = "isSF_loose==1&&leptonsMatchChi2ToChi1_loose==1";
  TCut ee = "abs(leptonFlavor1_loose)==11";
  TCut mm = "abs(leptonFlavor1_loose)==13";
  TCut othercuts = somecuts.Data();

  TTree* signal_true_ee = signal->CopyTree(TCut(ee && edge_signal && othercuts).GetTitle());
  TTree* signal_true_mm = signal->CopyTree(TCut(mm && edge_signal && othercuts).GetTitle());

  RooRealVar mll("mll_loose","mll",20,100);
  //maybe i could get the same results with the cutVar argument here. whatever
  RooDataSet r_signal_true_ee("r_signal_true_ee","true edge ee",signal_true_ee,RooArgSet(mll));
  RooDataSet r_signal_true_mm("r_signal_true_mm","true edge mm",signal_true_mm,RooArgSet(mll));

  //strangely, the fit converges better when the edge needs to drift up a bit rather than starting bang on
  RooRealVar edge("edge","edge",67,0,100);
  RooRealVar sigma_ee("sigma_ee","sigma_ee",2,0.01,10);
  RooRealVar sigma_mm("sigma_mm","sigma_mm",2,0.01,10);

  //  sigma_ee.setConstant();
  //  sigma_mm.setConstant();

  RooEdge edge_ee("edge_ee","edge pdf ee",mll,edge,sigma_ee);
  RooFitResult* eefit =   edge_ee.fitTo(r_signal_true_ee,RooFit::Save(1));
  eefit->Print();
  RooPlot* eeframe = mll.frame() ;
  r_signal_true_ee.plotOn(eeframe) ;
  edge_ee.plotOn(eeframe) ;

  edge.setVal(68);
  RooEdge edge_mm("edge_mm","edge pdf mm",mll,edge,sigma_mm);
  RooFitResult* mmfit =   edge_mm.fitTo(r_signal_true_mm,RooFit::Save(1));
  mmfit->Print();
  RooPlot* mmframe = mll.frame() ;
  r_signal_true_mm.plotOn(mmframe) ;
  edge_mm.plotOn(mmframe) ;

  TString filename = "EdgeFit_signalonly";
  if (cuts_description!="") {filename += "_"; filename += cuts_description;}
  filename+=".root";
  TFile fout(filename,"recreate");
  eeframe->SetName("edge_ee_plot");  eeframe->Write();
  mmframe->SetName("edge_mm_plot");  mmframe->Write();
  fout.Close();


}

void fitOF(TString option="datastats") { //option is datastats or mcstats

  RooDataSet * rds_of=0;

  RooRealVar * mll_loose=0;
  RooRealVar * weight=0;
  bool useweighting=false;
  //use combinefs or combinefssig
  if (option.Contains("mcstats")) {
    useweighting=true;
    initSamples("combinefssig skimmed"); //proxy for only FS events, although not perfect
    
    TCut dileptons_loose="mll_loose>20 && mll_loose<=150";
    TCut of_loose = "isSF_loose==0";
    TCut selection = dileptons_loose && of_loose && TCut("njets40>=6") && TCut("MET>450") && TCut("nbjets40medium>=3")&&TCut("HT>1250");
    
    TTree* of = getTree("tt-4p-0-600-v1510_14TEV")->CopyTree(selection);
    
    mll_loose = new RooRealVar("mll_loose","mll_loose",20,150);
    weight = new RooRealVar("weight","weight",0,100);
    rds_of = new RooDataSet("rds_of","of DY",of,RooArgSet(*mll_loose,*weight),0,"weight");
  }
  else if (option.Contains("datastats")) {
    mll_loose = new RooRealVar("mll","mll_loose",20,150); //need name 'mll' for this dataset
    TFile fin("DelphesEdge/datasets.root");
    RooDataSet* data_of_all = (RooDataSet*) fin.Get("data_of");
    rds_of = (RooDataSet*) data_of_all->reduce("mll<=150");
    fin.Close();
  }
  else assert(0);


  
  //luke's pdf
  RooRealVar turn_on_mass("turn_on_mass","turn on",0);
  RooRealVar dummy1("dummy1","",0);
  RooRealVar dummy4("dummy4","",0);
  RooRealVar dummy5("dummy5","",0);
  RooRealVar dummy6("dummy6","",0);
  RooRealVar dummy7("dummy7","",0);
  RooRealVar peakmass("peakmass","peak",50,20,90);
  RooRealVar lognsigma("lognsigma","width",1.9847,0.0001,100);
  RooRealVar gtheta("gtheta","theta",80,0.01,200);
  RooRealVar fracvar("fracvar","f",-0.5,-10,10);
  RooFallingSpectrum of_pdf("of_pdf","falling spectrum",*mll_loose,turn_on_mass,dummy1,peakmass,lognsigma,gtheta,fracvar,dummy4,dummy5,dummy6,dummy7);

  //  RooRealVar turn_on_alpha("turn_on_alpha","turn_on_alpha",10,1,30);
  //  RooRealVar turn_on_beta("turn_on_beta","turn_on_beta",0.98,0.4,0.9999999);
  //  RooGenericPdf turn_on_pdf("turn_on_pdf","turn on","1.0-exp(-1.0*pow(mll_loose/turn_on_alpha,turn_on_beta))",RooArgList(*mll_loose,turn_on_alpha,turn_on_beta));
  //  RooProdPdf of_pdf("of_pdf","of_pdf",of_pdf_luke,turn_on_pdf);

  //  lognsigma.setConstant(true);

  //lognormal
  //RooLognormal of_pdf("of_pdf","log normal",*mll_loose,peakmass,lognsigma);

  //for multi-piece pdf
  RooRealVar c("c","b1Central",1000,-8000,1000000);
  RooRealVar d("d","b2Central",120,-800,800);
  RooRealVar f("f","b4Central",0.002,0.0000001,0.1);
  RooRealVar m1("m1","m1Central",35,22,50);
  RooRealVar m2("m2","m2Central",75,65,160);

  //cruijff
  RooRealVar cjmean("cjmean","",50,20,100);
  RooRealVar cjsl("cjsl","",20,0.1,100);
  RooRealVar cjsr("cjsr","",100,0.1,999);
  RooRealVar cjal("cjal","",0.1,-10,10);
  RooRealVar cjar("cjar","",0.1,-10,10);

  
  //crystal ball
  RooRealVar cb0("cb0","mean",60,30,100);
  RooRealVar cbsigma("cbsigma","sigma",30,1,200);
  RooRealVar cbalpha("cbalpha","alpha",-1,-5,5);
  RooRealVar cbn("cbn","n",1,0,9);
  RooCBShape cb_pdf("cb_pdf","of pdf",*mll_loose,cb0,cbsigma,cbalpha,cbn);

  RooRealVar fcb("fcb","cb const",7.7544e-01,0,1);

  //bifur cb -- don't share params with cb for now
  RooRealVar bcb0("bcb0","mean",5.7815e+01,30,100);
  RooRealVar bcbsigmal("cbsigmal","sigma left",3.3278e+01,1,200);
  RooRealVar bcbsigmar("cbsigmar","sigma right",3.3278e+01,1,800);
  RooRealVar bcbalpha("bcbalpha","alpha", -9.1561e-01,-5,5);
  RooRealVar bcbn("bcbn","n",7.7442e-01,-2,9);
  RooBifurCB bcb_pdf("bcb_pdf","of pdf",*mll_loose,bcb0,bcbsigmal,bcbsigmar,bcbalpha,bcbn);

  //for simpler pdf
  RooRealVar alpha("alpha","alpha",8.2828e-01,0.01,10);
  RooRealVar beta("beta","beta", 1.0366e-02,-5,5);
  //'simple' pdf
  RooEdgeFlavSym fs_pdf("fs_pdf","fs pdf",*mll_loose,alpha,beta);

  //  RooAddPdf of_pdf("of_pdf","",bcb_pdf,fs_pdf,fcb); //nice but not well-behaved with data stats
  //  RooCruijff of_pdf("of_pdf","",mll_loose,cjmean,cjsl,cjsr,cjal,cjar);
    //multi-piece pdf
  //  RooTopPairProductionSpline of_pdf("of_pdf","of pdf official",mll_loose,c,d,f,m1,m2);

  RooFitResult* offit = of_pdf.fitTo(*rds_of,RooFit::Save(1),RooFit::SumW2Error(useweighting));
  offit->SetName("offitResult");

  RooPlot* offrame = mll_loose->frame() ;
  rds_of->plotOn(offrame) ;
  of_pdf.plotOn(offrame) ;

  RooPlot* ofzoom = mll_loose->frame(20,120,50) ;
  rds_of->plotOn(ofzoom) ;
  of_pdf.plotOn(ofzoom) ;

  cout<<"chi2 values = "<<offrame->chiSquare()<<" zoom = "<<ofzoom->chiSquare()<<endl;

  TFile fout("EdgeFit_OF.root","recreate");
  offrame->SetName("edge_of_plot");  offrame->Write();
  ofzoom->SetName("edge_of_plot_zoom");  ofzoom->Write();
  offit->Write();
  fout.Close();

}

void fitDY(bool drawPictures=true) {

  //  TChain * dy = new TChain("simpleTree");
  //  dy->Add("/cu3/joshmt/Upgrade/PhaseII_Configuration4v2_140PileUp/v02/reducedTree.def.B-4p-0-1-v1510_14TEV.root");
  //  dy->Add("trees/reducedTree.def.Bj-*.root");
  //dy->Add("skimmed_mll_ht/reducedTree.def.*.root");
  initSamples("combinenonfs skimmed"); //proxy for only non-FS events, although not perfect
  //  initSamples("nm1"); //with my selection, only Z is in signal
  TCut dileptons_loose="mll_loose>20";
  TCut sf_loose = "isSF_loose==1";
  TCut selection = dileptons_loose && sf_loose && TCut("njets40>=6") && TCut("MET>450") && TCut("nbjets40medium>=1")&&TCut("HT>1250");

  //  TTree* sf = getTree("LL-4p-0-100-v1510_14TEV")->CopyTree(selection.GetTitle());
  TTree* sf = getTree("naturalModel1")->CopyTree(selection.GetTitle());

  RooRealVar mll_loose("mll_loose","mll_loose",20,300);
  RooRealVar weight("weight","weight",0,100);
  RooDataSet rds_sf("rds_sf","sf DY",sf,RooArgSet(mll_loose,weight),0,"weight");

  RooRealVar cexp("cexp","exponential const",0.1,-10,10); //for the RooExponential
  RooExponential pdf_exp("pdf_exp","exp pdf",mll_loose,cexp);

  //  RooRealVar cpol2("cpol2","pol2 coef",0,-100,100);
  //  RooPolynomial pdf_pol("pdf_pol","pol pdf",mll_loose,RooArgList(cexp,cpol2));

  //now the Voigtian
  RooRealVar zmass("zmass","z mass",91,80,100);
  RooRealVar zwidth("zwidth","z width", 2.5,0.1,10);
  RooRealVar zsigma("zsigma","z resolution",0.5,0.00001,20);
  RooVoigtian pdf_peak("pdf_peak","z peak pdf",mll_loose,zmass,zwidth,zsigma);

  RooRealVar fexp("fexp","non-bw fraction",0.1,0,1);
  RooAddPdf pdf_dy("pdf_dy","total dy pdf",pdf_exp,pdf_peak,fexp);
  //RooAddPdf pdf_dy("pdf_dy","total dy pdf",pdf_pol,pdf_peak,fexp);

  RooFitResult* fitresult =  pdf_dy.fitTo(rds_sf,RooFit::Save(1),RooFit::SumW2Error(kTRUE));
  fitresult->Print(); fitresult->SetName("Zfit");
  TFile fout("EdgeFit_DY.root","RECREATE");
  fitresult->Write();
  if (drawPictures) {
    RooPlot* dyframe = mll_loose.frame();
    RooPlot* dyframe2 = mll_loose.frame(80,100,30);
    rds_sf.plotOn(dyframe,RooFit::DataError(RooAbsData::SumW2));
    pdf_dy.plotOn(dyframe);
    rds_sf.plotOn(dyframe2,RooFit::DataError(RooAbsData::SumW2));
    pdf_dy.plotOn(dyframe2);
    TCanvas cfit("cfit","",800,800);
    cfit.Divide(2,2);
    cfit.cd(1); dyframe->Draw();
    cfit.cd(2); dyframe2->Draw();
    cfit.cd(3); dyframe->Draw();
    cfit.cd(4); dyframe2->Draw();
    cfit.GetPad(3)->SetLogy();
    cfit.GetPad(4)->SetLogy();
    cfit.Write();    
    dyframe->SetName("edge_dy_plot");  dyframe->Write();
  }
  fout.Close();

}

void sanityChecks() {
  initSamples("combinefssig skimmed"); //proxy for only FS events, although not perfect

  TCut dileptons_loose="mll_loose>20 && mll_loose<=150";
  TCut of_loose = "isSF_loose==0";
  TCut selection = dileptons_loose && of_loose && TCut("njets40eta3p0>=4") && TCut("MET>400") && TCut("nbjets40tight>=1")&&TCut("HT>1500");

  TTree* of = getTree("tt-4p-0-600-v1510_14TEV")->CopyTree(selection);

  TH1D h_direct("h_direct","direct",40,20,150);
  h_direct.Sumw2();
  of->Draw("mll_loose>>h_direct","weight","goff");

  TH1D h_rds("h_rds","rds",40,20,150);
  h_rds.Sumw2();

  // ====== now the roodataset part
  TFile fin("DelphesEdge/datasets_3000_102766.root");
  RooDataSet* data_dy = (RooDataSet*) fin.Get("data_dy");
  RooDataSet* data_of_all = (RooDataSet*) fin.Get("data_of");
  RooDataSet* data_sf_ee = (RooDataSet*) fin.Get("data_sf_ee");
  RooDataSet* data_sf_mm = (RooDataSet*) fin.Get("data_sf_mm");

  //for now, combine ee and mm!
  data_sf_ee->append(*data_sf_mm); 

  //our pdf is validated to 150
  RooDataSet* data_sf = (RooDataSet*) data_sf_ee->reduce("mll<=150");
  RooDataSet* data_of = (RooDataSet*) data_of_all->reduce("mll<=150");

  for (int ii=0; ii<data_of->numEntries(); ii++) {
    h_rds.Fill( ((RooRealVar*)data_of_all->get(ii)->find("mll"))->getVal());
  }

  cout<<"KS = "<<  h_rds.KolmogorovTest(&h_direct)<<endl;

  h_rds.Scale(1.0/h_rds.Integral());
  h_direct.Scale(1.0/h_direct.Integral());

  cout<<"KS = "<< h_rds.KolmogorovTest(&h_direct)<<endl;

  TFile fout("sanity.root","RECREATE");
  h_rds.Write();
  h_direct.Write();
  fout.Close();


}

void generateDatasets( const float lumi=3000, UInt_t seed = 123456, const bool backgroundonly=false ) {

  float scale = lumi/3000;

  gRandom->SetSeed(seed);

  TString thepath = "DelphesEdge/";
  if ( TString(gSystem->Getenv("TOYPATH"))!="") thepath = TString(gSystem->Getenv("TOYPATH"));
  TFile fin(thepath+"templates.root");

  //need a:
  //DY dataset
  //OF dataset
  //SF dataset including signal

  TString histostub = backgroundonly ? "smsusynoedge" : "smsusy";

  TH1D* dy_template = backgroundonly ? 0 : (TH1D*) fin.Get(histostub+"_mll_DY");
  TH1D* of_template = (TH1D*) fin.Get(histostub+"_mll_OF");
  TH1D* sfee_template = (TH1D*) fin.Get(histostub+"_mll_ee");
  TH1D* sfmm_template = (TH1D*) fin.Get(histostub+"_mll_mm");


  int n_dy = backgroundonly ? 0 : gRandom->Poisson(dy_template->Integral()*scale);
  int n_of = gRandom->Poisson(of_template->Integral()*scale);
  int n_sf_ee = gRandom->Poisson(sfee_template->Integral()*scale);
  int n_sf_mm = gRandom->Poisson(sfmm_template->Integral()*scale);

  RooRealVar mll("mll","mll",0,300);
  RooDataSet data_dy("data_dy","dy data",RooArgSet(mll));
  for (int i=0;i<n_dy;i++) {
    mll.setVal( dy_template->GetRandom());
    data_dy.add(RooArgSet(mll));
  }
  RooDataSet data_of("data_of","of data",RooArgSet(mll));
  for (int i=0;i<n_of;i++) {
    mll.setVal( of_template->GetRandom());
    data_of.add(RooArgSet(mll));
  }
  RooDataSet data_sf_ee("data_sf_ee","sfee data",RooArgSet(mll));
  for (int i=0;i<n_sf_ee;i++) {
    mll.setVal( sfee_template->GetRandom());
    data_sf_ee.add(RooArgSet(mll));
  }
  RooDataSet data_sf_mm("data_sf_mm","sfmm data",RooArgSet(mll));
  for (int i=0;i<n_sf_mm;i++) {
    mll.setVal( sfmm_template->GetRandom());
    data_sf_mm.add(RooArgSet(mll));
  }

  if (backgroundonly) {
    thepath += "NoEdge/";
    thepath += (Long_t) seed % 1000;
    thepath+="/";
  }
  else {
    thepath += (Long_t) TMath::Nint(lumi);
    thepath+="/";
  }
  gSystem->mkdir(thepath.Data());
  TString fname;
  fname.Form("%sdatasets_%d_%d.root",thepath.Data(), TMath::Nint(lumi),seed);
  TFile fout(fname,"recreate");
  data_dy.Write();
  data_of.Write();
  data_sf_ee.Write();
  data_sf_mm.Write();
  fout.Close();
  fin.Close();

}

void bulkGenerate_forbinnedfits() {

  const  int ndatasetsperpoint=1000;

  UInt_t aseed = 100001;
  int lumi[]={300,500,1000,1500,2000,2500,3000};
  int nlumi = 7;
  for ( int i=0; i<nlumi; i++) {
    for (int jj=0;jj<ndatasetsperpoint; jj++) {
      generateDatasets(lumi[i],aseed);
      aseed++;
    }
  }
  //  generateDatasets(30000,200001);

}

void bulkGenerate_backgroundonly() {

  const  int ndatasetsperpoint=50000;

  UInt_t aseed = 200001;
  int lumi[]={3000};
  int nlumi = 1;
  for ( int i=0; i<nlumi; i++) {
    for (int jj=0;jj<ndatasetsperpoint; jj++) {
      generateDatasets(lumi[i],aseed,true); //last argument = background only
      aseed++;
    }
  }

}

void full_fit(float lumi,TString inputfile,float constraint_syst=-1,bool background_only=false,float fixed_edge_mass=-1) {

  cout<<" ~~~~~~~ full fit begin ~~~~~~~"<<endl;

  //fixed_edge_mass >0 means fix the edge mass to that position in the fit

  Long_t  thisseed = extractSeedFromFilename(inputfile);
  gRandom->SetSeed(thisseed);

  cout<<" constraint term set to "<<constraint_syst<<endl;

  //input file should point to an input dataset
  //e.g. DelphesEdge/datasets_3000_108980.root

  //load datasets
  TFile fin(inputfile);
  RooDataSet* data_dy = (RooDataSet*) fin.Get("data_dy");
  RooDataSet* data_of_all = (RooDataSet*) fin.Get("data_of");
  RooDataSet* data_sf_ee = (RooDataSet*) fin.Get("data_sf_ee");
  RooDataSet* data_sf_mm = (RooDataSet*) fin.Get("data_sf_mm");

  //for now, combine ee and mm!
  data_sf_ee->append(*data_sf_mm); 

  //our pdf is validated to 150
  RooDataSet* data_sf = (RooDataSet*) data_sf_ee->reduce("mll<=150");
  RooDataSet* data_of = (RooDataSet*) data_of_all->reduce("mll<=150");

  RooRealVar* mll = (RooRealVar*) data_dy->get()->find("mll");
  mll->setMin(20);
  mll->setMax(150); //keep in sync with above

  //fit for Z peak
  RooRealVar cexp("cexp","exponential const",0.1,-10,10); //for the RooExponential
  RooExponential pdf_exp("pdf_exp","exp pdf",*mll,cexp);
  //now the Voigtian
  RooRealVar zmass("zmass","z mass",91,80,100);
  RooRealVar zwidth("zwidth","z width", 2.5,0.1,10);
  RooRealVar zsigma("zsigma","z resolution",0.5,0.00001,20);
  RooVoigtian pdf_peak("pdf_peak","z peak pdf",*mll,zmass,zwidth,zsigma,kTRUE);

  RooRealVar fexp("fexp","exp fraction",0.1,0,1);
  RooAddPdf pdf_dy("pdf_dy","total dy pdf",pdf_exp,pdf_peak,fexp);

  //  RooFitResult* fit_dy=  pdf_dy.fitTo(*data_dy,RooFit::Save(1));
  TFile fdyfit("EdgeFit_DY.root");
  RooFitResult * Zfit = (RooFitResult*) fdyfit.Get("Zfit");
  RooArgList dyfitargs=Zfit->floatParsFinal();
  cexp.setVal(((RooRealVar*)dyfitargs.find("cexp"))->getVal());
  fexp.setVal(((RooRealVar*)dyfitargs.find("fexp"))->getVal());
  zmass.setVal(((RooRealVar*)dyfitargs.find("zmass"))->getVal());
  zsigma.setVal(((RooRealVar*)dyfitargs.find("zsigma"))->getVal());
  zwidth.setVal(((RooRealVar*)dyfitargs.find("zwidth"))->getVal());
  //these fit parameters are then fixed
  cexp.setConstant();
  zmass.setConstant();
  zwidth.setConstant();
  zsigma.setConstant();
  fexp.setConstant();

  //FS PDF 

  //luke's pdf
  RooRealVar turn_on_mass("turn_on_mass","turn on",0); //fix
  RooRealVar dummy1("dummy1","",0);
  RooRealVar dummy4("dummy4","",0);
  RooRealVar dummy5("dummy5","",0);
  RooRealVar dummy6("dummy6","",0);
  RooRealVar dummy7("dummy7","",0);
  RooRealVar peakmass("peakmass","peak",50,20,90);
  RooRealVar lognsigma("lognsigma","width",1.9847); //fix
  RooRealVar gtheta("gtheta","theta",80,0,1e3);
  RooRealVar fracvar("fracvar","f",-0.5,-10,10);
  RooFallingSpectrum pdf_fs("pdf_fs","falling spectrum",*mll,
			    turn_on_mass,dummy1,peakmass,lognsigma,gtheta,fracvar,
			    dummy4,dummy5,dummy6,dummy7);
  /*
  //bifur cb
  RooRealVar bcb0("bcb0","mean",5.7815e+01,30,100);
  RooRealVar bcbsigmal("cbsigmal","sigma left",3.3278e+01,1,200);
  RooRealVar bcbsigmar("cbsigmar","sigma right",3.3278e+01,1,800);
  RooRealVar bcbalpha("bcbalpha","alpha", -9.1561e-01,-5,5);
  RooRealVar bcbn("bcbn","n",7.7442e-01,-2,9);
  RooBifurCB bcb_pdf("bcb_pdf","of pdf",*mll,bcb0,bcbsigmal,bcbsigmar,bcbalpha,bcbn);
  //for simpler pdf
  RooRealVar alpha("alpha","alpha",8.2828e-01,0.01,10);
  RooRealVar beta("beta","beta", 1.0366e-02,-5,5);
  //'simple' pdf
  RooEdgeFlavSym fs_pdf("fs_pdf","fs pdf",*mll,alpha,beta);
  */
  //RooRealVar fcb("fcb","cb const",7.7544e-01,0,1);
  //RooAddPdf pdf_fs("pdf_fs","",bcb_pdf,fs_pdf,fcb);
  //OF PDF
  //  RooEdgeFlavSym pdf_fs_of("pdf_fs_of","FS pdf",*mll,alpha,beta);

  //create combined of/sf dataset
  RooCategory sample("sample","sample");
  sample.defineType("SF");
  sample.defineType("OF");
  RooDataSet data("data","sf+of",RooArgSet(*mll),RooFit::Index(sample),RooFit::Import("SF",*data_sf),RooFit::Import("OF",*data_of));

  //fit to OF sample 
  //this is a 'prefit' to get the parameters close to the best values
  RooFitResult* fit_of = pdf_fs.fitTo(*data_of,RooFit::Save(1));
  fit_of->Print();

  RooPlot* ofplot =  mll->frame();
  data_of->plotOn(ofplot);
  pdf_fs.plotOn(ofplot);
  ofplot->SetName("ofplot");

  //fix the parameters
  //alpha.setConstant();
  //beta.setConstant();

  //signal
  /*
    explanation: if we're fitting for a signal we need to start near the true edge or we end up in local minima
    if we're not, then we want to just start randomly
   */
  float starting_mass=67; float upper_bound=100; float lower_bound=25;
  if (background_only) starting_mass = gRandom->Uniform(lower_bound+5,upper_bound-5);

  RooRealVar m_edge("m_edge","edge mass",starting_mass,lower_bound,upper_bound);
  if (fixed_edge_mass>0) {
    m_edge.setVal(fixed_edge_mass);
    m_edge.setConstant();
  }
  cout<<"Starting mass for edge fit = "<<m_edge.getVal()<<endl;

  RooRealVar edge_sigma("edge_sigma","edge resolution",2); //used to be 1.3385
  RooEdge pdf_edge("pdf_edge","edge pdf",*mll,m_edge,edge_sigma);
  edge_sigma.setConstant();

  float scale = lumi/3000;

  //yields
  float nominal_yield = background_only ? 0 : 170; //yield in 3000 fb-1
  RooRealVar n_edge("n_edge","edge yield",nominal_yield*scale,-80,500*scale);
  RooRealVar n_fs("n_fs","FS yield",data_of->numEntries(),-15,data_of->numEntries()*3);
  RooRealVar n_dy("n_dy","Z yield",30*scale,-80,400*scale);
  n_fs.setConstant(kFALSE); //not fixed

  //constrain the fs yield to OF yield +/- error
  RooRealVar n_of_events("n_of_events","",data_of->numEntries());
  float e_val = constraint_syst>0 ? constraint_syst : 0.1;
  RooRealVar e_of_events("e_of_events","",e_val*data_of->numEntries() );
  RooGaussian fs_constraint("fs_constraint","fs constraint",n_fs,n_of_events,e_of_events);

  RooArgList pdfs(pdf_fs,pdf_dy,pdf_edge);
  RooArgList yields(n_fs,n_dy,n_edge);

  RooAddPdf pdf_total("pdf_total","total pdf",pdfs,yields);

  //hack!
  RooRealVar n_of("n_of","OF yield",data_of->numEntries(),0,data_of->numEntries()*2); 
  RooAddPdf pdf_total_of("pdf_total_of","",RooArgList(pdf_fs),RooArgList(n_of)); //can fix yields to each other by using n_fs

  RooSimultaneous pdf_sim("pdf_sim","",sample);
  pdf_sim.addPdf(pdf_total,"SF");
  pdf_sim.addPdf(pdf_total_of,"OF");

  RooFitResult * rfr=0;
  //RooFitResult* rfr=  pdf_total.fitTo(*data_sf,RooFit::Save(1)); //no simultaneous
  if (constraint_syst<0) rfr = pdf_sim.fitTo(data,RooFit::Save(1));  // no constraint
  else rfr = pdf_sim.fitTo(data,RooFit::Save(1),RooFit::ExternalConstraints(fs_constraint)); //with constraint

  rfr->Print();

  int nbins_plot=26;
  RooPlot* fullfitframe = mll->frame(RooFit::Bins(nbins_plot));
  data.plotOn(fullfitframe,RooFit::Cut("sample==sample::SF"));
  pdf_sim.plotOn(fullfitframe,RooFit::Slice(sample,"SF"),RooFit::ProjWData(sample,data),RooFit::Components(RooArgSet(pdf_fs,pdf_dy)),RooFit::LineStyle(kDotted),RooFit::LineColor(kMagenta));
  pdf_sim.plotOn(fullfitframe,RooFit::Slice(sample,"SF"),RooFit::ProjWData(sample,data));
  pdf_sim.plotOn(fullfitframe,RooFit::Slice(sample,"SF"),RooFit::ProjWData(sample,data),RooFit::Components(pdf_fs),RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen+3));
  fullfitframe->SetName("fullfitframe");

  RooPlot* fullfitframe_of = mll->frame(RooFit::Bins(nbins_plot));
  data.plotOn(fullfitframe_of,RooFit::Cut("sample==sample::OF"));
  pdf_sim.plotOn(fullfitframe_of,RooFit::Slice(sample,"OF"),RooFit::ProjWData(sample,data));
  fullfitframe_of->SetName("fullfitframe_of");

  // ===== now that fit is done!
  // make a fit to the no-signal hypothesis!
  RooFitResult* rfr0=0;
  //could disable this fit in fixed edge mass case (fixed_edge_mass>0) because it is redundant,
  //but for simplicity leave it "on" for now
  n_edge.setVal(0);
  n_edge.setConstant(true);
  m_edge.setConstant(true);
  if (constraint_syst<0) rfr0 = pdf_sim.fitTo(data,RooFit::Save(1));  // no constraint
  else rfr0 = pdf_sim.fitTo(data,RooFit::Save(1),RooFit::ExternalConstraints(fs_constraint)); //with constraint
  
  cout<<"2xDeltaNLL = "<<2 * (rfr0->minNll() - rfr->minNll())<<endl;
  TString outputfile =inputfile;
  TString replacestring="unbinned";
  if (constraint_syst<0) {
    if (fixed_edge_mass>0) {
      replacestring +="Fixedmass";
      replacestring += (Long_t) fixed_edge_mass;
    }
    replacestring+= "fit_" ;
  }
  else {
    TString fm="";
    if (fixed_edge_mass>0) {
      fm="Fixedmass";
      fm+=(Long_t) fixed_edge_mass;
    }
    replacestring.Form("unbinned%sWithConstraint%.2f_",fm.Data(),constraint_syst);
    replacestring= jmt::fortranize(replacestring);
  }
  outputfile.ReplaceAll("datasets_",replacestring);
  TFile fout(outputfile,"recreate");
  fullfitframe->Write();
  fullfitframe_of->Write();
  ofplot->Write();
  rfr->SetName("fitresult");
  rfr->Write();
  rfr0->SetName("fitresult_nosignal");
  rfr0->Write();
  fout.Close();

  fin.Close();
  cout<<" ~~~ end"<<endl;
}

std::pair<RooFitResult*,RooPlot*> do_binned_fit(const TH1D & subtracted,float lumi,bool secondedge=false) {
  float scale = lumi/3000;

  //construct PDF of N_sig*P(Edge) + N_Z*P(Z)
  //binned fit from 20<mll<100
  float upperlimit = secondedge ? 130 : 100;
  RooRealVar mll("mll","mll",20,upperlimit);
  RooDataHist data("data","sub data",RooArgList(mll),&subtracted);

  //Z PDF
  RooRealVar cexp("cexp","exponential const",0.01,-10,10); //for the RooExponential
  RooExponential pdf_exp("pdf_exp","exp pdf",mll,cexp);
  //now the Voigtian
  RooRealVar zmass("zmass","z mass",91,80,100);
  RooRealVar zwidth("zwidth","z width", 2.5,0.1,10);
  RooRealVar zsigma("zsigma","z resolution",1e-2,0.00001,20);
  RooVoigtian pdf_peak("pdf_peak","z peak pdf",mll,zmass,zwidth,zsigma,kTRUE);

  RooRealVar fexp("fexp","exp fraction",1e-2);
  RooAddPdf pdf_dy("pdf_dy","total dy pdf",pdf_exp,pdf_peak,fexp);

  TFile fdyfit("EdgeFit_DY.root");
  RooFitResult * Zfit = (RooFitResult*) fdyfit.Get("Zfit");
  RooArgList dyfitargs=Zfit->floatParsFinal();
  cexp.setVal(((RooRealVar*)dyfitargs.find("cexp"))->getVal());
  fexp.setVal(((RooRealVar*)dyfitargs.find("fexp"))->getVal());
  zmass.setVal(((RooRealVar*)dyfitargs.find("zmass"))->getVal());
  zsigma.setVal(((RooRealVar*)dyfitargs.find("zsigma"))->getVal());
  zwidth.setVal(((RooRealVar*)dyfitargs.find("zwidth"))->getVal());
  //these fit parameters are then fixed
  cexp.setConstant();
  zmass.setConstant();
  zwidth.setConstant();
  zsigma.setConstant();
  fexp.setConstant();

  fdyfit.Close();

  //signal
  float startingpoint = 67;//gRandom->Uniform(25,80);
  RooRealVar m_edge("m_edge","edge mass",startingpoint,20,100);
  RooRealVar edge_sigma("edge_sigma","edge resolution",2,0.01,9);
  RooEdge pdf_edge("pdf_edge","edge pdf",mll,m_edge,edge_sigma);
  //m_edge.setConstant();
  edge_sigma.setConstant();

  RooRealVar m_edge2("m_edge2","2nd edge mass",110,100,130);
  RooEdge pdf_edge2("pdf_edge2","2nd edge pdf",mll,m_edge2,edge_sigma);

  //yields
  RooRealVar n_edge("n_edge","edge yield",100*scale,-40,500*scale);
  RooRealVar n_edge2("n_edge2","2nd edge yield",50*scale,-40,200*scale);
  RooRealVar n_dy("n_dy","Z yield",50*scale,-10,1000*scale);

  RooArgList pdfs(pdf_dy,pdf_edge);
  RooArgList yields(n_dy,n_edge);
  if (secondedge) {
    pdfs.add(pdf_edge2);
    yields.add(n_edge2);
  }
  RooAddPdf pdf_total("pdf_total","total pdf",pdfs,yields);
  RooFitResult* rfr=  pdf_total.fitTo(data,RooFit::Save(1),RooFit::SumW2Error(kTRUE));

  //  cout<<"Pre-fit edge mass = "<<startingpoint<<endl;
  cout<<"NLL = "<<rfr->minNll()<<endl;
  cout<<"Edge mass = "<<m_edge.getVal()<<endl;
  cout<<"Edge yield = "<<n_edge.getVal()<<endl;
  if (secondedge)  cout<<"Edge mass  2 = "<<m_edge2.getVal()<<endl;
  if (secondedge)  cout<<"Edge yield 2 = "<<n_edge2.getVal()<<endl;

  RooPlot* fullfitframe = mll.frame();
  data.plotOn(fullfitframe,RooFit::DataError(RooAbsData::SumW2));
  pdf_total.plotOn(fullfitframe);

  //redo the fit with no-edge hypothesis
  //  n_edge.setVal(0);
  //  n_edge.setConstant();
  //TO MAKE THIS WORK -- at a minimum need to fix the edge mass.
  //AT that point there is barely a fit!
  //  RooFitResult* rfr0=  pdf_total.fitTo(data,RooFit::Save(1),RooFit::SumW2Error(kTRUE));
  //  cout<<"Nominal, 0 signal NLL = "<<rfr->minNll()<<" "<<rfr0->minNll()<<" n sigma="<<sqrt(2*(rfr0->minNll()-rfr->minNll()))<<endl;

  return make_pair(rfr,fullfitframe);

}


void collect_bulk_fit_results_unbinned(TString which) {
  TString outputfile;
  vector<TString> inputfiles;
  if (which=="unbinned") {
    inputfiles.push_back("DelphesEdge/3000/unbinnedfit_*.root");
    inputfiles.push_back("DelphesEdge/2000/unbinnedfit_*.root");
    inputfiles.push_back("DelphesEdge/1000/unbinnedfit_*.root");
    inputfiles.push_back("DelphesEdge/300/unbinnedfit_*.root");
    outputfile="DelphesEdge/unbinnedfits_results.root";
  }
  else if (which=="fixed0.05") {
    inputfiles.push_back("DelphesEdge/300/unbinnedFixedmass70WithConstraint0p05_*.root");
    inputfiles.push_back("DelphesEdge/1000/unbinnedFixedmass70WithConstraint0p05_*.root");
    inputfiles.push_back("DelphesEdge/2000/unbinnedFixedmass70WithConstraint0p05_*.root");
    inputfiles.push_back("DelphesEdge/3000/unbinnedFixedmass70WithConstraint0p05_*.root");
    outputfile="DelphesEdge/unbinnedFixedmass70WithConstraint0p05_results.root";
  }
  else if (which=="0.1"){
    //    inputfiles="DelphesEdge/unbinnedWithConstraint0p1_*.root";
    //outputfile="DelphesEdge/unbinnedWithConstraint0p1_results.root";
    assert(0);
  }
  else if (which=="0.05"){
    inputfiles.push_back("DelphesEdge/3000/unbinnedWithConstraint0p05_*.root");
    inputfiles.push_back("DelphesEdge/2000/unbinnedWithConstraint0p05_*.root");
    inputfiles.push_back("DelphesEdge/1000/unbinnedWithConstraint0p05_*.root");
    inputfiles.push_back("DelphesEdge/300/unbinnedWithConstraint0p05_*.root");
    outputfile="DelphesEdge/unbinnedWithConstraint0p05_results.root";
  }
  else if (which=="noedge0.05") {
    //could generalize this but for now i won't
    TString inpath="DelphesEdge/NoEdge/";
    TSystemDirectory topleveldir(inpath,inpath);
    for (int ifile=0; ifile<topleveldir.GetListOfFiles()->GetEntries(); ++ifile) {
      TSystemFile * adir = (TSystemFile*)topleveldir.GetListOfFiles()->At(ifile);
      if (adir->IsDirectory() && TString(adir->GetName()).IsDigit() ) {
	TString fullpath = inpath+"/";
	fullpath+= TString(adir->GetName());
	fullpath+="/unbinnedWithConstraint0p05*.root";
	inputfiles.push_back(fullpath);
      }
    }
    outputfile="DelphesEdge/NoEdge/unbinnedfits_results_0p05.root";
  }
  else if (which=="noedge0.05fixed70") {
    //could generalize this but for now i won't
    TString inpath="DelphesEdge/NoEdge/";
    TSystemDirectory topleveldir(inpath,inpath);
    for (int ifile=0; ifile<topleveldir.GetListOfFiles()->GetEntries(); ++ifile) {
      TSystemFile * adir = (TSystemFile*)topleveldir.GetListOfFiles()->At(ifile);
      if (adir->IsDirectory() && TString(adir->GetName()).IsDigit() ) {
	TString fullpath = inpath+"/";
	fullpath+= TString(adir->GetName());
	fullpath+="/unbinnedFixedmass70WithConstraint0p05*.root";
	inputfiles.push_back(fullpath);
	cout<<fullpath<<endl;
      }
    }
    outputfile="DelphesEdge/NoEdge/unbinnedfits_results_0p05_fixed70.root";
  }
  else assert(0);
  //load all of the unbinnedfit_xxxx_yyyyyy.root files are translate the RooFitResults into a TTree

  float edge_yield_err,edge_yield,edge_mass_err,edge_mass;
  float lumi; int status; int seed;
  TTree* resultsTree = new TTree("resultsTree","fit results");
  resultsTree->Branch("lumi",&lumi,"lumi/F");
  resultsTree->Branch("status",&status,"status/I");
  resultsTree->Branch("seed",&seed,"seed/I");
  resultsTree->Branch("edge_yield",&edge_yield,"edge_yield/F");
  resultsTree->Branch("edge_yield_err",&edge_yield_err,"edge_yield_err/F");
  resultsTree->Branch("edge_mass",&edge_mass,"edge_mass/F");
  resultsTree->Branch("edge_mass_err",&edge_mass_err,"edge_mass_err/F");
  //background parameters
  float fracvar,gtheta,n_fs,n_of,peakmass;
  resultsTree->Branch("fs_fracvar",&fracvar,"fs_fracvar/F");
  resultsTree->Branch("fs_gtheta",&gtheta,"fs_gtheta/F");
  resultsTree->Branch("fs_peakmass",&peakmass,"fs_peakmass/F");
  resultsTree->Branch("n_fs",&n_fs,"n_fs/F");
  resultsTree->Branch("n_of",&n_of,"n_of/F");
  //signif
  double minNll,minNll0;
  resultsTree->Branch("minNll",&minNll,"minNll/D");
  resultsTree->Branch("minNll0",&minNll0,"minNll0/D");
  double chiSquared;
  resultsTree->Branch("chiSquared",&chiSquared,"chiSquared/D");
  float naiveSignificance,significance2dof;
  resultsTree->Branch("naiveSignificance",&naiveSignificance,"naiveSignificance/F");
  resultsTree->Branch("significance2dof",&significance2dof,"significance2dof/F");

  //my old trick
  TChain inputdatasets("blah");
  for (vector<TString>::iterator it=inputfiles.begin(); it!=inputfiles.end(); ++it)   inputdatasets.Add(*it);
  TObjArray* inputlist=inputdatasets.GetListOfFiles();

  int weird=0;
  int ntotal= inputlist->GetSize();
  for (int ifile=0; ifile< ntotal;ifile++) {
    TObject* o = inputlist->At(ifile);
    if (o==0) continue;
    TString filename = o->GetTitle();

    TRegexp re("_[0-9]+_");//look for the _NUMBER_ pattern
    //then get rid of the _ and _
    lumi=float(TString(filename(filename.Index(re)+1,filename(re).Length()-2)).Atoi());

    seed = extractSeedFromFilename(filename);

    cout<<filename<<endl;
    TFile fin(filename);
    if (fin.IsZombie()) {
      cout<<filename<<" is bad"<<endl;
      fin.Close();
      continue;
    }
    RooFitResult * rfr=(RooFitResult*) fin.Get("fitresult");
    RooFitResult * rfr0=(RooFitResult*) fin.Get("fitresult_nosignal");
    edge_yield = ((RooRealVar*)rfr->floatParsFinal().find("n_edge"))->getVal();
    edge_yield_err = ((RooRealVar*)rfr->floatParsFinal().find("n_edge"))->getError();
    if (!which.Contains("fixed")) {
      edge_mass = ((RooRealVar*)rfr->floatParsFinal().find("m_edge"))->getVal();
      edge_mass_err = ((RooRealVar*)rfr->floatParsFinal().find("m_edge"))->getError();
    }
    fracvar = ((RooRealVar*)rfr->floatParsFinal().find("fracvar"))->getVal();
    gtheta = ((RooRealVar*)rfr->floatParsFinal().find("gtheta"))->getVal();
    peakmass = ((RooRealVar*)rfr->floatParsFinal().find("peakmass"))->getVal();
    n_fs = ((RooRealVar*)rfr->floatParsFinal().find("n_fs"))->getVal();
    n_of = ((RooRealVar*)rfr->floatParsFinal().find("n_of"))->getVal();

    minNll=rfr->minNll();
    minNll0=rfr0->minNll();
    chiSquared = 2*(minNll0-minNll);

    fin.Close();

    if (chiSquared<-0.02) { //throw these out (no, don't)
      cout<<"strange chiSquared = "<<chiSquared<<endl;
      weird++;
      naiveSignificance = -1;
      significance2dof=-1;
    }
    else {
      if (chiSquared<0) chiSquared=0; //assume these are rounding problems
      naiveSignificance = sqrt( chiSquared);
      //hard-coded solution to a numerical precision problem
      significance2dof = chiSquared<=68 ? TMath::NormQuantile(1 - 0.5* TMath::Prob(chiSquared,2)) : 8;
        
    }
    status = rfr->covQual();
    resultsTree->Fill();//keep the weird guys now too

  }

  TFile fout(outputfile,"recreate");
  resultsTree->Write();
  fout.Close();

  cout<<"weird / total = "<<weird<<" / "<<ntotal<<endl;

}

void doLeeScan(int seed) {

  const float lumi=1000;

  TString filename;
  filename.Form("DelphesEdge/NoEdge/%d/datasets_%d_%d.root",seed%1000,(int)lumi,seed);

  for (float mass = 25; mass <=95; mass+=2) {
    full_fit(lumi,filename,
	     0.05,true,mass);
  }

}

void run_Lee_scans(int min,int max) {
  for (int iseed=min; iseed<=max;iseed++) {
    doLeeScan(iseed);
  }

}

int plotLeeScan(Long_t seed) {
  //for now i am doing both of the required fits in every scan, so can load them both from each fit result file
  
  //should check that NLL0 is ~equal in every file [done by making a plot that confirms it]

  //plot -2DeltaNLL as a function of mass

  TGraph chi2;
  chi2.SetName("chi2");
  TGraph g_minNll0;
  g_minNll0.SetName("g_minNll0");

  const  TString inputpathbase = "DelphesEdge/NoEdge/";
  float lumi=1000;

  TString filenames;
  filenames.Form("%s/%d/unbinnedFixedmass*.root",inputpathbase.Data(),(int)seed%1000);
  TChain achain("dummy");
  achain.Add(filenames);
  TObjArray* inputlist=achain.GetListOfFiles();

  int ip=0;
  int n_up_crossings=0;
  const float threshold=0.5; float lastTwoDeltaNll=-1;
  for (int ifile=0; ifile< inputlist->GetSize();ifile++) {
    TObject* o = inputlist->At(ifile);
    if (o==0) continue;
    TString filename = o->GetTitle();

    TRegexp re("_[0-9]+_");//look for the _NUMBER_ pattern
    //then get rid of the _ and _
    float    l=float(TString(filename(filename.Index(re)+1,filename(re).Length()-2)).Atoi());
    if (l!=lumi) continue;
    if (extractSeedFromFilename(filename)!= seed) continue;

    TFile fin(filename);
    if (fin.IsZombie()) { fin.Close(); continue;}

    RooFitResult * rfr = (RooFitResult*) fin.Get("fitresult");
    RooFitResult* rfr0 = (RooFitResult*) fin.Get("fitresult_nosignal");

    float twoDeltaNll = 2*(rfr0->minNll() - rfr->minNll());

    TRegexp rm("Fixedmass[0-9]+");
    float themass=    TString(filename(filename.Index(rm)+9,2)).Atoi(); //this only works because all masses are 2 digits!
    g_minNll0.SetPoint(ip,themass,rfr0->minNll());
    chi2.SetPoint(ip++,themass,twoDeltaNll);


    if ( (twoDeltaNll >= threshold ) && (lastTwoDeltaNll<threshold)) n_up_crossings++;
    lastTwoDeltaNll = twoDeltaNll;
  }

  cout<<" Found n_up_crossings = "<<n_up_crossings<<endl;

  TString outfilename;
  outfilename.Form("%s/Lee%d.root",inputpathbase.Data(),(int)seed);
  TFile fout(outfilename,"recreate");
  g_minNll0.Write();
  chi2.Write();
  fout.Close();
  return n_up_crossings;
}

void bulk_plot_Lee_scan() {

  TH1D h_up_crossings("h_up_crossings","",10,0,10);
  for (int i=1;i<=32;i++) {
    h_up_crossings.Fill(    plotLeeScan(i) );
  }

  TFile f("Lee.root","recreate");
  h_up_crossings.Write();
  f.Close();

}


void plot_bulk_fit_results(TString which) { //for binned or unbinned fits

  TString inputfiles;
  TString edge_mass_dist_name;
  TString summaryfile;
  TString stub; 
  if (which=="binned") {
    inputfiles="DelphesEdge/binned_fits_bulk_*.root";
    edge_mass_dist_name = "DelphesEdge/edge_mass_distributions.pdf";
    summaryfile = "DelphesEdge/binned_fits_summary.root";
    stub="DelphesEdge/binned_";
  }
  else if (which=="unbinned") {
    inputfiles="DelphesEdge/unbinnedfits_results.root";
    edge_mass_dist_name = "DelphesEdge/unbinned_edge_mass_distributions.pdf";
    summaryfile = "DelphesEdge/unbinned_fits_summary.root";
    stub="DelphesEdge/unbinned_";
  }
  //could improve this by coding more generically e.g. which.Contains(".") and then generate the strings from the value automatically
  else if (which=="0.1") {
    inputfiles="DelphesEdge/unbinnedWithConstraint0p1_results.root";
    edge_mass_dist_name = "DelphesEdge/unbinned0p1_edge_mass_distributions.pdf";
    summaryfile = "DelphesEdge/unbinned0p1_fits_summary.root";
    stub="DelphesEdge/unbinned0p1_";
  }
  else if (which=="0.05") {
    inputfiles="DelphesEdge/unbinnedWithConstraint0p05_results.root";
    edge_mass_dist_name = "DelphesEdge/unbinned0p05_edge_mass_distributions.pdf";
    summaryfile = "DelphesEdge/unbinned0p05_fits_summary.root";
    stub="DelphesEdge/unbinned0p05_";
  }
  else if (which=="0.05fixed") {
    inputfiles="DelphesEdge/unbinnedFixedmass70WithConstraint0p05_results.root";
    edge_mass_dist_name = "";
    summaryfile = "DelphesEdge/unbinned0p05fixed70_fits_summary.root";
    stub="DelphesEdge/unbinned0p05fixed70_";
  }
  else assert(0);

  TChain r("resultsTree");
  r.Add(inputfiles);

  TGraph mass;
  TGraph mass_error;
  TGraph mass_sigma;
  TGraph mass_mean;
  TGraphErrors yield;
  TGraph yield_error;
  TGraph yield_fracerror;
  TGraph signif;
  int ip=0;
  //for each mass, draw histogram of the mass fit results and get the rms
  //as a sanity check, also plot the edge yield
 
  Long_t *lumis=0;  int nlumi = 0;
  if (which=="binned") {
    //nlumi=9;    Long_t mylumis[]={300,500,750,1000,1250,1500,2000,2500,3000};
    nlumi=7;    Long_t mylumis[]={300,500,1000,1500,2000,2500,3000};
    lumis =  new Long_t[nlumi];
    for (int ii=0;ii<nlumi;ii++) lumis[ii]=mylumis[ii];
  }
  else if (which=="unbinned") {
    nlumi=4;
    //    Long_t mylumis[]={300,500,1000,1500,2000,2500,3000};
    Long_t mylumis[]={300,1000,2000,3000};
    lumis =  new Long_t[nlumi];
    for (int ii=0;ii<nlumi;ii++) lumis[ii]=mylumis[ii];
  }
  else if (which=="0.1"){
    nlumi=4;
    Long_t mylumis[]={300,1000,2000,3000};
    lumis =  new Long_t[nlumi];
    for (int ii=0;ii<nlumi;ii++) lumis[ii]=mylumis[ii];
  }
  else if (which=="0.05" || which=="0.05fixed"){
    nlumi=4;
    Long_t mylumis[]={300,1000,2000,3000};
    lumis =  new Long_t[nlumi];
    for (int ii=0;ii<nlumi;ii++) lumis[ii]=mylumis[ii];
  }
  TFile fout(summaryfile,"recreate");

  //  int colors[]={kRed,kOrange+10,kYellow,kGreen+1,kCyan+1,kBlue,kViolet+10,kViolet,kBlack};
  int colors[]={kRed,kRed-9,kGreen+1,kGreen+4,kCyan,kBlue,kBlue+4,kMagenta,kBlack};
  for (int i=0;i<nlumi;i++) {
    Long_t lumi=lumis[i];
    TH1D hmass(TString("hmass")+lumi,"",200,0,100);
    TH1D hyield(TString("hyield")+lumi,"",10,0,100000);
    TH1D hsig(TString("hsig")+lumi,"",50,0,10);
    TH1D hChiSq(TString("hChiSq")+lumi,"2 x DeltaLogL",100,0,100);
    r.Project(TString("hmass")+lumi,"edge_mass",TString("lumi==")+lumi);
    r.Project(TString("hyield")+lumi,"edge_yield",TString("lumi==")+lumi);
    if (which == "0.05fixed") {
      r.Project(TString("hsig")+lumi,"naiveSignificance",TString("lumi==")+lumi);
      TFitResultPtr mysigfit = hsig.Fit("gaus","s");
      double signif_1dof =   mysigfit->Parameter(1);
      signif_1dof = signif_1dof*signif_1dof; //get chi2
      //now get LEE-adjusted significance
      double N = 3; //use <N>=3
      signif_1dof = TMath::NormQuantile( 1- N*exp(-0.5*signif_1dof) - 0.5*TMath::Prob(signif_1dof,1));
      signif.SetPoint(ip,lumi,signif_1dof);
    }
    else   if (which!="binned")  {
      r.Project(TString("hsig")+lumi,"significance2dof",TString("lumi==")+lumi);
      TFitResultPtr mysigfit = hsig.Fit("gaus","s");
      signif.SetPoint(ip,lumi,mysigfit->Parameter(1));

      r.Project(TString("hChiSq")+lumi,"chiSquared",TString("lumi==")+lumi);
    }

    yield.SetPoint(ip,lumi,hyield.GetMean());
    yield.SetPointError(ip,0,hyield.GetRMS());
    yield_error.SetPoint(ip,lumi,hyield.GetRMS());
    yield_fracerror.SetPoint(ip,lumi,hyield.GetRMS()/hyield.GetMean());
    mass.SetPoint(ip,lumi,hmass.GetMean());
    mass_error.SetPoint(ip,lumi,hmass.GetRMS());
    TFitResultPtr myfit = hmass.Fit("gaus","s");
    double gaus_sigma=myfit->Parameter(2);
    double gaus_mean=myfit->Parameter(1);
    mass_mean.SetPoint(ip,lumi,gaus_mean);
    mass_sigma.SetPoint(ip++,lumi,gaus_sigma);

    hmass.SetLineColor(colors[i]);
    //    if (hmass.GetMaximum()>max)    max=hmass.GetMaximum();
    //else hmass.SetMaximum(max*1.1);
    //if (i==0)  hmass.SetMaximum(hmass.GetMaximum()*2);
    //hmass.DrawCopy(drawopt); drawopt="same";
    fout.cd();    hmass.Write(); 
    if (which!="binned") {
      hsig.Write();
      hChiSq.Write();
    }
  }

  mass.SetName("mass");
  mass_error.SetName("mass_error");
  mass_sigma.SetName("mass_sigma");
  mass_mean.SetName("mass_mean");
  yield_error.SetName("yield_error");
  yield_fracerror.SetName("yield_fracerror");
  yield.SetName("yield");
  signif.SetName("significance");

  yield.GetHistogram()->SetXTitle("Integrated luminosity (fb^{-1})");
  yield_error.GetHistogram()->SetXTitle("Integrated luminosity (fb^{-1})");
  yield_fracerror.GetHistogram()->SetXTitle("Integrated luminosity (fb^{-1})");
  mass_error.GetHistogram()->SetXTitle("Integrated luminosity (fb^{-1})");
  mass_sigma.GetHistogram()->SetXTitle("Integrated luminosity (fb^{-1})");
  mass_mean.GetHistogram()->SetXTitle("Integrated luminosity (fb^{-1})");
  mass.GetHistogram()->SetXTitle("Integrated luminosity (fb^{-1})");
  signif.GetHistogram()->SetXTitle("Integrated luminosity (fb^{-1})");

  mass.GetHistogram()->SetYTitle("Edge mass (GeV)");
  mass_error.GetHistogram()->SetYTitle("Error on edge mass (GeV)");
  mass_sigma.GetHistogram()->SetYTitle("Error on edge mass (GeV)");
  mass_mean.GetHistogram()->SetYTitle("Mean edge mass (GeV)");
  yield_error.GetHistogram()->SetYTitle("Error on edge yield");
  yield.GetHistogram()->SetYTitle("Edge yield");
  yield_fracerror.GetHistogram()->SetYTitle("Fractional error on edge yield");
  if (which=="0.05fixed")
    signif.GetHistogram()->SetYTitle("LEE-corrected Significance");
  else 
    signif.GetHistogram()->SetYTitle("Significance for 2 d.o.f.");


  //add error bars to yield plot
  mass.Write();
  mass_error.Write();
  mass_sigma.Write();
  mass_mean.Write();
  yield.Write();
  yield_error.Write();
  yield_fracerror.Write();
  if (which!="binned") signif.Write();

  resetPadDimensions();//i've included my usual drawReducedTrees.h so I can use the utilities there

  renewCanvas(); 
  mass.Draw("*AL");
  thecanvas->SaveAs(stub+"fits_Mass.pdf");
  renewCanvas();
  mass_error.Draw("*AL");
  thecanvas->SaveAs(stub+"fits_MassError.pdf");
  renewCanvas();
  mass_sigma.Draw("*AL");
  thecanvas->SaveAs(stub+"fits_MassGausSigma.pdf");
  renewCanvas();
  mass_mean.Draw("*AL");
  thecanvas->SaveAs(stub+"fits_MassGausMean.pdf");
  renewCanvas();
  yield.Draw("*AL");
  thecanvas->SaveAs(stub+"fits_EdgeYield.pdf");
  renewCanvas(); 
  yield_fracerror.Draw("*AL");
  thecanvas->SaveAs(stub+"fits_EdgeYieldFracError.pdf");
  if (which!="binned") {
    renewCanvas(); 
    signif.Draw("*AL");
    thecanvas->SaveAs(stub+"fits_Significance.pdf");
  }
  fout.Close();

}

void full_fit_unbinned_bulk(float lumiToFit=-1,float constraint=-1, const bool backgroundonly=false ) {

  //my old trick
  TChain inputdatasets("blah");
  TString inpath = backgroundonly ? "DelphesEdgeNoEdge/" : "DelphesEdge/";
  inputdatasets.Add(inpath+"datasets_*.root");
  TObjArray* inputlist=inputdatasets.GetListOfFiles();

  for (int ifile=0; ifile< inputlist->GetSize();ifile++) {
    TObject* o = inputlist->At(ifile);
    if (o==0) continue;
    TString filename = o->GetTitle();

    TRegexp re("_[0-9]+_");//look for the _NUMBER_ pattern
    //then get rid of the _ and _
    float    lumi=float(TString(filename(filename.Index(re)+1,filename(re).Length()-2)).Atoi());

    if (lumiToFit>0 && (lumi!=lumiToFit)) continue;
    //    maxFits--;
    //    if (maxFits<0) break;
    
    cout<<" ======= Fitting dataset from "<<filename<<endl;
    full_fit(lumi,filename,constraint);
  }

}

void full_fit_subtracted_bulk(float lumiToFit=-1,int maxFits=9999999) {
  //hard-code for now
  bool fit_for_second_edge=false;

  TString thepath = "DelphesEdge/";
  if ( TString(gSystem->Getenv("TOYPATH"))!="") thepath = TString(gSystem->Getenv("TOYPATH"));


  bool runOnce=false;
  /*
    make a tree for the fit results

    load a dataset
    bin the data into OF and SF histograms with the appropriate Sumw2() errors
    get the SF-OF histogram
    fit the SF-OF histogram
    store the fit results to the tree (fit values, errors, lumi)
  */
  float edge_yield_err,edge_yield,edge_mass_err,edge_mass;
  float edge_yield_2,edge_mass_2;
  float lumi; float minNll;
  TTree* resultsTree = new TTree("resultsTree","fit results");
  resultsTree->Branch("lumi",&lumi,"lumi/F");
  resultsTree->Branch("minNll",&minNll,"minNll/F");
  resultsTree->Branch("edge_yield",&edge_yield,"edge_yield/F");
  resultsTree->Branch("edge_yield_err",&edge_yield_err,"edge_yield_err/F");
  resultsTree->Branch("edge_mass",&edge_mass,"edge_mass/F");
  resultsTree->Branch("edge_mass_err",&edge_mass_err,"edge_mass_err/F");
  if (fit_for_second_edge) {
    resultsTree->Branch("edge_mass_2",&edge_mass_2,"edge_mass_2/F");
    resultsTree->Branch("edge_yield_2",&edge_yield_2,"edge_yield_2/F");
  }

  //my old trick
  TChain inputdatasets("blah");
  TString inpath=thepath;
  inpath += (Long_t) TMath::Nint(lumiToFit);
  inpath += "/datasets_*.root";
  inputdatasets.Add(inpath);
  TObjArray* inputlist=inputdatasets.GetListOfFiles();

  for (int ifile=0; ifile< inputlist->GetSize();ifile++) {
    TObject* o = inputlist->At(ifile);
    if (o==0) continue;
    TString filename = o->GetTitle();

    TRegexp re("_[0-9]+_");//look for the _NUMBER_ pattern
    //then get rid of the _ and _
    lumi=float(TString(filename(filename.Index(re)+1,filename(re).Length()-2)).Atoi());

    if (lumiToFit>0 && (lumi!=lumiToFit)) continue;
    maxFits--;
    if (maxFits<0) break;
    
    cout<<" ======= Fitting dataset from "<<filename<<endl;

    TFile frds(filename);
    RooDataSet* data_of = (RooDataSet*) frds.Get("data_of");
    RooDataSet* data_sf_ee = (RooDataSet*) frds.Get("data_sf_ee");
    RooDataSet* data_sf_mm = (RooDataSet*) frds.Get("data_sf_mm");

    float hmin=20;
    float hmax=100;
    int nbins = (hmax-hmin)/2; //2 gev bins
    TH1D h_of("h_of","of",nbins,hmin,hmax);
    TH1D h_sf("h_sf","sf",nbins,hmin,hmax);
    TH1D h_subtracted("h_subtracted","sf-of",nbins,hmin,hmax);
    h_of.Sumw2();
    h_sf.Sumw2();
    h_subtracted.Sumw2();
    //i can't find a better way to do this
    for (int ie=0;ie<data_of->numEntries();ie++) {
      double val = ((RooRealVar*)data_of->get(ie)->find("mll"))->getVal();
      if (val>=hmin && val<hmax)  h_of.Fill(val );
    }

    for (int ie=0;ie<data_sf_ee->numEntries();ie++) {
      double val = ((RooRealVar*)data_sf_ee->get(ie)->find("mll"))->getVal();
      if (val>=hmin && val<hmax)  h_sf.Fill(val );
    }
    for (int ie=0;ie<data_sf_mm->numEntries();ie++) {
      double val = ((RooRealVar*)data_sf_mm->get(ie)->find("mll"))->getVal();
      if (val>=hmin && val<hmax)  h_sf.Fill(val );
    }

    h_subtracted.Add(&h_sf,&h_of,1,-1);

    std::pair<RooFitResult*,RooPlot*> fitpair = do_binned_fit(h_subtracted,lumi,fit_for_second_edge);
    RooFitResult* fitresult = fitpair.first;
    RooPlot* fitplot = fitpair.second;
    edge_yield = ((RooRealVar*)fitresult->floatParsFinal().find("n_edge"))->getVal();
    edge_yield_err = ((RooRealVar*)fitresult->floatParsFinal().find("n_edge"))->getError();
    edge_mass = ((RooRealVar*)fitresult->floatParsFinal().find("m_edge"))->getVal();
    edge_mass_err = ((RooRealVar*)fitresult->floatParsFinal().find("m_edge"))->getError();
    if (fit_for_second_edge) {
      edge_yield_2 = ((RooRealVar*)fitresult->floatParsFinal().find("n_edge2"))->getVal();
      edge_mass_2 = ((RooRealVar*)fitresult->floatParsFinal().find("m_edge2"))->getVal();
    }
    minNll =fitresult->minNll() ;
    resultsTree->Fill();

    TString plotname = filename;
    TString replacestring = fit_for_second_edge ? "binnedfitplot2_" : "binnedfitplot_";
    plotname.ReplaceAll("datasets_",replacestring);
    plotname = plotname(plotname.Last('/')+1,plotname.Length()-6-plotname.Last('/')); //assumes ends in .root
    cout<<"plotname = "<<plotname<<endl;
    fitplot->SetName(plotname);
    TString outfilename;
    if (lumiToFit>0)   outfilename.Form("%sbinnedfitplots%s_%d.root",thepath.Data(),fit_for_second_edge ? "2" : "",TMath::Nint(lumiToFit));
    else               outfilename= fit_for_second_edge ?thepath+"binnedfitplots2.root" : thepath+"binnedfitplots.root";
    TString filemode = runOnce ? "update":"recreate";
    TFile fp(outfilename,filemode);
    fitplot->Write();
    fp.Close();
    frds.Close();
    runOnce=true;//switch to using 'update' for subsequent iterations
  }

  TString treefilename;
  if (lumiToFit>0)   treefilename.Form("%sbinned_fits_bulk%s_%d.root",thepath.Data(),fit_for_second_edge ? "2" : "",TMath::Nint(lumiToFit));
  else               treefilename=fit_for_second_edge ? thepath+"binned_fits_bulk2.root" : thepath+"binned_fits_bulk.root";
  TFile fout(treefilename,"recreate");
  resultsTree->Write();
  fout.Close();

}

void full_fit_subtracted_test() {

  /*
fit the SF-OF data
this is necessarily a binned fit
  */

  TFile fin("DelphesEdge/sfMinusOf-fine.root");
  TH1D* total_sf = (TH1D*) fin.Get("total_sf");
  TH1D* total_of = (TH1D*) fin.Get("total_of");
  //these histos are event counts
  //reset errors on these histos to be sqrt(N)
  //  TH1D* subtracted = new TH1D("subtracted","",total_sf->GetNbinsX(),total_sf->GetXaxis()->GetXmin(),total_sf->GetXaxis()->GetXmax());
  {for (int i=1; i<=total_sf->GetNbinsX(); i++) {
    total_sf->SetBinError(i,sqrt(total_sf->GetBinContent(i)));
    total_of->SetBinError(i,sqrt(total_of->GetBinContent(i)));
    }}
  //now subtract them to get SF-OF histogram
  //  gROOT->cd();
  TH1D* subtracted= (TH1D*)  total_sf->Clone("subtracted");
  subtracted->Reset(); subtracted->Sumw2();
  subtracted->Add(total_sf,total_of,1,-1);

  //construct PDF of N_sig*P(Edge) + N_Z*P(Z)
  //binned fit from 20<mll<100
  RooRealVar mll("mll","mll",20,100);
  RooDataHist data("data","sub data",RooArgList(mll),subtracted);

  //Z PDF
  RooRealVar cexp("cexp","exponential const",0.01,-10,10); //for the RooExponential
  RooExponential pdf_exp("pdf_exp","exp pdf",mll,cexp);
  //now the Voigtian
  RooRealVar zmass("zmass","z mass",91,80,100);
  RooRealVar zwidth("zwidth","z width", 2.5,0.1,10);
  RooRealVar zsigma("zsigma","z resolution",1e-2,0.00001,20);
  RooVoigtian pdf_peak("pdf_peak","z peak pdf",mll,zmass,zwidth,zsigma,kTRUE);

  RooRealVar fexp("fexp","exp fraction",1e-2);
  RooAddPdf pdf_dy("pdf_dy","total dy pdf",pdf_exp,pdf_peak,fexp);

  TFile fdyfit("EdgeFit_DY.root");
  RooFitResult * Zfit = (RooFitResult*) fdyfit.Get("Zfit");
  RooArgList dyfitargs=Zfit->floatParsFinal();
  cexp.setVal(((RooRealVar*)dyfitargs.find("cexp"))->getVal());
  fexp.setVal(((RooRealVar*)dyfitargs.find("fexp"))->getVal());
  zmass.setVal(((RooRealVar*)dyfitargs.find("zmass"))->getVal());
  zsigma.setVal(((RooRealVar*)dyfitargs.find("zsigma"))->getVal());
  zwidth.setVal(((RooRealVar*)dyfitargs.find("zwidth"))->getVal());
  //these fit parameters are then fixed
  cexp.setConstant();
  zmass.setConstant();
  zwidth.setConstant();
  zsigma.setConstant();
  fexp.setConstant();

  //signal
  RooRealVar m_edge("m_edge","edge mass",68,20,100);
  RooRealVar edge_sigma("edge_sigma","edge resolution",2,0.01,9); //update to new value
  RooEdge pdf_edge("pdf_edge","edge pdf",mll,m_edge,edge_sigma);
  //m_edge.setConstant();
  edge_sigma.setConstant();

  //yields
  RooRealVar n_edge("n_edge","edge yield",100,-5,1000);
  RooRealVar n_dy("n_dy","Z yield",50,-5,1000);

  RooArgList pdfs(pdf_dy,pdf_edge);
  RooArgList yields(n_dy,n_edge);
  RooAddPdf pdf_total("pdf_total","total pdf",pdfs,yields);
  RooFitResult* rfr=  pdf_total.fitTo(data,RooFit::Save(1),RooFit::SumW2Error(kTRUE));
  rfr->Print();

  RooPlot* fullfitframe = mll.frame();
  data.plotOn(fullfitframe);
  pdf_total.plotOn(fullfitframe);
  fullfitframe->SetName("fullfitframe");


  fin.Close();
}

void test(){

  double edge = 70;
  double sigma = 2;

  TGraph g;
  int i=1;
  //loop over x (mass to draw)
  {
  for (double x = 10; x<=100; x+=1) {
    TF1 F("F","exp(- ([0]-[1]*sqrt(0.5*(1+x)))*([0]-[1]*sqrt(0.5*(1+x))) * (1/(2*([2]*[2])) ))",-1,1);
    
    F.SetParameter(0,x);
    F.SetParameter(1,edge);
    F.SetParameter(2,sigma);

    double P = F.Integral(-1,1);
    g.SetPoint(i++,x,P);
  }
  }

}

void testOF() {

  double m1 = 50;
  double m2 = 120;

  double c20 = 200;
  double c21 = 120;
  double c22 = -1.3;
  double c23 = 0.01;

  double alpha = (m1*c21 + 2*c22*m1*m1 + 3*c23*m1*m1*m1)/ (c20 + c21*m1+c22*m1*m1+c23*m1*m1*m1);
  double beta = -(c21+2*c22*m2+3*c23*m2*m2)/(c20+c21*m2+c22*m2*m2+c23*m2*m2*m2); 

  double c1 = (c20 + c21*m1 + c22*m1*m1 + c23*m1*m1*m1) / pow(m1,alpha);
  double c3 = (c20 + c21*m2 + c22*m2*m2 + c23*m2*m2*m2) *exp(beta*m2);

  TGraph g; int i=1;
  {
  for (double m= 20; m<= 300; m+= 1) {

    double P=0;
    if (m<m1) {
      P = c1*pow(m,alpha);
    }
    else if (m>=m1 && m<m2) {
      P = c20 + c21*m + c22*m*m + c23*m*m*m;
    }
    else {
      P = c3 *exp(-beta*m);
    }
    g.SetPoint(i++,m,P);
  }
  }


}

void testOF2() {

  double alpha = 1;
  double beta = 0.01;

  TGraph g; int i=1;
  {
  for (double m= 20; m<= 300; m+= 1) {

    double P= pow(m,alpha)* exp(-beta*m);
    g.SetPoint(i++,m,P);

  }}


}


