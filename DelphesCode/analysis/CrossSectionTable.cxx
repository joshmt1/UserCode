#include "CrossSectionTable.h"

/*
Class that loads a text file containing NLO cross sections across the mSugra plane
format is:

m0 m12 s0 s1 s2 s3 s4 s5 s6 s7 s8 s9

where the sN are the cross-sections for the production sub-processes given in the 
enum defined in the header file

-- new addition: handle SMS cross sections in a format such as seen here:
https://twiki.cern.ch/twiki/pub/LHCPhysics/SUSYCrossSections/stst_decoupled7TeV.txt

also handle SMS cross sections in a ROOT input file

-- another addition: handle pMSSM cross sections in the format seen here:
http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/pMSSMInterpretation2013/parameters/

*/

#include "TObjArray.h"
#include "TFile.h"
#include "TH1.h"

#include <string>

#include "TString.h"

CrossSectionTable::CrossSectionTable(const TString & inputFile, const TString & format,const TString & histoname) //:
//  filename_(inputFile) {
{
  //load the contents of file into the database
  if (format=="smstext") loadFileToDatabaseSMS(inputFile);
  else if (format=="smsroot") loadFileToDatabaseSMSRoot(inputFile,histoname);
  else if (format=="CMSSM") loadFileToDatabase(inputFile);
  else if (format=="pMSSM") loadFileToDatabasePMSSM(inputFile,false);
  else if (format=="simplesms") loadFileToDatabaseSMS_simple(inputFile);
  else assert(0);

}

CrossSectionTable::~CrossSectionTable() {}

std::map<SUSYProcess,double>  CrossSectionTable::operator[] (const std::pair<int,int> & scanpoint) const {

  std::map<std::pair<int, int>, std::map<SUSYProcess, double> >::const_iterator mapval = database_.find(scanpoint);

  if (mapval == database_.end() ) {
    std::cout<<"[CrossSectionTable] did not find cross section for scan point "<<scanpoint.first<<" "<<scanpoint.second<<std::endl;
  }
  return mapval->second;

}

double CrossSectionTable::getCrossSection(const int m0, const int m12, const SUSYProcess process) const {

  std::map<SUSYProcess,double> forThisPoint = (*this)[std::make_pair(m0,m12)];
  return forThisPoint[process];

}

double CrossSectionTable::getSMSCrossSection(const int m0) const {

  std::map<SUSYProcess,double> forThisPoint = (*this)[std::make_pair(m0,0)];
  return forThisPoint[NotFound];

}

void CrossSectionTable::loadFileToDatabasePMSSM(const TString & filename,const bool append) {
  using namespace std;

  cout<<"Loading file "<<filename<<" into CrossSectionTable"<<endl;
  int idnumber;
  string pointname;
  double xsec; TString xsecin;
  string dummy1,dummy2;

  ifstream file10(filename.Data());
  if (!file10.good() ) {
    cout<<"Problem with file! Terminating..."<<endl;
    assert(0);
  }

  if (!append)  database_.clear();

  //load the header row
  file10>>dummy1>>pointname>>dummy2;
  int nloaded=0;
  while (file10>>idnumber>>pointname>>xsecin ) {

    if (xsecin.IsFloat() ) xsec=xsecin.Atof();
    else {
      cout<<" invalid cross-section "<<xsecin<<" for point "<<idnumber<<" "<<pointname<<endl;
      xsec=0;
    }

    //this class was built for mass pairs; here we just have an index
    pair<int,int> thispoint = make_pair(idnumber,0);
    map<SUSYProcess, double> theseCrossSections;
    //use NotFound 
    theseCrossSections[NotFound] = xsec;
    database_[thispoint] = theseCrossSections;
    nloaded++;
  }
  file10.close();
  cout<<"Loaded cross sections for "<<nloaded<<" scan points. Table size = "<<database_.size()<<endl;

}

void CrossSectionTable::loadFileToDatabase(const TString & filename) {
  using namespace std;

  cout<<"Loading file "<<filename<<" into CrossSectionTable"<<endl;


  int m0,m12;
  double c1,c2,c3,c4,c5,c6,c7,c8,c9,c10;
  ifstream file10(filename.Data());
  if (!file10.good() ) {
    cout<<"Problem with file! Terminating..."<<endl;
    assert(0);
  }
  database_.clear(); //just to be sure
  while ( file10>>m0>>m12>>c1>>c2>>c3>>c4>>c5>>c6>>c7>>c8>>c9>>c10 ) {
    pair<int, int> scanpoint = make_pair(m0,m12);

    map<SUSYProcess, double> theseCrossSections;
    theseCrossSections[ng] = c1;
    theseCrossSections[ns] = c2;
    theseCrossSections[nn] = c3;
    theseCrossSections[ll] = c4;
    theseCrossSections[sb] = c5;
    theseCrossSections[ss] = c6;
    theseCrossSections[tb] = c7;
    theseCrossSections[bb] = c8;
    theseCrossSections[gg] = c9;
    theseCrossSections[sg] = c10;
    database_[scanpoint] = theseCrossSections;
  }
  file10.close();
  cout<<"Loaded cross sections for "<<database_.size()<<" scan points"<<endl;

}

void CrossSectionTable::loadFileToDatabaseSMSRoot(const TString & filename,const TString & histoname) { //SMS from root file
  using namespace std;
  const int m12 = 0; //no change as a function of m12 axis; used fixed value. 0 seems logical

  TFile fxs(filename);
  TH1D* h = (TH1D*) fxs.Get(histoname);
  if (h!=0) {

    for (int ibin=1; ibin<=h->GetNbinsX(); ibin++) {
      int m0 = h->GetBinCenter(ibin); //careful that this choice matches what is used in the source histogram
      double xs = h->GetBinContent(ibin);

      pair<int, int> scanpoint = make_pair(m0,m12);
      map<SUSYProcess, double> theseCrossSections;
      //i could do this not as NotFound but as the actual SUSY process (tb for stop-stop, gg for gluino-gluino)
      //but i think that is not needed
      //this choice is connected to the hard-coded use of NotFound in getSMSCrossSection()
      theseCrossSections[NotFound] = xs;
      database_[scanpoint] = theseCrossSections;
    }
  }
  else cout << "Problem with file or histogram!"<<endl; 

  fxs.Close();

}

void CrossSectionTable::loadFileToDatabaseSMS_simple(const TString & filename) {

  using namespace std;
  const int m12 = 0; //no change as a function of m12 axis; used fixed value. 0 seems logical

  ifstream file10(filename.Data());
  //simple 2 column format
  // mass cross-section

  if (!file10.good() ) {
    cout<<"Problem with file! Terminating..."<<endl;
    assert(0);
  }

  database_.clear(); //just to be sure
  int mass;
  double xs;
  while ( file10>>mass>>xs ) {
    pair<int, int> scanpoint = make_pair(mass,m12);
    map<SUSYProcess, double> theseCrossSections;
    //see comments elsewhere about choice of NotFound
    theseCrossSections[NotFound] = xs;
    database_[scanpoint] = theseCrossSections;
  }

  file10.close();
  cout<<"Loaded cross sections for "<<database_.size()<<" scan points"<<endl;

}

void CrossSectionTable::loadFileToDatabaseSMS(const TString & filename) {
  using namespace std;

  cout<<"Loading SMS file "<<filename<<" into CrossSectionTable"<<endl;

  const int m12 = 0; //no change as a function of m12 axis; used fixed value. 0 seems logical
  //nb -- connect to value hard-coded in getSMSCrossSection()

  //these files have a rather annoying format
  ifstream file10(filename.Data());
  string line;
  int linenum=0;
  if (file10.is_open()) {
    while ( file10.good() )   {
      getline (file10,line);
      //      cout<<" got a line: "<<line<<endl;
      if (linenum>0 && line.length() > 1) { //skip the first line and blank lines
	TString thisline = line.c_str(); //TStrings have Tokenize, etc
	TString mass = thisline.Tokenize("|")->At(1)->GetName(); // " 100 GeV "
	TString xs_pluscrap = thisline.Tokenize("|")->At(2)->GetName();
	//now get the numbers alone
	TString massval=mass.Tokenize(" ")->At(0)->GetName();
	TString xsval=xs_pluscrap.Tokenize(" ")->At(0)->GetName();
	//convert to numeric format
	int   m0  = massval.Atoi();
	double xs = xsval.Atof();

	pair<int, int> scanpoint = make_pair(m0,m12);
	map<SUSYProcess, double> theseCrossSections;
	//i could do this not as NotFound but as the actual SUSY process (tb for stop-stop, gg for gluino-gluino)
	//but i think that is not needed
	//this choice is connected to the hard-coded use of NotFound in getSMSCrossSection()
	theseCrossSections[NotFound] = xs;
	database_[scanpoint] = theseCrossSections;
      }
      ++linenum;
    }
    file10.close();
  }
  else cout << "Unable to open file"<<endl; 

  file10.close();
  cout<<"Loaded cross sections for "<<database_.size()<<" scan points"<<endl;

}
