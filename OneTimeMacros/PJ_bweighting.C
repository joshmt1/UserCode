#include "btagToy.C"

double jetTagEff_ewk(const pair<jetFlavor, tagStatus> & jet,const tagStatus workingPoint) {

  //this should be eff * SF for a jet

  //i think this is right. neither of these depends on how the jet was actually tagged
  return getEffMC( jet.first, workingPoint) * getSF(jet.first, workingPoint);

}


void calculateTagProb_ewk(vector<pair<jetFlavor,tagStatus> > & event, float &Prob2b, float &Prob3b, float &Prob4b){ 
  //this should take in an 'event' (vector of jets) and return the prob* variables

  //  float pTthresh = 20;
  //Init
  Prob2b = 0;
  Prob3b = 0;
  Prob4b = 0;
  
  /* comment out not needed parts 
  if (f_tageff_ == 0) {
    Prob2b = -1;
    Prob3b = -1;
    Prob4b = -1;
  } 
  */

  /*
  char btageffname[200], ctageffname[200], ltageffname[200];
  std::string sbtageff = "h_btageff";  std::string sctageff = "h_ctageff";  std::string sltageff = "h_ltageff";
  sprintf(btageffname,"%s",sbtageff.c_str());   
  sprintf(ctageffname,"%s",sctageff.c_str());   
  sprintf(ltageffname,"%s",sltageff.c_str());   
  TH1D * h_btageff  = (TH1D *)f_tageff_->Get(btageffname);
  TH1D * h_ctageff  = (TH1D *)f_tageff_->Get(ctageffname);
  TH1D * h_ltageff  = (TH1D *)f_tageff_->Get(ltageffname);
  */

  unsigned int jetsize = event.size() ;

  if ( jetsize < 1 ) return;

  //cout << "jetsize " << jetsize << endl;
  //2b
  for (unsigned int ijet=0; ijet<jetsize-1; ++ijet) {
    //    if(!isGoodJet(ijet,pTthresh)) continue; //switch to 20 GeV threshold for b jets
    double effi = jetTagEff_ewk( event.at(ijet), kT);
    for (unsigned int jjet=ijet+1; jjet<jetsize; ++jjet) {
      if (jjet == ijet ) continue; 
      //      if(!isGoodJet(jjet,pTthresh)) continue;
      double effj = jetTagEff_ewk( event.at(jjet), kT);
      double prod2b = 1;
      for (unsigned int kjet=0; kjet<jetsize; ++kjet) {
        if( (kjet == jjet) || (kjet == ijet) ) continue;
	//        if(!isGoodJet(kjet,pTthresh)) continue;
        double effk = jetTagEff_ewk( event.at(kjet), kM);
        double effkt = jetTagEff_ewk( event.at(kjet), kL);
	prod2b *= (1-effk);
	//PJ already commented this
	//double prod3b = 1; 
        //for (unsigned int ljet=0; ljet<jetsize; ++ljet) {
        //  if( (ljet == kjet) || (ljet == ijet) || (ljet == jjet) ) continue;
        //  if(isGoodJet(ljet,pTthresh)) continue;
        //  double effl = jetTagEff_ewk(ljet, kL); // Loose
	  //prod3b *= (1-effl);
        //} //l loop
        //Prob3b += prod3b*effi*effj*(effk-effkt);
      }// k loop
      Prob2b += prod2b*effi*effj;
    }//j loop
  }//i loop


  if ( jetsize >= 4 ) {
 
  //3b_part1
  for (unsigned int ijet=0; ijet<jetsize-1; ++ijet) {
    //    if(!isGoodJet(ijet,pTthresh)) continue; //switch to 20 GeV threshold for b jets
    double effi = jetTagEff_ewk( event.at(ijet), kT);
    for (unsigned int jjet=ijet+1; jjet<jetsize; ++jjet) {
      if (jjet == ijet ) continue; 
      //      if(!isGoodJet(jjet,pTthresh)) continue;
      double effj = jetTagEff_ewk(event.at(jjet), kT);
      for (unsigned int kjet=0; kjet<jetsize; ++kjet) {
        if( (kjet == jjet) || (kjet == ijet) ) continue;
	//        if(!isGoodJet(kjet,pTthresh)) continue;
        double effk = jetTagEff_ewk( event.at(kjet), kM);
        double effkt = jetTagEff_ewk( event.at(kjet), kL);
	double prod3b = 1; 
        for (unsigned int ljet=0; ljet<jetsize; ++ljet) {
          if( (ljet == kjet) || (ljet == ijet) || (ljet == jjet) ) continue;
	  //          if(isGoodJet(ljet,pTthresh)) continue;
          double effl = jetTagEff_ewk( event.at(ljet), kL);
	  prod3b *= (1-effl);
        } //l loop
        Prob3b += prod3b*effi*effj*(effk-effkt);
      }// k loop
    }//j loop
  }//i loop

  //3b_part2
  for (unsigned int ijet=0; ijet<jetsize-2; ++ijet) {
    //    if(!isGoodJet(ijet,pTthresh)) continue; //switch to 20 GeV threshold for b jets
    double effi = jetTagEff_ewk( event.at(ijet), kT);
    for (unsigned int jjet=ijet+1; jjet<jetsize-1; ++jjet) {
      if (jjet == ijet ) continue; 
      //      if(!isGoodJet(jjet,pTthresh)) continue;
      double effj = jetTagEff_ewk( event.at(jjet), kT);
      for (unsigned int kjet=jjet+1; kjet<jetsize; ++kjet) {
        if( (kjet == jjet) || (kjet == ijet) ) continue;
	//        if(!isGoodJet(kjet,pTthresh)) continue;
        double effk = jetTagEff_ewk( event.at(kjet), kT);
	double prod3b = 1; 
        for (unsigned int ljet=0; ljet<jetsize; ++ljet) {
          if( (ljet == kjet) || (ljet == ijet) || (ljet == jjet) ) continue;
	  //          if(isGoodJet(ljet,pTthresh)) continue;
          double effl = jetTagEff_ewk( event.at( ljet), kL);
	  prod3b *= (1-effl);
        } //l loop
        Prob3b += prod3b*effi*effj*effk;
      }// k loop
    }//j loop
  }//i loop
  } //ensure 4jets

  if ( jetsize == 4 ) {

  //4T
  for (unsigned int ijet=0; ijet<jetsize-3; ++ijet) {
    //    if(!isGoodJet(ijet,pTthresh)) continue; //switch to 20 GeV threshold for b jets
    double effi = jetTagEff_ewk( event.at(ijet), kT);
    for (unsigned int jjet=ijet+1; jjet<jetsize-2; ++jjet) {
      //      if(!isGoodJet(jjet,pTthresh)) continue;
      double effj = jetTagEff_ewk(event.at(jjet), kT);
      for (unsigned int kjet=jjet+1; kjet<jetsize-1; ++kjet) {
	//        if(!isGoodJet(kjet,pTthresh)) continue;
        double effk = jetTagEff_ewk( event.at(kjet), kT);
        for (unsigned int ljet=kjet+1; ljet<jetsize; ++ljet) {
	  //          if(isGoodJet(ljet,pTthresh)) continue;
          double effl = jetTagEff_ewk(event.at(ljet), kT);
          Prob4b += effi*effj*effk*effl;
        } //l loop
      }// k loop
    }//j loop
  }//i loop

  //3T, 1(L-T)
  for (unsigned int ijet=0; ijet<jetsize-2; ++ijet) {
    //    if(!isGoodJet(ijet,pTthresh)) continue; //switch to 20 GeV threshold for b jets
    double effi = jetTagEff_ewk(event.at(ijet), kT);
    for (unsigned int jjet=ijet+1; jjet<jetsize-1; ++jjet) {
      //      if(!isGoodJet(jjet,pTthresh)) continue;
      double effj = jetTagEff_ewk(event.at(jjet), kT);
      for (unsigned int kjet=jjet+1; kjet<jetsize; ++kjet) {
	//        if(!isGoodJet(kjet,pTthresh)) continue;
        double effk = jetTagEff_ewk(event.at(kjet), kT);
        double prod4b = 1;
        for (unsigned int ljet=0; ljet<jetsize; ++ljet) {
	  //          if(isGoodJet(ljet,pTthresh)) continue;
          if( (ljet == kjet) || (ljet == ijet) || (ljet == jjet) ) continue;
          double effl  = jetTagEff_ewk(event.at(ljet), kL);
          double efflt = jetTagEff_ewk(event.at(ljet), kT);
	  prod4b *= (effl-efflt);
        } //l loop
        Prob4b += prod4b*effi*effj*effk;
      }// k loop
    }//j loop
  }//i loop

  //2T, 2(M-T)
  for (unsigned int ijet=0; ijet<jetsize-1; ++ijet) {
    //    if(!isGoodJet(ijet,pTthresh)) continue; //switch to 20 GeV threshold for b jets
    double effi = jetTagEff_ewk(event.at(ijet),kT);
    for (unsigned int jjet=ijet+1; jjet<jetsize; ++jjet) {
      //      if(!isGoodJet(jjet,pTthresh)) continue;
      double effj = jetTagEff_ewk(event.at(jjet), kT);
      double prod1 = 1;
      for (unsigned int kjet=0; kjet<jetsize-1; ++kjet) {
	//        if(!isGoodJet(kjet,pTthresh)) continue;
        if( (kjet == ijet) || (kjet == jjet) ) continue;
        double effkt = jetTagEff_ewk(event.at(kjet),kT);
        double effkm = jetTagEff_ewk(event.at(kjet),kM);
	prod1 *= (effkm-effkt);
        double prod4b = 1;
        for (unsigned int ljet=kjet+1; ljet<jetsize; ++ljet) {
	  //          if(isGoodJet(ljet,pTthresh)) continue;
          if( (ljet == ijet) || (ljet == jjet) ) continue;
          double efflt = jetTagEff_ewk(event.at(ljet), kT); // Tight
          double efflm = jetTagEff_ewk(event.at(ljet), kM); // Medium
	  prod4b *= (efflm-efflt);
        } //l loop
        Prob4b += prod1*prod4b*effi*effj;
      }// k loop
    }//j loop
  }//i loop

  //2T, 1(M-T), 1(L-M)
  for (unsigned int ijet=0; ijet<jetsize-1; ++ijet) {
    //    if(!isGoodJet(ijet,pTthresh)) continue; //switch to 20 GeV threshold for b jets
    double effi = jetTagEff_ewk(event.at(ijet), kT); // Tight
    for (unsigned int jjet=ijet+1; jjet<jetsize; ++jjet) {
      //      if(!isGoodJet(jjet,pTthresh)) continue;
      double effj = jetTagEff_ewk(event.at(jjet), kT); // Tight
      double prod1 = 1;
      for (unsigned int kjet=0; kjet<jetsize; ++kjet) {
	//        if(!isGoodJet(kjet,pTthresh)) continue;
        if( (kjet == ijet) || (kjet == jjet) ) continue;
        double effkt = jetTagEff_ewk(event.at(kjet), kT); // Tight
        double effkm = jetTagEff_ewk(event.at(kjet), kM); // Medium
	prod1 *= (effkm-effkt);
        double prod4b = 1;
        for (unsigned int ljet=0; ljet<jetsize; ++ljet) {
	  //          if(isGoodJet(ljet,pTthresh)) continue;
          if( (ljet == ijet) || (ljet == jjet) || (ljet == kjet) ) continue;
          double efflt = jetTagEff_ewk(event.at(ljet), kT); // Tight
          double efflm = jetTagEff_ewk(event.at(ljet), kM); // Medium
	  prod4b *= (efflt-efflm);
        } //l loop
        Prob4b += prod4b*effi*effj;
      }// k loop
    }//j loop
  }//i loop

  } //4jets only

  if ( jetsize >= 5 ) {

  //5T
  for (unsigned int ijet=0; ijet<jetsize-4; ++ijet) {
    //    if(!isGoodJet(ijet,pTthresh)) continue; //switch to 20 GeV threshold for b jets
    double effi = jetTagEff_ewk(event.at(ijet), kT); // Tight
    for (unsigned int jjet=ijet+1; jjet<jetsize-3; ++jjet) {
      //      if(!isGoodJet(jjet,pTthresh)) continue;
      double effj = jetTagEff_ewk(event.at(jjet), kT); // Tight
      for (unsigned int kjet=jjet+1; kjet<jetsize-2; ++kjet) {
	//        if(!isGoodJet(kjet,pTthresh)) continue;
        double effk = jetTagEff_ewk(event.at(kjet), kT); // Tight
        for (unsigned int ljet=kjet+1; ljet<jetsize-1; ++ljet) {
	  //          if(isGoodJet(ljet,pTthresh)) continue;
          double effl = jetTagEff_ewk(event.at(ljet), kT); // Tight
          for (unsigned int mjet=ljet+1; mjet<jetsize; ++mjet) {
	    //            if(isGoodJet(mjet,pTthresh)) continue;
            double effm = jetTagEff_ewk(event.at(mjet), kT); // Tight
            Prob4b += effi*effj*effk*effl*effm;
          } //m loop
        } //l loop
      }// k loop
    }//j loop
  }//i loop

  //4T (1-T)
  for (unsigned int ijet=0; ijet<jetsize-3; ++ijet) {
    //    if(!isGoodJet(ijet,pTthresh)) continue; //switch to 20 GeV threshold for b jets
    double effi = jetTagEff_ewk(event.at(ijet), kT); // Tight
    for (unsigned int jjet=ijet+1; jjet<jetsize-2; ++jjet) {
      //      if(!isGoodJet(jjet,pTthresh)) continue;
      double effj = jetTagEff_ewk(event.at(jjet), kT); // Tight
      for (unsigned int kjet=jjet+1; kjet<jetsize-1; ++kjet) {
	//        if(!isGoodJet(kjet,pTthresh)) continue;
        double effk = jetTagEff_ewk(event.at(kjet), kT); // Tight
        for (unsigned int ljet=kjet+1; ljet<jetsize; ++ljet) {
	  //          if(isGoodJet(ljet,pTthresh)) continue;
          double effl = jetTagEff_ewk(event.at(ljet), kT); // Tight
	  double prod1 = 1;
          for (unsigned int mjet=0; mjet<jetsize; ++mjet) {
	    //            if(isGoodJet(mjet,pTthresh)) continue;
	    if ( (mjet==ijet)||(mjet==jjet)||(mjet==kjet)||(mjet==ljet) ) continue; 
            double effm = jetTagEff_ewk(event.at(mjet), kT); // Tight
	    prod1 *= (1-effm);
          } //m loop
          Prob4b += effi*effj*effk*effl*prod1;
        } //l loop
      }// k loop
    }//j loop
  }//i loop

  //3T (M-T) (M-T)
  //3T (L-M) (L-M)
  for (unsigned int ijet=0; ijet<jetsize-2; ++ijet) {
    //    if(!isGoodJet(ijet,pTthresh)) continue; //switch to 20 GeV threshold for b jets
    double effi = jetTagEff_ewk(event.at(ijet), kT); // Tight
    for (unsigned int jjet=ijet+1; jjet<jetsize-1; ++jjet) {
      //      if(!isGoodJet(jjet,pTthresh)) continue;
      double effj = jetTagEff_ewk(event.at(jjet), kT); // Tight
      for (unsigned int kjet=jjet+1; kjet<jetsize; ++kjet) {
	//        if(!isGoodJet(kjet,pTthresh)) continue;
        double effk = jetTagEff_ewk(event.at(kjet), kT); // Tight
        double prod1 = 1;
        double prod3 = 1;
        for (unsigned int ljet=0; ljet<jetsize-1; ++ljet) {
	  //          if(isGoodJet(ljet,pTthresh)) continue;
	  if ( (ljet==ijet)||(ljet==jjet)||(ljet==kjet) ) continue; 
          double efflt = jetTagEff_ewk(event.at(ljet), kT); // Tight
          double efflm = jetTagEff_ewk(event.at(ljet), kM); // Medium 
          double effll = jetTagEff_ewk(event.at(ljet), kL); // Loose 
	  prod1 *= (efflm-efflt);
	  prod3 *= (effll-efflm);
	  double prod2 = 1;
	  double prod4 = 1;
          for (unsigned int mjet=ljet+1; mjet<jetsize; ++mjet) {
	    //            if(isGoodJet(mjet,pTthresh)) continue;
	    if ( (mjet==ijet)||(mjet==jjet)||(mjet==kjet) ) continue; 
            double effmt = jetTagEff_ewk(event.at(mjet), kT); // Tight
            double effmm = jetTagEff_ewk(event.at(mjet), kM); // Medium 
            double effml = jetTagEff_ewk(event.at(mjet), kL); // Loose 
	    prod2 *= (effmm-effmt);
	    prod4 *= (effml-effmm);
          } //m loop
          Prob4b += effi*effj*effk*prod1*prod2;
          Prob4b += effi*effj*effk*prod3*prod4;
        } //l loop
      }// k loop
    }//j loop
  }//i loop

  //3T (M-T) (1-M)
  //3T (L-M) (1-L)
  for (unsigned int ijet=0; ijet<jetsize-2; ++ijet) {
    //    if(!isGoodJet(ijet,pTthresh)) continue; //switch to 20 GeV threshold for b jets
    double effi = jetTagEff_ewk(event.at(ijet), kT); // Tight
    for (unsigned int jjet=ijet+1; jjet<jetsize-1; ++jjet) {
      //      if(!isGoodJet(jjet,pTthresh)) continue;
      double effj = jetTagEff_ewk(event.at(jjet), kT); // Tight
      for (unsigned int kjet=jjet+1; kjet<jetsize; ++kjet) {
	//        if(!isGoodJet(kjet,pTthresh)) continue;
        double effk = jetTagEff_ewk(event.at(kjet), kT); // Tight
        double prod1 = 1;
        double prod3 = 1;
        for (unsigned int ljet=0; ljet<jetsize; ++ljet) {
	  //          if(isGoodJet(ljet,pTthresh)) continue;
	  if ( (ljet==ijet)||(ljet==jjet)||(ljet==kjet) ) continue; 
          double efflt = jetTagEff_ewk(event.at(ljet), kT); // Tight
          double efflm = jetTagEff_ewk(event.at(ljet), kM); // Medium 
          double effll = jetTagEff_ewk(event.at(ljet), kL); // Loose 
	  prod1 *= (efflm-efflt);
	  prod3 *= (effll-efflm);
	  double prod2 = 1;
	  double prod4 = 1;
          for (unsigned int mjet=0; mjet<jetsize; ++mjet) {
	    //            if(isGoodJet(mjet,pTthresh)) continue;
	    if ( (mjet==ijet)||(mjet==jjet)||(mjet==kjet)||(mjet==ljet) ) continue; 
            double effmm = jetTagEff_ewk(event.at(mjet), kM); // Medium 
            double effml = jetTagEff_ewk(event.at(mjet), kL); // Loose 
	    prod2 *= (1-effmm);
	    prod4 *= (1-effml);
          } //m loop
          Prob4b += effi*effj*effk*prod1*prod2;
          Prob4b += effi*effj*effk*prod3*prod4;
        } //l loop
      }// k loop
    }//j loop
  }//i loop

  //2T (M-T) (M-T) (M-T)
  for (unsigned int ijet=0; ijet<jetsize-1; ++ijet) {
    //    if(!isGoodJet(ijet,pTthresh)) continue; //switch to 20 GeV threshold for b jets
    double effi = jetTagEff_ewk(event.at(ijet), kT); // Tight
    for (unsigned int jjet=ijet+1; jjet<jetsize; ++jjet) {
      //      if(!isGoodJet(jjet,pTthresh)) continue;
      double effj = jetTagEff_ewk(event.at(jjet), kT); // Tight
      double prod = 1;
      for (unsigned int kjet=0; kjet<jetsize-2; ++kjet) {
	//        if(!isGoodJet(kjet,pTthresh)) continue;
	if ( (kjet==ijet)||(kjet==jjet) ) continue; 
        double effkt = jetTagEff_ewk(event.at(kjet), kT); // Tight
        double effkm = jetTagEff_ewk(event.at(kjet), kM); // Medium
        prod *= (effkm-effkt);
        double prod1 = 1;
        for (unsigned int ljet=kjet+1; ljet<jetsize-1; ++ljet) {
	  //          if(isGoodJet(ljet,pTthresh)) continue;
	  if ( (ljet==ijet)||(ljet==jjet) ) continue; 
          double efflt = jetTagEff_ewk(event.at(ljet), kT); // Tight
          double efflm = jetTagEff_ewk(event.at(ljet), kM); // Medium 
	  prod1 *= (efflm-efflt);
	  double prod2 = 1;
          for (unsigned int mjet=ljet+1; mjet<jetsize; ++mjet) {
	    //            if(isGoodJet(mjet,pTthresh)) continue;
	    if ( (mjet==ijet)||(mjet==jjet) ) continue; 
            double effmt = jetTagEff_ewk(event.at(mjet), kT); // Tight 
            double effmm = jetTagEff_ewk(event.at(mjet), kM); // Medium 
	    prod2 *= (effmm-effmt);
          } //m loop
          Prob4b += effi*effj*prod*prod1*prod2;
        } //l loop
      }// k loop
    }//j loop
  }//i loop

  //2T (M-T) (M-T) (1-M)
  for (unsigned int ijet=0; ijet<jetsize-1; ++ijet) {
//    if(!isGoodJet(ijet,pTthresh)) continue; //switch to 20 GeV threshold for b jets
    double effi = jetTagEff_ewk(event.at(ijet), kT); // Tight
    for (unsigned int jjet=ijet+1; jjet<jetsize; ++jjet) {
      //      if(!isGoodJet(jjet,pTthresh)) continue;
      double effj = jetTagEff_ewk(event.at(jjet), kT); // Tight
      double prod = 1;
      for (unsigned int kjet=0; kjet<jetsize-1; ++kjet) {
	//        if(!isGoodJet(kjet,pTthresh)) continue;
	if ( (kjet==ijet)||(kjet==jjet) ) continue; 
        double effkt = jetTagEff_ewk(event.at(kjet), kT); // Tight
        double effkm = jetTagEff_ewk(event.at(kjet), kM); // Medium
        prod *= (effkm-effkt);
        double prod1 = 1;
        for (unsigned int ljet=kjet+1; ljet<jetsize; ++ljet) {
	  //          if(isGoodJet(ljet,pTthresh)) continue;
	  if ( (ljet==ijet)||(ljet==jjet) ) continue; 
          double efflt = jetTagEff_ewk(event.at(ljet), kT); // Tight
          double efflm = jetTagEff_ewk(event.at(ljet), kM); // Medium 
	  prod1 *= (efflm-efflt);
	  double prod2 = 1;
          for (unsigned int mjet=0; mjet<jetsize; ++mjet) {
	    //            if(isGoodJet(mjet,pTthresh)) continue;
	    if ( (mjet==ijet)||(mjet==jjet)||(mjet==kjet)||(mjet==ljet) ) continue; 
            double effmm = jetTagEff_ewk(event.at(mjet), kM); // Medium 
	    prod2 *= (1-effmm);
          } //m loop
          Prob4b += effi*effj*prod*prod1*prod2;
        } //l loop
      }// k loop
    }//j loop
  }//i loop

  //2T (M-T) (L-M) (L-M)
  for (unsigned int ijet=0; ijet<jetsize-1; ++ijet) {
    //    if(!isGoodJet(ijet,pTthresh)) continue; //switch to 20 GeV threshold for b jets
    double effi = jetTagEff_ewk(event.at(ijet), kT); // Tight
    for (unsigned int jjet=ijet+1; jjet<jetsize; ++jjet) {
      //      if(!isGoodJet(jjet,pTthresh)) continue;
      double effj = jetTagEff_ewk(event.at(jjet), kT); // Tight
      double prod = 1;
      for (unsigned int kjet=0; kjet<jetsize; ++kjet) {
	//        if(!isGoodJet(kjet,pTthresh)) continue;
	if ( (kjet==ijet)||(kjet==jjet) ) continue; 
        double effkt = jetTagEff_ewk(event.at(kjet), kT); // Tight
        double effkm = jetTagEff_ewk(event.at(kjet), kM); // Medium
        prod *= (effkm-effkt);
        double prod1 = 1;
        for (unsigned int ljet=0; ljet<jetsize-1; ++ljet) {
	  //          if(isGoodJet(ljet,pTthresh)) continue;
	  if ( (ljet==ijet)||(ljet==jjet) ) continue; 
          double efflm = jetTagEff_ewk(event.at(ljet), kM); // Medium 
          double effll = jetTagEff_ewk(event.at(ljet), kL); // Loose 
	  prod1 *= (effll-efflm);
	  double prod2 = 1;
          for (unsigned int mjet=ljet+1; mjet<jetsize; ++mjet) {
	    //            if(isGoodJet(mjet,pTthresh)) continue;
	    if ( (mjet==ijet)||(mjet==jjet)||(mjet==kjet) ) continue; 
            double effmm = jetTagEff_ewk(event.at(mjet), kM); // Medium 
            double effml = jetTagEff_ewk(event.at(mjet), kL); // Loose 
	    prod2 *= (effml-effmm);
          } //m loop
          Prob4b += effi*effj*prod*prod1*prod2;
        } //l loop
      }// k loop
    }//j loop
  }//i loop

  //2T (M-T) (L-M) (1-L)
  for (unsigned int ijet=0; ijet<jetsize-1; ++ijet) {
    //    if(!isGoodJet(ijet,pTthresh)) continue; //switch to 20 GeV threshold for b jets
    double effi = jetTagEff_ewk(event.at(ijet), kT); // Tight
    for (unsigned int jjet=ijet+1; jjet<jetsize; ++jjet) {
      //      if(!isGoodJet(jjet,pTthresh)) continue;
      double effj = jetTagEff_ewk(event.at(jjet), kT); // Tight
      double prod = 1;
      for (unsigned int kjet=0; kjet<jetsize; ++kjet) {
	//        if(!isGoodJet(kjet,pTthresh)) continue;
	if ( (kjet==ijet)||(kjet==jjet) ) continue; 
        double effkt = jetTagEff_ewk(event.at(kjet), kT); // Tight
        double effkm = jetTagEff_ewk(event.at(kjet), kM); // Medium
        prod *= (effkm-effkt);
        double prod1 = 1;
        for (unsigned int ljet=0; ljet<jetsize; ++ljet) {
	  //          if(isGoodJet(ljet,pTthresh)) continue;
	  if ( (ljet==ijet)||(ljet==jjet)||(ljet==kjet) ) continue; 
          double efflm = jetTagEff_ewk(event.at(ljet), kM); // Medium 
          double effll = jetTagEff_ewk(event.at(ljet), kL); // Loose 
	  prod1 *= (effll-efflm);
	  double prod2 = 1;
          for (unsigned int mjet=0; mjet<jetsize; ++mjet) {
	    //            if(isGoodJet(mjet,pTthresh)) continue;
	    if ( (mjet==ijet)||(mjet==jjet)||(mjet==kjet)||(mjet==ljet) ) continue; 
            double effml = jetTagEff_ewk(event.at(mjet), kL); // Loose 
	    prod2 *= (1-effml);
          } //m loop
          Prob4b += effi*effj*prod*prod1*prod2;
        } //l loop
      }// k loop
    }//j loop
  }//i loop

  } //5jets or more 

}//end method


void countWeighted( const vector< vector<pair<jetFlavor,tagStatus> >  > & dataset,double & n2b,double & n3b,double & n4b) {


  n2b=0;
  n3b=0;
  n4b=0;

  double  n2be=0;
  double  n3be=0;
  double  n4be=0;

  //loop over events
  for (size_t iev = 0; iev<dataset.size(); ++iev) {
    vector<pair<jetFlavor,tagStatus> > event = dataset.at(iev);
    float prob2,prob3,prob4;
    calculateTagProb_ewk(event, prob2,prob3,prob4);
    n2b+= prob2;
    n3b+= prob3;
    n4b+=prob4;
    //sum of the squares of the weights
    n2be+= prob2*prob2;
    n3be+= prob3*prob3;
    n4be+=prob4*prob4;
  }

  cout<<"2b 3b 4b  "<<n2b<<" +- "<<sqrt(n2be)<<" "<<n3b<<" +- "<<sqrt(n3be)<<" "<<n4b<<" +- "<<sqrt(n4be)<<endl;


}

void PJgo(const ULong64_t ngen = 1000000) {
  rnd_ = new TRandom3(223449876); //some seed
  closureTestMode_=true; //set SFs to 1
  chanceOfJet5_=0.5; //turn off 5th jet


  vector< vector<pair<jetFlavor,tagStatus> >  > MC;
  //  vector< vector<pair<jetFlavor,tagStatus> >  > data;

  for (ULong64_t iev = 0; iev <ngen; ++iev) {
    //    vector<pair<jetFlavor,tagStatus> > dataevent  = generateTtbar(false);
    //    data.push_back(dataevent);

    vector<pair<jetFlavor,tagStatus> > mcevent  = generateTtbar(true);
    MC.push_back(mcevent);
  }

  double n2b_mc,n3b_mc,n4b_mc;
  double n2b_data,n3b_data,n4b_data;
  cout<<" == MC   =="<<endl;
  countInBins(MC,n2b_mc,n3b_mc,n4b_mc);
  cout<<" == MC after weights  =="<<endl;
  countWeighted(MC, n2b_mc,n3b_mc,n4b_mc);
  //  cout<<" == data =="<<endl;
  //  countInBins(data,n2b_data,n3b_data,n4b_data);


}
