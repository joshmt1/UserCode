//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 19 14:02:26 2010 by ROOT version 5.22/00d
// from TChain makeTopologyNtuple/tree/
//////////////////////////////////////////////////////////

#ifndef eventLoop_h
#define eventLoop_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <TH1D.h>

#include <iostream>

class eventLoop {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           numEle;
   Float_t         eleE[3];   //[numEle]
   Float_t         eleET[3];   //[numEle]
   Float_t         elePX[3];   //[numEle]
   Float_t         elePY[3];   //[numEle]
   Float_t         elePZ[3];   //[numEle]
   Float_t         elePhi[3];   //[numEle]
   Float_t         eleTheta[3];   //[numEle]
   Float_t         eleEta[3];   //[numEle]
   Int_t           eleCharge[3];   //[numEle]
   Int_t           eleIDQuality[3];   //[numEle]
   Float_t         eleTrackPt[3];   //[numEle]
   Float_t         eleTrackPhi[3];   //[numEle]
   Float_t         eleTrackEta[3];   //[numEle]
   Float_t         eleTrackChi2[3];   //[numEle]
   Float_t         eleTrackNDOF[3];   //[numEle]
   Float_t         eleTrackD0[3];   //[numEle]
   Float_t         eleBeamSpotCorrectedTrackD0[3];   //[numEle]
   Float_t         eleTrackDz[3];   //[numEle]
   Float_t         eleSCEta[3];   //[numEle]
   Float_t         eleSCE[3];   //[numEle]
   Float_t         eleSCPhi[3];   //[numEle]
   Float_t         eleSCSigmaEtaEta[3];   //[numEle]
   Float_t         eleSCSigmaIEtaIEta[3];   //[numEle]
   Float_t         eleSCE1x5[3];   //[numEle]
   Float_t         eleSCE5x5[3];   //[numEle]
   Float_t         eleSCE2x5max[3];   //[numEle]
   Float_t         eleTrackIso[3];   //[numEle]
   Float_t         eleEcalIso[3];   //[numEle]
   Float_t         eledr04EcalRecHitSumEt[3];   //[numEle]
   Float_t         eledr03EcalRecHitSumEt[3];   //[numEle]
   Float_t         eleHcalIso[3];   //[numEle]
   Float_t         eleEcalIsoDeposit[3];   //[numEle]
   Float_t         eleHcalIsoDeposit[3];   //[numEle]
   Float_t         eleComRelIso[3];   //[numEle]
   Int_t           elePhotonConversionTag[3];   //[numEle]
   Float_t         elePhotonConversionDist[3];   //[numEle]
   Float_t         elePhotonConversionDcot[3];   //[numEle]
   Float_t         eleTriggerMatch[3];   //[numEle]
   Float_t         eleJetOverlap[3];   //[numEle]
   Float_t         genEleET[3];   //[numEle]
   Float_t         genElePX[3];   //[numEle]
   Float_t         genElePY[3];   //[numEle]
   Float_t         genElePZ[3];   //[numEle]
   Float_t         genElePhi[3];   //[numEle]
   Float_t         genEleTheta[3];   //[numEle]
   Float_t         genEleEta[3];   //[numEle]
   Int_t           genEleCharge[3];   //[numEle]
   Int_t           numMuo;
   Float_t         muoE[3];   //[numMuo]
   Float_t         muoET[3];   //[numMuo]
   Float_t         muoPt[3];   //[numMuo]
   Float_t         muoPX[3];   //[numMuo]
   Float_t         muoPY[3];   //[numMuo]
   Float_t         muoPZ[3];   //[numMuo]
   Float_t         muoPhi[3];   //[numMuo]
   Float_t         muoTheta[3];   //[numMuo]
   Float_t         muoEta[3];   //[numMuo]
   Int_t           muoCharge[3];   //[numMuo]
   Float_t         muonChi2[3];   //[numMuo]
   Float_t         muonD0[3];   //[numMuo]
   Float_t         muonBeamSpotCorrectedD0[3];   //[numMuo]
   Int_t           muonTrackNHits[3];   //[numMuo]
   Float_t         muonNDOF[3];   //[numMuo]
   Float_t         muonTrackIso[3];   //[numMuo]
   Float_t         muonEcalIso[3];   //[numMuo]
   Float_t         muonHcalIso[3];   //[numMuo]
   Float_t         muonComRelIso[3];   //[numMuo]
   Float_t         genMuoET[3];   //[numMuo]
   Float_t         genMuoPX[3];   //[numMuo]
   Float_t         genMuoPY[3];   //[numMuo]
   Float_t         genMuoPZ[3];   //[numMuo]
   Float_t         genMuoPhi[3];   //[numMuo]
   Float_t         genMuoTheta[3];   //[numMuo]
   Float_t         genMuoEta[3];   //[numMuo]
   Int_t           genMuoCharge[3];   //[numMuo]
   Int_t           nT;
   Int_t           nThadronic;
   Float_t         T_hadronicMCTruthE[3];   //[nThadronic]
   Float_t         T_hadronicMCTruthEt[3];   //[nThadronic]
   Float_t         T_hadronicMCTruthPx[3];   //[nThadronic]
   Float_t         T_hadronicMCTruthPy[3];   //[nThadronic]
   Float_t         T_hadronicMCTruthPz[3];   //[nThadronic]
   Int_t           T_hadronicMCMotherIndex[3];   //[nThadronic]
   Int_t           nTleptonic;
   Float_t         T_leptonicMCTruthE[2];   //[nTleptonic]
   Float_t         T_leptonicMCTruthEt[2];   //[nTleptonic]
   Float_t         T_leptonicMCTruthPx[2];   //[nTleptonic]
   Float_t         T_leptonicMCTruthPy[2];   //[nTleptonic]
   Float_t         T_leptonicMCTruthPz[2];   //[nTleptonic]
   Int_t           T_leptonicMCMotherIndex[2];   //[nTleptonic]
   Int_t           nb;
   Float_t         bMCTruthE[4];   //[nb]
   Float_t         bMCTruthEt[4];   //[nb]
   Float_t         bMCTruthPx[4];   //[nb]
   Float_t         bMCTruthPy[4];   //[nb]
   Float_t         bMCTruthPz[4];   //[nb]
   Int_t           bMCTruthMother[4];   //[nb]
   Int_t           nWhadronic;
   Float_t         W_hadronicMCTruthE[4];   //[nWhadronic]
   Float_t         W_hadronicMCTruthEt[4];   //[nWhadronic]
   Float_t         W_hadronicMCTruthPx[4];   //[nWhadronic]
   Float_t         W_hadronicMCTruthPy[4];   //[nWhadronic]
   Float_t         W_hadronicMCTruthPz[4];   //[nWhadronic]
   Int_t           W_hadronicMCTruthPID[4];   //[nWhadronic]
   Int_t           W_hadronicMCTruthMother[4];   //[nWhadronic]
   Int_t           nWleptonic;
   Float_t         W_leptonicMCTruthE[3];   //[nWleptonic]
   Float_t         W_leptonicMCTruthEt[3];   //[nWleptonic]
   Float_t         W_leptonicMCTruthPx[3];   //[nWleptonic]
   Float_t         W_leptonicMCTruthPy[3];   //[nWleptonic]
   Float_t         W_leptonicMCTruthPz[3];   //[nWleptonic]
   Int_t           W_leptonicMCTruthPID[3];   //[nWleptonic]
   Int_t           W_leptonicMCTruthMother[3];   //[nWleptonic]
   Int_t           isElePlusJets;
   Int_t           VQQBosonAbsId;
   Int_t           numJet;
   Float_t         jetE[15];   //[numJet]
   Float_t         jetEt[15];   //[numJet]
   Float_t         jetPt[15];   //[numJet]
   Float_t         jetCorEt[15];   //[numJet]
   Float_t         jetEta[15];   //[numJet]
   Float_t         jetTheta[15];   //[numJet]
   Float_t         jetPhi[15];   //[numJet]
   Float_t         jetPx[15];   //[numJet]
   Float_t         jetPy[15];   //[numJet]
   Float_t         jetPz[15];   //[numJet]
   Int_t           jetNtracksInJet[15];   //[numJet]
   Float_t         jetJetCharge[15];   //[numJet]
   Float_t         jetMuEnergy[15];   //[numJet]
   Float_t         jetMuEnergyFraction[15];   //[numJet]
   Float_t         jetChargedMultiplicity[15];   //[numJet]
   Float_t         jetNeutralHadEnergy[15];   //[numJet]
   Float_t         jetEMEnergyInEB[15];   //[numJet]
   Float_t         jetEMEnergyInEE[15];   //[numJet]
   Float_t         jetEMEnergyFraction[15];   //[numJet]
   Float_t         jetEMEnergyInHF[15];   //[numJet]
   Float_t         jetHadEnergyInHB[15];   //[numJet]
   Float_t         jetHadEnergyInHE[15];   //[numJet]
   Float_t         jetHadEnergyInHF[15];   //[numJet]
   Float_t         jetHadEnergyInHO[15];   //[numJet]
   Float_t         jetBtagTrackCountHighPurity[15];   //[numJet]
   Float_t         jetBtagTrackCountHighEff[15];   //[numJet]
   Float_t         jetBtagProbability[15];   //[numJet]
   Float_t         jetBtagSoftElectron[15];   //[numJet]
   Float_t         jetBtagSoftMuon[15];   //[numJet]
   Float_t         jetBtagSoftMuonNoIP[15];   //[numJet]
   Float_t         jetBtagSoftMuonPtRel[15];   //[numJet]
   Float_t         jetBtagSoftMuonQuality[15];   //[numJet]
   Float_t         jetBtagSecondaryVertex[15];   //[numJet]
   Float_t         jetBtagSecondaryVertexNegative[15];   //[numJet]
   Float_t         jetBtagCombinedSVLL[15];   //[numJet]
   Float_t         jetBtagCombinedSVMVA[15];   //[numJet]
   Float_t         jetCorrFactor[15];   //[numJet]
   Float_t         jetN60[15];   //[numJet]
   Float_t         jetN90[15];   //[numJet]
   Float_t         jetNeutralEmEnergy[15];   //[numJet]
   Float_t         jetTriggered[15];   //[numJet]
   Float_t         jetSVPT[15];   //[numJet]
   Float_t         jetSVL2D[15];   //[numJet]
   Float_t         jetSVL2Dxy[15];   //[numJet]
   Float_t         jetSVL2DxyErr[15];   //[numJet]
   Float_t         jetSVL2DxySig[15];   //[numJet]
   Float_t         jetSVL3D[15];   //[numJet]
   Float_t         jetSVL3DErr[15];   //[numJet]
   Float_t         jetSVL3DSig[15];   //[numJet]
   Float_t         jetSVMass[15];   //[numJet]
   Int_t           jetSVNtracks[15];   //[numJet]
   Float_t         genJetET[15];   //[numJet]
   Float_t         genJetPX[15];   //[numJet]
   Float_t         genJetPY[15];   //[numJet]
   Float_t         genJetPZ[15];   //[numJet]
   Float_t         genJetPhi[15];   //[numJet]
   Float_t         genJetTheta[15];   //[numJet]
   Float_t         genJetEta[15];   //[numJet]
   Int_t           genJetPID[15];   //[numJet]
   Int_t           jetPID[15];   //[numJet]
   Float_t         jetClosestBPartonDeltaR[15];   //[numJet]
   Float_t         jetClosestCPartonDeltaR[15];   //[numJet]
   Float_t         btagParamDiscCut_MISTAGSSVM;
   Float_t         btagParamDiscCut_PTRELSSVM;
   Float_t         btagParamDiscCut_PTRELTCHEL;
   Float_t         btagParamDiscCut_PTRELTCHEM;
   Float_t         btagParamDiscCut_PTRELTCHET;
   Float_t         btagParamDiscCut_PTRELTCHPL;
   Float_t         btagParamDiscCut_PTRELTCHPM;
   Float_t         btagParamDiscCut_PTRELTCHPT;
   Float_t         btagParamDiscCut_SYSTEM8SSVM;
   Float_t         jetBtagParam_MISTAGSSVM_BTAGLEFF[15];   //[numJet]
   Float_t         jetBtagParam_MISTAGSSVM_BTAGLERR[15];   //[numJet]
   Float_t         jetBtagParam_MISTAGSSVM_BTAGLEFFCORR[15];   //[numJet]
   Float_t         jetBtagParam_MISTAGSSVM_BTAGLERRCORR[15];   //[numJet]
   Float_t         jetBtagParam_PTRELSSVM_BTAGBEFFCORR[15];   //[numJet]
   Float_t         jetBtagParam_PTRELSSVM_BTAGBERRCORR[15];   //[numJet]
   Float_t         jetBtagParam_SYSTEM8SSVM_BTAGBEFF[15];   //[numJet]
   Float_t         jetBtagParam_SYSTEM8SSVM_BTAGBERR[15];   //[numJet]
   Float_t         jetBtagParam_PTRELTCHEL_BTAGBEFFCORR[15];   //[numJet]
   Float_t         jetBtagParam_PTRELTCHEM_BTAGBEFFCORR[15];   //[numJet]
   Float_t         jetBtagParam_PTRELTCHET_BTAGBEFFCORR[15];   //[numJet]
   Float_t         jetBtagParam_PTRELTCHPL_BTAGBEFFCORR[15];   //[numJet]
   Float_t         jetBtagParam_PTRELTCHPM_BTAGBEFFCORR[15];   //[numJet]
   Float_t         jetBtagParam_PTRELTCHPT_BTAGBEFFCORR[15];   //[numJet]
   Int_t           numGeneralTracks;
   Float_t         generalTracksPt[456];   //[numGeneralTracks]
   Float_t         generalTracksEta[456];   //[numGeneralTracks]
   Float_t         generalTracksTheta[456];   //[numGeneralTracks]
   Float_t         generalTracksBeamSpotCorrectedD0[456];   //[numGeneralTracks]
   Float_t         generalTracksPhi[456];   //[numGeneralTracks]
   Int_t           generalTracksCharge[456];   //[numGeneralTracks]
   Int_t           processId;
   Float_t         processPtHat;
   Float_t         processMCWeight;
   Float_t         beamSpotX;
   Float_t         beamSpotY;
   Float_t         beamSpotZ;
   Float_t         topo_sphericity;
   Float_t         topo_aplanarity;
   Float_t         topo_sphericity_e;
   Float_t         topo_aplanarity_e;
   Float_t         topo_ht;
   Float_t         topo_ht_e;
   Float_t         topo_sqrts;
   Float_t         topo_sqrts_e;
   Float_t         topo_oblateness;
   Float_t         metEt;
   Float_t         metPhi;
   Float_t         metPt;
   Float_t         metPx;
   Float_t         metPy;
   Float_t         genMetEt;
   Float_t         genMetPhi;
   Float_t         genMetPt;
   Float_t         genMetPx;
   Float_t         genMetPy;
   Float_t         metMaxEtEM;
   Float_t         metMaxEtHad;
   Float_t         metEtFracHad;
   Float_t         metEtFracEM;
   Float_t         metHadEtHB;
   Float_t         metHadEtHO;
   Float_t         metHadEtHF;
   Float_t         metHadEtHE;
   Float_t         metEmEtHF;
   Float_t         metEmEtEE;
   Float_t         metEmEtEB;
   Float_t         metSignificance;
   Float_t         metScalarEt;
   Float_t         metEtUncorrected;
   Float_t         metPhiUncorrected;
   Float_t         mhtPt;
   Float_t         mhtPy;
   Float_t         mhtPx;
   Float_t         mhtPhi;
   Float_t         mhtSumEt;
   Float_t         mhtSignif;
   Int_t           nZCandidates;
   Float_t         ZCandidates[3];   //[nZCandidates]
   Int_t           nTriggerBits;
   Int_t           TriggerBits[79];   //[nTriggerBits]
   Int_t           HLT_Ele15_LW_L1R;
   Int_t           flavorhistory;
   Int_t           myProcess;
   Float_t         genParEta[2];
   Float_t         genParPhi[2];
   Float_t         genParTheta[2];
   Float_t         genParE[2];
   Float_t         genParEt[2];
   Float_t         genParP[2];
   Float_t         genParPt[2];
   Int_t           genParId[2];
   Int_t           genParStat[2];
   Int_t           ndaughters[2];
   Int_t           genParCharge[2];
   Float_t         genEta[2][2];
   Float_t         genPhi[2][2];
   Float_t         genTheta[2][2];
   Float_t         genE[2][2];
   Float_t         genEt[2][2];
   Float_t         genP[2][2];
   Float_t         genPt[2][2];
   Int_t           genId[2][2];
   Int_t           genStat[2][2];
   Int_t           genCharge[2][2];
   Int_t           numTau;
   Float_t         tauE[5];   //[numTau]
   Float_t         tauPt[5];   //[numTau]
   Float_t         tauPhi[5];   //[numTau]
   Float_t         tauEta[5];   //[numTau]
   Int_t           numPhoton;
   Float_t         photonE[13];   //[numPhoton]
   Float_t         photonPt[13];   //[numPhoton]
   Float_t         photonPhi[13];   //[numPhoton]
   Float_t         photonEta[13];   //[numPhoton]
   Int_t           eventRun;
   Int_t           eventNum;
   Float_t         eventLumiblock;

   // List of branches
   TBranch        *b_numEle;   //!
   TBranch        *b_eleE;   //!
   TBranch        *b_eleET;   //!
   TBranch        *b_elePX;   //!
   TBranch        *b_elePY;   //!
   TBranch        *b_elePZ;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleTheta;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleIDQuality;   //!
   TBranch        *b_eleTrackPt;   //!
   TBranch        *b_eleTrackPhi;   //!
   TBranch        *b_eleTrackEta;   //!
   TBranch        *b_eleTrackChi2;   //!
   TBranch        *b_eleTrackNDOF;   //!
   TBranch        *b_eleTrackD0;   //!
   TBranch        *b_eleBeamSpotCorrectedTrackD0;   //!
   TBranch        *b_eleTrackDz;   //!
   TBranch        *b_eleSCEta;   //!
   TBranch        *b_eleSCE;   //!
   TBranch        *b_eleSCPhi;   //!
   TBranch        *b_eleSCSigmaEtaEta;   //!
   TBranch        *b_eleSCSigmaIEtaIEta;   //!
   TBranch        *b_eleSCE1x5;   //!
   TBranch        *b_eleSCE5x5;   //!
   TBranch        *b_eleSCE2x5max;   //!
   TBranch        *b_eleTrackIso;   //!
   TBranch        *b_eleEcalIso;   //!
   TBranch        *b_eledr04EcalRecHitSumEt;   //!
   TBranch        *b_eledr03EcalRecHitSumEt;   //!
   TBranch        *b_eleHcalIso;   //!
   TBranch        *b_eleEcalIsoDeposit;   //!
   TBranch        *b_eleHcalIsoDeposit;   //!
   TBranch        *b_eleComRelIso;   //!
   TBranch        *b_elePhotonConversionTag;   //!
   TBranch        *b_elePhotonConversionDist;   //!
   TBranch        *b_elePhotonConversionDcot;   //!
   TBranch        *b_eleTriggerMatch;   //!
   TBranch        *b_eleJetOverlap;   //!
   TBranch        *b_genEleET;   //!
   TBranch        *b_genElePX;   //!
   TBranch        *b_genElePY;   //!
   TBranch        *b_genElePZ;   //!
   TBranch        *b_genElePhi;   //!
   TBranch        *b_genEleTheta;   //!
   TBranch        *b_genEleEta;   //!
   TBranch        *b_genEleCharge;   //!
   TBranch        *b_numMuo;   //!
   TBranch        *b_muoE;   //!
   TBranch        *b_muoET;   //!
   TBranch        *b_muoPt;   //!
   TBranch        *b_muoPX;   //!
   TBranch        *b_muoPY;   //!
   TBranch        *b_muoPZ;   //!
   TBranch        *b_muoPhi;   //!
   TBranch        *b_muoTheta;   //!
   TBranch        *b_muoEta;   //!
   TBranch        *b_muoCharge;   //!
   TBranch        *b_muonChi2;   //!
   TBranch        *b_muonD0;   //!
   TBranch        *b_muonBeamSpotCorrectedD0;   //!
   TBranch        *b_muonTrackNHits;   //!
   TBranch        *b_muonNDOF;   //!
   TBranch        *b_muonTrackIso;   //!
   TBranch        *b_muonEcalIso;   //!
   TBranch        *b_muonHcalIso;   //!
   TBranch        *b_muonComRelIso;   //!
   TBranch        *b_genMuoET;   //!
   TBranch        *b_genMuoPX;   //!
   TBranch        *b_genMuoPY;   //!
   TBranch        *b_genMuoPZ;   //!
   TBranch        *b_genMuoPhi;   //!
   TBranch        *b_genMuoTheta;   //!
   TBranch        *b_genMuoEta;   //!
   TBranch        *b_genMuoCharge;   //!
   TBranch        *b_nT;   //!
   TBranch        *b_nThadronic;   //!
   TBranch        *b_T_hadronicMCTruthE;   //!
   TBranch        *b_T_hadronicMCTruthEt;   //!
   TBranch        *b_T_hadronicMCTruthPx;   //!
   TBranch        *b_T_hadronicMCTruthPy;   //!
   TBranch        *b_T_hadronicMCTruthPz;   //!
   TBranch        *b_T_hadronicMCMotherIndex;   //!
   TBranch        *b_nTleptonic;   //!
   TBranch        *b_T_leptonicMCTruthE;   //!
   TBranch        *b_T_leptonicMCTruthEt;   //!
   TBranch        *b_T_leptonicMCTruthPx;   //!
   TBranch        *b_T_leptonicMCTruthPy;   //!
   TBranch        *b_T_leptonicMCTruthPz;   //!
   TBranch        *b_T_leptonicMCMotherIndex;   //!
   TBranch        *b_nb;   //!
   TBranch        *b_bMCTruthE;   //!
   TBranch        *b_bMCTruthEt;   //!
   TBranch        *b_bMCTruthPx;   //!
   TBranch        *b_bMCTruthPy;   //!
   TBranch        *b_bMCTruthPz;   //!
   TBranch        *b_bMCTruthMother;   //!
   TBranch        *b_nWhadronic;   //!
   TBranch        *b_W_hadronicMCTruthE;   //!
   TBranch        *b_W_hadronicMCTruthEt;   //!
   TBranch        *b_W_hadronicMCTruthPx;   //!
   TBranch        *b_W_hadronicMCTruthPy;   //!
   TBranch        *b_W_hadronicMCTruthPz;   //!
   TBranch        *b_W_hadronicMCTruthPID;   //!
   TBranch        *b_W_hadronicMCTruthMother;   //!
   TBranch        *b_nWleptonic;   //!
   TBranch        *b_W_leptonicMCTruthE;   //!
   TBranch        *b_W_leptonicMCTruthEt;   //!
   TBranch        *b_W_leptonicMCTruthPx;   //!
   TBranch        *b_W_leptonicMCTruthPy;   //!
   TBranch        *b_W_leptonicMCTruthPz;   //!
   TBranch        *b_W_leptonicMCTruthPID;   //!
   TBranch        *b_W_leptonicMCTruthMother;   //!
   TBranch        *b_isElePlusJets;   //!
   TBranch        *b_VQQBosonAbsId;   //!
   TBranch        *b_numJet;   //!
   TBranch        *b_jetE;   //!
   TBranch        *b_jetEt;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetCorEt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetTheta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetPx;   //!
   TBranch        *b_jetPy;   //!
   TBranch        *b_jetPz;   //!
   TBranch        *b_jetNtracksInJet;   //!
   TBranch        *b_jetJetCharge;   //!
   TBranch        *b_jetMuEnergy;   //!
   TBranch        *b_jetMuEnergyFraction;   //!
   TBranch        *b_jetChargedMultiplicity;   //!
   TBranch        *b_jetNeutralHadEnergy;   //!
   TBranch        *b_jetEMEnergyInEB;   //!
   TBranch        *b_jetEMEnergyInEE;   //!
   TBranch        *b_jetEMEnergyFraction;   //!
   TBranch        *b_jetEMEnergyInHF;   //!
   TBranch        *b_jetHadEnergyInHB;   //!
   TBranch        *b_jetHadEnergyInHE;   //!
   TBranch        *b_jetHadEnergyInHF;   //!
   TBranch        *b_jetHadEnergyInHO;   //!
   TBranch        *b_jetBtagTrackCountHighPurity;   //!
   TBranch        *b_jetBtagTrackCountHighEff;   //!
   TBranch        *b_jetBtagProbability;   //!
   TBranch        *b_jetBtagSoftElectron;   //!
   TBranch        *b_jetBtagSoftMuon;   //!
   TBranch        *b_jetBtagSoftMuonNoIP;   //!
   TBranch        *b_jetBtagSoftMuonPtRel;   //!
   TBranch        *b_jetBtagSoftMuonQuality;   //!
   TBranch        *b_jetBtagSecondaryVertex;   //!
   TBranch        *b_jetBtagSecondaryVertexNegative;   //!
   TBranch        *b_jetBtagCombinedSVLL;   //!
   TBranch        *b_jetBtagCombinedSVMVA;   //!
   TBranch        *b_jetCorrFactor;   //!
   TBranch        *b_jetN60;   //!
   TBranch        *b_jetN90;   //!
   TBranch        *b_jetNeutralEmEnergy;   //!
   TBranch        *b_jetTriggered;   //!
   TBranch        *b_jetSVPT;   //!
   TBranch        *b_jetSVL2D;   //!
   TBranch        *b_jetSVL2Dxy;   //!
   TBranch        *b_jetSVL2DxyErr;   //!
   TBranch        *b_jetSVL2DxySig;   //!
   TBranch        *b_jetSVL3D;   //!
   TBranch        *b_jetSVL3DErr;   //!
   TBranch        *b_jetSVL3DSig;   //!
   TBranch        *b_jetSVMass;   //!
   TBranch        *b_jetSVNtracks;   //!
   TBranch        *b_genJetET;   //!
   TBranch        *b_genJetPX;   //!
   TBranch        *b_genJetPY;   //!
   TBranch        *b_genJetPZ;   //!
   TBranch        *b_genJetPhi;   //!
   TBranch        *b_genJetTheta;   //!
   TBranch        *b_genJetEta;   //!
   TBranch        *b_genJetPID;   //!
   TBranch        *b_jetPID;   //!
   TBranch        *b_jetClosestBPartonDeltaR;   //!
   TBranch        *b_jetClosestCPartonDeltaR;   //!
   TBranch        *b_btagParamDiscCut_MISTAGSSVM;   //!
   TBranch        *b_btagParamDiscCut_PTRELSSVM;   //!
   TBranch        *b_btagParamDiscCut_PTRELTCHEL;   //!
   TBranch        *b_btagParamDiscCut_PTRELTCHEM;   //!
   TBranch        *b_btagParamDiscCut_PTRELTCHET;   //!
   TBranch        *b_btagParamDiscCut_PTRELTCHPL;   //!
   TBranch        *b_btagParamDiscCut_PTRELTCHPM;   //!
   TBranch        *b_btagParamDiscCut_PTRELTCHPT;   //!
   TBranch        *b_btagParamDiscCut_SYSTEM8SSVM;   //!
   TBranch        *b_jetBtagParam_MISTAGSSVM_BTAGLEFF;   //!
   TBranch        *b_jetBtagParam_MISTAGSSVM_BTAGLERR;   //!
   TBranch        *b_jetBtagParam_MISTAGSSVM_BTAGLEFFCORR;   //!
   TBranch        *b_jetBtagParam_MISTAGSSVM_BTAGLERRCORR;   //!
   TBranch        *b_jetBtagParam_PTRELSSVM_BTAGBEFFCORR;   //!
   TBranch        *b_jetBtagParam_PTRELSSVM_BTAGBERRCORR;   //!
   TBranch        *b_jetBtagParam_SYSTEM8SSVM_BTAGBEFF;   //!
   TBranch        *b_jetBtagParam_SYSTEM8SSVM_BTAGBERR;   //!
   TBranch        *b_jetBtagParam_PTRELTCHEL_BTAGBEFFCORR;   //!
   TBranch        *b_jetBtagParam_PTRELTCHEM_BTAGBEFFCORR;   //!
   TBranch        *b_jetBtagParam_PTRELTCHET_BTAGBEFFCORR;   //!
   TBranch        *b_jetBtagParam_PTRELTCHPL_BTAGBEFFCORR;   //!
   TBranch        *b_jetBtagParam_PTRELTCHPM_BTAGBEFFCORR;   //!
   TBranch        *b_jetBtagParam_PTRELTCHPT_BTAGBEFFCORR;   //!
   TBranch        *b_numGeneralTracks;   //!
   TBranch        *b_generalTracksPt;   //!
   TBranch        *b_generalTracksEta;   //!
   TBranch        *b_generalTracksTheta;   //!
   TBranch        *b_generalTracksBeamSpotCorrectedD0;   //!
   TBranch        *b_generalTracksPhi;   //!
   TBranch        *b_generalTracksCharge;   //!
   TBranch        *b_processId;   //!
   TBranch        *b_processPtHat;   //!
   TBranch        *b_processMCWeight;   //!
   TBranch        *b_beamSpotX;   //!
   TBranch        *b_beamSpotY;   //!
   TBranch        *b_beamSpotZ;   //!
   TBranch        *b_topo_sphericity;   //!
   TBranch        *b_topo_aplanarity;   //!
   TBranch        *b_topo_sphericity_e;   //!
   TBranch        *b_topo_aplanarity_e;   //!
   TBranch        *b_topo_ht;   //!
   TBranch        *b_top_ht_e;   //!
   TBranch        *b_topo_sqrts;   //!
   TBranch        *b_top_sqrts_e;   //!
   TBranch        *b_topo_oblateness;   //!
   TBranch        *b_metEt;   //!
   TBranch        *b_metPhi;   //!
   TBranch        *b_metPt;   //!
   TBranch        *b_metPx;   //!
   TBranch        *b_metPy;   //!
   TBranch        *b_genMetEt;   //!
   TBranch        *b_genMetPhi;   //!
   TBranch        *b_genMetPt;   //!
   TBranch        *b_genMetPx;   //!
   TBranch        *b_genMetPy;   //!
   TBranch        *b_metMaxEtEM;   //!
   TBranch        *b_metMaxEtHad;   //!
   TBranch        *b_metEtFracHad;   //!
   TBranch        *b_metEtFracEM;   //!
   TBranch        *b_metHadEtHB;   //!
   TBranch        *b_metHadEtHO;   //!
   TBranch        *b_metHadEtHF;   //!
   TBranch        *b_metHadEtHE;   //!
   TBranch        *b_metEmEtHF;   //!
   TBranch        *b_metEmEtEE;   //!
   TBranch        *b_metEmEtEB;   //!
   TBranch        *b_metSignificance;   //!
   TBranch        *b_metScalarEt;   //!
   TBranch        *b_metEtUncorrected;   //!
   TBranch        *b_metPhiUncorrected;   //!
   TBranch        *b_mhtPt;   //!
   TBranch        *b_mhtPy;   //!
   TBranch        *b_mhtPx;   //!
   TBranch        *b_mhtPhi;   //!
   TBranch        *b_mhtSumEt;   //!
   TBranch        *b_mhtSignif;   //!
   TBranch        *b_nZCandidates;   //!
   TBranch        *b_ZCandidates;   //!
   TBranch        *b_nTriggerBits;   //!
   TBranch        *b_TriggerBits;   //!
   TBranch        *b_HLT_Ele15_LW_L1R;   //!
   TBranch        *b_flavorhistory;   //!
   TBranch        *b_myProcess;   //!
   TBranch        *b_genParEta;   //!
   TBranch        *b_genParPhi;   //!
   TBranch        *b_genParTheta;   //!
   TBranch        *b_genParE;   //!
   TBranch        *b_genParEt;   //!
   TBranch        *b_genParP;   //!
   TBranch        *b_genParPt;   //!
   TBranch        *b_genParId;   //!
   TBranch        *b_genParStat;   //!
   TBranch        *b_ndaughters;   //!
   TBranch        *b_genParCharge;   //!
   TBranch        *b_genEta;   //!
   TBranch        *b_genPhi;   //!
   TBranch        *b_genTheta;   //!
   TBranch        *b_genE;   //!
   TBranch        *b_genEt;   //!
   TBranch        *b_genP;   //!
   TBranch        *b_genPt;   //!
   TBranch        *b_genId;   //!
   TBranch        *b_genStat;   //!
   TBranch        *b_genCharge;   //!
   TBranch        *b_numTau;   //!
   TBranch        *b_tauE;   //!
   TBranch        *b_tauPt;   //!
   TBranch        *b_tauPhi;   //!
   TBranch        *b_tauEta;   //!
   TBranch        *b_numPhoton;   //!
   TBranch        *b_photonE;   //!
   TBranch        *b_photonPt;   //!
   TBranch        *b_photonPhi;   //!
   TBranch        *b_photonEta;   //!
   TBranch        *b_eventRun;   //!
   TBranch        *b_eventNum;   //!
   TBranch        *b_eventLumiblock;   //!

   //for changing configuration
   bool writeASCII_;

   //variables for defining cuts!
   UInt_t minGoodJets_;
   UInt_t goodBJets_;
   Float_t minJetPt_Good_;
   Float_t maxJetEta_Good_;

   eventLoop(TTree *tree=0);
   virtual ~eventLoop();
   virtual Int_t    Cut();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString outfile,TString weightfile="",bool isSignal=false);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void writeASCII(bool writeASCII) {writeASCII_=writeASCII;}
   Bool_t isJetGood(UInt_t index) ;
   UInt_t nGoodJets();
   UInt_t nGoodElectrons();
   UInt_t nGoodMuons();
   double getDeltaPhi(double phi1,double phi2);
};

#endif

#ifdef eventLoop_cxx
eventLoop::eventLoop(TTree *tree) :
  writeASCII_(false),
  //define default values of selection cuts
  minGoodJets_(3),
  goodBJets_(0),
  minJetPt_Good_(50),
  maxJetEta_Good_(2.4)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f) {
         f = new TFile("Memory Directory");
         f->cd("Rint:/");
      }
      tree = (TTree*)gDirectory->Get("makeTopologyNtuple/tree");

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("makeTopologyNtuple/tree","");
      chain->Add("data/Ntuple_LM13-7TeV_1.root/makeTopologyNtuple/tree");
      chain->Add("data/Ntuple_LM13-7TeV_2.root/makeTopologyNtuple/tree");
      tree = chain;
#endif // SINGLE_TREE

   }

   Init(tree);
}

eventLoop::~eventLoop()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t eventLoop::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t eventLoop::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void eventLoop::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("numEle", &numEle, &b_numEle);
   fChain->SetBranchAddress("eleE", eleE, &b_eleE);
   fChain->SetBranchAddress("eleET", eleET, &b_eleET);
   fChain->SetBranchAddress("elePX", elePX, &b_elePX);
   fChain->SetBranchAddress("elePY", elePY, &b_elePY);
   fChain->SetBranchAddress("elePZ", elePZ, &b_elePZ);
   fChain->SetBranchAddress("elePhi", elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleTheta", eleTheta, &b_eleTheta);
   fChain->SetBranchAddress("eleEta", eleEta, &b_eleEta);
   fChain->SetBranchAddress("eleCharge", eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleIDQuality", eleIDQuality, &b_eleIDQuality);
   fChain->SetBranchAddress("eleTrackPt", eleTrackPt, &b_eleTrackPt);
   fChain->SetBranchAddress("eleTrackPhi", eleTrackPhi, &b_eleTrackPhi);
   fChain->SetBranchAddress("eleTrackEta", eleTrackEta, &b_eleTrackEta);
   fChain->SetBranchAddress("eleTrackChi2", eleTrackChi2, &b_eleTrackChi2);
   fChain->SetBranchAddress("eleTrackNDOF", eleTrackNDOF, &b_eleTrackNDOF);
   fChain->SetBranchAddress("eleTrackD0", eleTrackD0, &b_eleTrackD0);
   fChain->SetBranchAddress("eleBeamSpotCorrectedTrackD0", eleBeamSpotCorrectedTrackD0, &b_eleBeamSpotCorrectedTrackD0);
   fChain->SetBranchAddress("eleTrackDz", eleTrackDz, &b_eleTrackDz);
   fChain->SetBranchAddress("eleSCEta", eleSCEta, &b_eleSCEta);
   fChain->SetBranchAddress("eleSCE", eleSCE, &b_eleSCE);
   fChain->SetBranchAddress("eleSCPhi", eleSCPhi, &b_eleSCPhi);
   fChain->SetBranchAddress("eleSCSigmaEtaEta", eleSCSigmaEtaEta, &b_eleSCSigmaEtaEta);
   fChain->SetBranchAddress("eleSCSigmaIEtaIEta", eleSCSigmaIEtaIEta, &b_eleSCSigmaIEtaIEta);
   fChain->SetBranchAddress("eleSCE1x5", eleSCE1x5, &b_eleSCE1x5);
   fChain->SetBranchAddress("eleSCE5x5", eleSCE5x5, &b_eleSCE5x5);
   fChain->SetBranchAddress("eleSCE2x5max", eleSCE2x5max, &b_eleSCE2x5max);
   fChain->SetBranchAddress("eleTrackIso", eleTrackIso, &b_eleTrackIso);
   fChain->SetBranchAddress("eleEcalIso", eleEcalIso, &b_eleEcalIso);
   fChain->SetBranchAddress("eledr04EcalRecHitSumEt", eledr04EcalRecHitSumEt, &b_eledr04EcalRecHitSumEt);
   fChain->SetBranchAddress("eledr03EcalRecHitSumEt", eledr03EcalRecHitSumEt, &b_eledr03EcalRecHitSumEt);
   fChain->SetBranchAddress("eleHcalIso", eleHcalIso, &b_eleHcalIso);
   fChain->SetBranchAddress("eleEcalIsoDeposit", eleEcalIsoDeposit, &b_eleEcalIsoDeposit);
   fChain->SetBranchAddress("eleHcalIsoDeposit", eleHcalIsoDeposit, &b_eleHcalIsoDeposit);
   fChain->SetBranchAddress("eleComRelIso", eleComRelIso, &b_eleComRelIso);
   fChain->SetBranchAddress("elePhotonConversionTag", elePhotonConversionTag, &b_elePhotonConversionTag);
   fChain->SetBranchAddress("elePhotonConversionDist", elePhotonConversionDist, &b_elePhotonConversionDist);
   fChain->SetBranchAddress("elePhotonConversionDcot", elePhotonConversionDcot, &b_elePhotonConversionDcot);
   fChain->SetBranchAddress("eleTriggerMatch", eleTriggerMatch, &b_eleTriggerMatch);
   fChain->SetBranchAddress("eleJetOverlap", eleJetOverlap, &b_eleJetOverlap);
   fChain->SetBranchAddress("genEleET", genEleET, &b_genEleET);
   fChain->SetBranchAddress("genElePX", genElePX, &b_genElePX);
   fChain->SetBranchAddress("genElePY", genElePY, &b_genElePY);
   fChain->SetBranchAddress("genElePZ", genElePZ, &b_genElePZ);
   fChain->SetBranchAddress("genElePhi", genElePhi, &b_genElePhi);
   fChain->SetBranchAddress("genEleTheta", genEleTheta, &b_genEleTheta);
   fChain->SetBranchAddress("genEleEta", genEleEta, &b_genEleEta);
   fChain->SetBranchAddress("genEleCharge", genEleCharge, &b_genEleCharge);
   fChain->SetBranchAddress("numMuo", &numMuo, &b_numMuo);
   fChain->SetBranchAddress("muoE", muoE, &b_muoE);
   fChain->SetBranchAddress("muoET", muoET, &b_muoET);
   fChain->SetBranchAddress("muoPt", muoPt, &b_muoPt);
   fChain->SetBranchAddress("muoPX", muoPX, &b_muoPX);
   fChain->SetBranchAddress("muoPY", muoPY, &b_muoPY);
   fChain->SetBranchAddress("muoPZ", muoPZ, &b_muoPZ);
   fChain->SetBranchAddress("muoPhi", muoPhi, &b_muoPhi);
   fChain->SetBranchAddress("muoTheta", muoTheta, &b_muoTheta);
   fChain->SetBranchAddress("muoEta", muoEta, &b_muoEta);
   fChain->SetBranchAddress("muoCharge", muoCharge, &b_muoCharge);
   fChain->SetBranchAddress("muonChi2", muonChi2, &b_muonChi2);
   fChain->SetBranchAddress("muonD0", muonD0, &b_muonD0);
   fChain->SetBranchAddress("muonBeamSpotCorrectedD0", muonBeamSpotCorrectedD0, &b_muonBeamSpotCorrectedD0);
   fChain->SetBranchAddress("muonTrackNHits", muonTrackNHits, &b_muonTrackNHits);
   fChain->SetBranchAddress("muonNDOF", muonNDOF, &b_muonNDOF);
   fChain->SetBranchAddress("muonTrackIso", muonTrackIso, &b_muonTrackIso);
   fChain->SetBranchAddress("muonEcalIso", muonEcalIso, &b_muonEcalIso);
   fChain->SetBranchAddress("muonHcalIso", muonHcalIso, &b_muonHcalIso);
   fChain->SetBranchAddress("muonComRelIso", muonComRelIso, &b_muonComRelIso);
   fChain->SetBranchAddress("genMuoET", genMuoET, &b_genMuoET);
   fChain->SetBranchAddress("genMuoPX", genMuoPX, &b_genMuoPX);
   fChain->SetBranchAddress("genMuoPY", genMuoPY, &b_genMuoPY);
   fChain->SetBranchAddress("genMuoPZ", genMuoPZ, &b_genMuoPZ);
   fChain->SetBranchAddress("genMuoPhi", genMuoPhi, &b_genMuoPhi);
   fChain->SetBranchAddress("genMuoTheta", genMuoTheta, &b_genMuoTheta);
   fChain->SetBranchAddress("genMuoEta", genMuoEta, &b_genMuoEta);
   fChain->SetBranchAddress("genMuoCharge", genMuoCharge, &b_genMuoCharge);
   fChain->SetBranchAddress("nT", &nT, &b_nT);
   fChain->SetBranchAddress("nThadronic", &nThadronic, &b_nThadronic);
   fChain->SetBranchAddress("T_hadronicMCTruthE", T_hadronicMCTruthE, &b_T_hadronicMCTruthE);
   fChain->SetBranchAddress("T_hadronicMCTruthEt", T_hadronicMCTruthEt, &b_T_hadronicMCTruthEt);
   fChain->SetBranchAddress("T_hadronicMCTruthPx", T_hadronicMCTruthPx, &b_T_hadronicMCTruthPx);
   fChain->SetBranchAddress("T_hadronicMCTruthPy", T_hadronicMCTruthPy, &b_T_hadronicMCTruthPy);
   fChain->SetBranchAddress("T_hadronicMCTruthPz", T_hadronicMCTruthPz, &b_T_hadronicMCTruthPz);
   fChain->SetBranchAddress("T_hadronicMCMotherIndex", T_hadronicMCMotherIndex, &b_T_hadronicMCMotherIndex);
   fChain->SetBranchAddress("nTleptonic", &nTleptonic, &b_nTleptonic);
   fChain->SetBranchAddress("T_leptonicMCTruthE", T_leptonicMCTruthE, &b_T_leptonicMCTruthE);
   fChain->SetBranchAddress("T_leptonicMCTruthEt", T_leptonicMCTruthEt, &b_T_leptonicMCTruthEt);
   fChain->SetBranchAddress("T_leptonicMCTruthPx", T_leptonicMCTruthPx, &b_T_leptonicMCTruthPx);
   fChain->SetBranchAddress("T_leptonicMCTruthPy", T_leptonicMCTruthPy, &b_T_leptonicMCTruthPy);
   fChain->SetBranchAddress("T_leptonicMCTruthPz", T_leptonicMCTruthPz, &b_T_leptonicMCTruthPz);
   fChain->SetBranchAddress("T_leptonicMCMotherIndex", T_leptonicMCMotherIndex, &b_T_leptonicMCMotherIndex);
   fChain->SetBranchAddress("nb", &nb, &b_nb);
   fChain->SetBranchAddress("bMCTruthE", bMCTruthE, &b_bMCTruthE);
   fChain->SetBranchAddress("bMCTruthEt", bMCTruthEt, &b_bMCTruthEt);
   fChain->SetBranchAddress("bMCTruthPx", bMCTruthPx, &b_bMCTruthPx);
   fChain->SetBranchAddress("bMCTruthPy", bMCTruthPy, &b_bMCTruthPy);
   fChain->SetBranchAddress("bMCTruthPz", bMCTruthPz, &b_bMCTruthPz);
   fChain->SetBranchAddress("bMCTruthMother", bMCTruthMother, &b_bMCTruthMother);
   fChain->SetBranchAddress("nWhadronic", &nWhadronic, &b_nWhadronic);
   fChain->SetBranchAddress("W_hadronicMCTruthE", W_hadronicMCTruthE, &b_W_hadronicMCTruthE);
   fChain->SetBranchAddress("W_hadronicMCTruthEt", W_hadronicMCTruthEt, &b_W_hadronicMCTruthEt);
   fChain->SetBranchAddress("W_hadronicMCTruthPx", W_hadronicMCTruthPx, &b_W_hadronicMCTruthPx);
   fChain->SetBranchAddress("W_hadronicMCTruthPy", W_hadronicMCTruthPy, &b_W_hadronicMCTruthPy);
   fChain->SetBranchAddress("W_hadronicMCTruthPz", W_hadronicMCTruthPz, &b_W_hadronicMCTruthPz);
   fChain->SetBranchAddress("W_hadronicMCTruthPID", W_hadronicMCTruthPID, &b_W_hadronicMCTruthPID);
   fChain->SetBranchAddress("W_hadronicMCTruthMother", W_hadronicMCTruthMother, &b_W_hadronicMCTruthMother);
   fChain->SetBranchAddress("nWleptonic", &nWleptonic, &b_nWleptonic);
   fChain->SetBranchAddress("W_leptonicMCTruthE", W_leptonicMCTruthE, &b_W_leptonicMCTruthE);
   fChain->SetBranchAddress("W_leptonicMCTruthEt", W_leptonicMCTruthEt, &b_W_leptonicMCTruthEt);
   fChain->SetBranchAddress("W_leptonicMCTruthPx", W_leptonicMCTruthPx, &b_W_leptonicMCTruthPx);
   fChain->SetBranchAddress("W_leptonicMCTruthPy", W_leptonicMCTruthPy, &b_W_leptonicMCTruthPy);
   fChain->SetBranchAddress("W_leptonicMCTruthPz", W_leptonicMCTruthPz, &b_W_leptonicMCTruthPz);
   fChain->SetBranchAddress("W_leptonicMCTruthPID", W_leptonicMCTruthPID, &b_W_leptonicMCTruthPID);
   fChain->SetBranchAddress("W_leptonicMCTruthMother", W_leptonicMCTruthMother, &b_W_leptonicMCTruthMother);
   fChain->SetBranchAddress("isElePlusJets", &isElePlusJets, &b_isElePlusJets);
   fChain->SetBranchAddress("VQQBosonAbsId", &VQQBosonAbsId, &b_VQQBosonAbsId);
   fChain->SetBranchAddress("numJet", &numJet, &b_numJet);
   fChain->SetBranchAddress("jetE", jetE, &b_jetE);
   fChain->SetBranchAddress("jetEt", jetEt, &b_jetEt);
   fChain->SetBranchAddress("jetPt", jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetCorEt", jetCorEt, &b_jetCorEt);
   fChain->SetBranchAddress("jetEta", jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetTheta", jetTheta, &b_jetTheta);
   fChain->SetBranchAddress("jetPhi", jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetPx", jetPx, &b_jetPx);
   fChain->SetBranchAddress("jetPy", jetPy, &b_jetPy);
   fChain->SetBranchAddress("jetPz", jetPz, &b_jetPz);
   fChain->SetBranchAddress("jetNtracksInJet", jetNtracksInJet, &b_jetNtracksInJet);
   fChain->SetBranchAddress("jetJetCharge", jetJetCharge, &b_jetJetCharge);
   fChain->SetBranchAddress("jetMuEnergy", jetMuEnergy, &b_jetMuEnergy);
   fChain->SetBranchAddress("jetMuEnergyFraction", jetMuEnergyFraction, &b_jetMuEnergyFraction);
   fChain->SetBranchAddress("jetChargedMultiplicity", jetChargedMultiplicity, &b_jetChargedMultiplicity);
   fChain->SetBranchAddress("jetNeutralHadEnergy", jetNeutralHadEnergy, &b_jetNeutralHadEnergy);
   fChain->SetBranchAddress("jetEMEnergyInEB", jetEMEnergyInEB, &b_jetEMEnergyInEB);
   fChain->SetBranchAddress("jetEMEnergyInEE", jetEMEnergyInEE, &b_jetEMEnergyInEE);
   fChain->SetBranchAddress("jetEMEnergyFraction", jetEMEnergyFraction, &b_jetEMEnergyFraction);
   fChain->SetBranchAddress("jetEMEnergyInHF", jetEMEnergyInHF, &b_jetEMEnergyInHF);
   fChain->SetBranchAddress("jetHadEnergyInHB", jetHadEnergyInHB, &b_jetHadEnergyInHB);
   fChain->SetBranchAddress("jetHadEnergyInHE", jetHadEnergyInHE, &b_jetHadEnergyInHE);
   fChain->SetBranchAddress("jetHadEnergyInHF", jetHadEnergyInHF, &b_jetHadEnergyInHF);
   fChain->SetBranchAddress("jetHadEnergyInHO", jetHadEnergyInHO, &b_jetHadEnergyInHO);
   fChain->SetBranchAddress("jetBtagTrackCountHighPurity", jetBtagTrackCountHighPurity, &b_jetBtagTrackCountHighPurity);
   fChain->SetBranchAddress("jetBtagTrackCountHighEff", jetBtagTrackCountHighEff, &b_jetBtagTrackCountHighEff);
   fChain->SetBranchAddress("jetBtagProbability", jetBtagProbability, &b_jetBtagProbability);
   fChain->SetBranchAddress("jetBtagSoftElectron", jetBtagSoftElectron, &b_jetBtagSoftElectron);
   fChain->SetBranchAddress("jetBtagSoftMuon", jetBtagSoftMuon, &b_jetBtagSoftMuon);
   fChain->SetBranchAddress("jetBtagSoftMuonNoIP", jetBtagSoftMuonNoIP, &b_jetBtagSoftMuonNoIP);
   fChain->SetBranchAddress("jetBtagSoftMuonPtRel", jetBtagSoftMuonPtRel, &b_jetBtagSoftMuonPtRel);
   fChain->SetBranchAddress("jetBtagSoftMuonQuality", jetBtagSoftMuonQuality, &b_jetBtagSoftMuonQuality);
   fChain->SetBranchAddress("jetBtagSecondaryVertex", jetBtagSecondaryVertex, &b_jetBtagSecondaryVertex);
   fChain->SetBranchAddress("jetBtagSecondaryVertexNegative", jetBtagSecondaryVertexNegative, &b_jetBtagSecondaryVertexNegative);
   fChain->SetBranchAddress("jetBtagCombinedSVLL", jetBtagCombinedSVLL, &b_jetBtagCombinedSVLL);
   fChain->SetBranchAddress("jetBtagCombinedSVMVA", jetBtagCombinedSVMVA, &b_jetBtagCombinedSVMVA);
   fChain->SetBranchAddress("jetCorrFactor", jetCorrFactor, &b_jetCorrFactor);
   fChain->SetBranchAddress("jetN60", jetN60, &b_jetN60);
   fChain->SetBranchAddress("jetN90", jetN90, &b_jetN90);
   fChain->SetBranchAddress("jetNeutralEmEnergy", jetNeutralEmEnergy, &b_jetNeutralEmEnergy);
   fChain->SetBranchAddress("jetTriggered", jetTriggered, &b_jetTriggered);
   fChain->SetBranchAddress("jetSVPT", jetSVPT, &b_jetSVPT);
   fChain->SetBranchAddress("jetSVL2D", jetSVL2D, &b_jetSVL2D);
   fChain->SetBranchAddress("jetSVL2Dxy", jetSVL2Dxy, &b_jetSVL2Dxy);
   fChain->SetBranchAddress("jetSVL2DxyErr", jetSVL2DxyErr, &b_jetSVL2DxyErr);
   fChain->SetBranchAddress("jetSVL2DxySig", jetSVL2DxySig, &b_jetSVL2DxySig);
   fChain->SetBranchAddress("jetSVL3D", jetSVL3D, &b_jetSVL3D);
   fChain->SetBranchAddress("jetSVL3DErr", jetSVL3DErr, &b_jetSVL3DErr);
   fChain->SetBranchAddress("jetSVL3DSig", jetSVL3DSig, &b_jetSVL3DSig);
   fChain->SetBranchAddress("jetSVMass", jetSVMass, &b_jetSVMass);
   fChain->SetBranchAddress("jetSVNtracks", jetSVNtracks, &b_jetSVNtracks);
   fChain->SetBranchAddress("genJetET", genJetET, &b_genJetET);
   fChain->SetBranchAddress("genJetPX", genJetPX, &b_genJetPX);
   fChain->SetBranchAddress("genJetPY", genJetPY, &b_genJetPY);
   fChain->SetBranchAddress("genJetPZ", genJetPZ, &b_genJetPZ);
   fChain->SetBranchAddress("genJetPhi", genJetPhi, &b_genJetPhi);
   fChain->SetBranchAddress("genJetTheta", genJetTheta, &b_genJetTheta);
   fChain->SetBranchAddress("genJetEta", genJetEta, &b_genJetEta);
   fChain->SetBranchAddress("genJetPID", genJetPID, &b_genJetPID);
   fChain->SetBranchAddress("jetPID", jetPID, &b_jetPID);
   fChain->SetBranchAddress("jetClosestBPartonDeltaR", jetClosestBPartonDeltaR, &b_jetClosestBPartonDeltaR);
   fChain->SetBranchAddress("jetClosestCPartonDeltaR", jetClosestCPartonDeltaR, &b_jetClosestCPartonDeltaR);
   fChain->SetBranchAddress("btagParamDiscCut_MISTAGSSVM", &btagParamDiscCut_MISTAGSSVM, &b_btagParamDiscCut_MISTAGSSVM);
   fChain->SetBranchAddress("btagParamDiscCut_PTRELSSVM", &btagParamDiscCut_PTRELSSVM, &b_btagParamDiscCut_PTRELSSVM);
   fChain->SetBranchAddress("btagParamDiscCut_PTRELTCHEL", &btagParamDiscCut_PTRELTCHEL, &b_btagParamDiscCut_PTRELTCHEL);
   fChain->SetBranchAddress("btagParamDiscCut_PTRELTCHEM", &btagParamDiscCut_PTRELTCHEM, &b_btagParamDiscCut_PTRELTCHEM);
   fChain->SetBranchAddress("btagParamDiscCut_PTRELTCHET", &btagParamDiscCut_PTRELTCHET, &b_btagParamDiscCut_PTRELTCHET);
   fChain->SetBranchAddress("btagParamDiscCut_PTRELTCHPL", &btagParamDiscCut_PTRELTCHPL, &b_btagParamDiscCut_PTRELTCHPL);
   fChain->SetBranchAddress("btagParamDiscCut_PTRELTCHPM", &btagParamDiscCut_PTRELTCHPM, &b_btagParamDiscCut_PTRELTCHPM);
   fChain->SetBranchAddress("btagParamDiscCut_PTRELTCHPT", &btagParamDiscCut_PTRELTCHPT, &b_btagParamDiscCut_PTRELTCHPT);
   fChain->SetBranchAddress("btagParamDiscCut_SYSTEM8SSVM", &btagParamDiscCut_SYSTEM8SSVM, &b_btagParamDiscCut_SYSTEM8SSVM);
   fChain->SetBranchAddress("jetBtagParam_MISTAGSSVM_BTAGLEFF", jetBtagParam_MISTAGSSVM_BTAGLEFF, &b_jetBtagParam_MISTAGSSVM_BTAGLEFF);
   fChain->SetBranchAddress("jetBtagParam_MISTAGSSVM_BTAGLERR", jetBtagParam_MISTAGSSVM_BTAGLERR, &b_jetBtagParam_MISTAGSSVM_BTAGLERR);
   fChain->SetBranchAddress("jetBtagParam_MISTAGSSVM_BTAGLEFFCORR", jetBtagParam_MISTAGSSVM_BTAGLEFFCORR, &b_jetBtagParam_MISTAGSSVM_BTAGLEFFCORR);
   fChain->SetBranchAddress("jetBtagParam_MISTAGSSVM_BTAGLERRCORR", jetBtagParam_MISTAGSSVM_BTAGLERRCORR, &b_jetBtagParam_MISTAGSSVM_BTAGLERRCORR);
   fChain->SetBranchAddress("jetBtagParam_PTRELSSVM_BTAGBEFFCORR", jetBtagParam_PTRELSSVM_BTAGBEFFCORR, &b_jetBtagParam_PTRELSSVM_BTAGBEFFCORR);
   fChain->SetBranchAddress("jetBtagParam_PTRELSSVM_BTAGBERRCORR", jetBtagParam_PTRELSSVM_BTAGBERRCORR, &b_jetBtagParam_PTRELSSVM_BTAGBERRCORR);
   fChain->SetBranchAddress("jetBtagParam_SYSTEM8SSVM_BTAGBEFF", jetBtagParam_SYSTEM8SSVM_BTAGBEFF, &b_jetBtagParam_SYSTEM8SSVM_BTAGBEFF);
   fChain->SetBranchAddress("jetBtagParam_SYSTEM8SSVM_BTAGBERR", jetBtagParam_SYSTEM8SSVM_BTAGBERR, &b_jetBtagParam_SYSTEM8SSVM_BTAGBERR);
   fChain->SetBranchAddress("jetBtagParam_PTRELTCHEL_BTAGBEFFCORR", jetBtagParam_PTRELTCHEL_BTAGBEFFCORR, &b_jetBtagParam_PTRELTCHEL_BTAGBEFFCORR);
   fChain->SetBranchAddress("jetBtagParam_PTRELTCHEM_BTAGBEFFCORR", jetBtagParam_PTRELTCHEM_BTAGBEFFCORR, &b_jetBtagParam_PTRELTCHEM_BTAGBEFFCORR);
   fChain->SetBranchAddress("jetBtagParam_PTRELTCHET_BTAGBEFFCORR", jetBtagParam_PTRELTCHET_BTAGBEFFCORR, &b_jetBtagParam_PTRELTCHET_BTAGBEFFCORR);
   fChain->SetBranchAddress("jetBtagParam_PTRELTCHPL_BTAGBEFFCORR", jetBtagParam_PTRELTCHPL_BTAGBEFFCORR, &b_jetBtagParam_PTRELTCHPL_BTAGBEFFCORR);
   fChain->SetBranchAddress("jetBtagParam_PTRELTCHPM_BTAGBEFFCORR", jetBtagParam_PTRELTCHPM_BTAGBEFFCORR, &b_jetBtagParam_PTRELTCHPM_BTAGBEFFCORR);
   fChain->SetBranchAddress("jetBtagParam_PTRELTCHPT_BTAGBEFFCORR", jetBtagParam_PTRELTCHPT_BTAGBEFFCORR, &b_jetBtagParam_PTRELTCHPT_BTAGBEFFCORR);
   fChain->SetBranchAddress("numGeneralTracks", &numGeneralTracks, &b_numGeneralTracks);
   fChain->SetBranchAddress("generalTracksPt", generalTracksPt, &b_generalTracksPt);
   fChain->SetBranchAddress("generalTracksEta", generalTracksEta, &b_generalTracksEta);
   fChain->SetBranchAddress("generalTracksTheta", generalTracksTheta, &b_generalTracksTheta);
   fChain->SetBranchAddress("generalTracksBeamSpotCorrectedD0", generalTracksBeamSpotCorrectedD0, &b_generalTracksBeamSpotCorrectedD0);
   fChain->SetBranchAddress("generalTracksPhi", generalTracksPhi, &b_generalTracksPhi);
   fChain->SetBranchAddress("generalTracksCharge", generalTracksCharge, &b_generalTracksCharge);
   fChain->SetBranchAddress("processId", &processId, &b_processId);
   fChain->SetBranchAddress("processPtHat", &processPtHat, &b_processPtHat);
   fChain->SetBranchAddress("processMCWeight", &processMCWeight, &b_processMCWeight);
   fChain->SetBranchAddress("beamSpotX", &beamSpotX, &b_beamSpotX);
   fChain->SetBranchAddress("beamSpotY", &beamSpotY, &b_beamSpotY);
   fChain->SetBranchAddress("beamSpotZ", &beamSpotZ, &b_beamSpotZ);
   fChain->SetBranchAddress("topo_sphericity", &topo_sphericity, &b_topo_sphericity);
   fChain->SetBranchAddress("topo_aplanarity", &topo_aplanarity, &b_topo_aplanarity);
   fChain->SetBranchAddress("topo_sphericity_e", &topo_sphericity_e, &b_topo_sphericity_e);
   fChain->SetBranchAddress("topo_aplanarity_e", &topo_aplanarity_e, &b_topo_aplanarity_e);
   fChain->SetBranchAddress("topo_ht", &topo_ht, &b_topo_ht);
   fChain->SetBranchAddress("topo_ht_e", &topo_ht_e, &b_top_ht_e);
   fChain->SetBranchAddress("topo_sqrts", &topo_sqrts, &b_topo_sqrts);
   fChain->SetBranchAddress("topo_sqrts_e", &topo_sqrts_e, &b_top_sqrts_e);
   fChain->SetBranchAddress("topo_oblateness", &topo_oblateness, &b_topo_oblateness);
   fChain->SetBranchAddress("metEt", &metEt, &b_metEt);
   fChain->SetBranchAddress("metPhi", &metPhi, &b_metPhi);
   fChain->SetBranchAddress("metPt", &metPt, &b_metPt);
   fChain->SetBranchAddress("metPx", &metPx, &b_metPx);
   fChain->SetBranchAddress("metPy", &metPy, &b_metPy);
   fChain->SetBranchAddress("genMetEt", &genMetEt, &b_genMetEt);
   fChain->SetBranchAddress("genMetPhi", &genMetPhi, &b_genMetPhi);
   fChain->SetBranchAddress("genMetPt", &genMetPt, &b_genMetPt);
   fChain->SetBranchAddress("genMetPx", &genMetPx, &b_genMetPx);
   fChain->SetBranchAddress("genMetPy", &genMetPy, &b_genMetPy);
   fChain->SetBranchAddress("metMaxEtEM", &metMaxEtEM, &b_metMaxEtEM);
   fChain->SetBranchAddress("metMaxEtHad", &metMaxEtHad, &b_metMaxEtHad);
   fChain->SetBranchAddress("metEtFracHad", &metEtFracHad, &b_metEtFracHad);
   fChain->SetBranchAddress("metEtFracEM", &metEtFracEM, &b_metEtFracEM);
   fChain->SetBranchAddress("metHadEtHB", &metHadEtHB, &b_metHadEtHB);
   fChain->SetBranchAddress("metHadEtHO", &metHadEtHO, &b_metHadEtHO);
   fChain->SetBranchAddress("metHadEtHF", &metHadEtHF, &b_metHadEtHF);
   fChain->SetBranchAddress("metHadEtHE", &metHadEtHE, &b_metHadEtHE);
   fChain->SetBranchAddress("metEmEtHF", &metEmEtHF, &b_metEmEtHF);
   fChain->SetBranchAddress("metEmEtEE", &metEmEtEE, &b_metEmEtEE);
   fChain->SetBranchAddress("metEmEtEB", &metEmEtEB, &b_metEmEtEB);
   fChain->SetBranchAddress("metSignificance", &metSignificance, &b_metSignificance);
   fChain->SetBranchAddress("metScalarEt", &metScalarEt, &b_metScalarEt);
   fChain->SetBranchAddress("metEtUncorrected", &metEtUncorrected, &b_metEtUncorrected);
   fChain->SetBranchAddress("metPhiUncorrected", &metPhiUncorrected, &b_metPhiUncorrected);
   fChain->SetBranchAddress("mhtPt", &mhtPt, &b_mhtPt);
   fChain->SetBranchAddress("mhtPy", &mhtPy, &b_mhtPy);
   fChain->SetBranchAddress("mhtPx", &mhtPx, &b_mhtPx);
   fChain->SetBranchAddress("mhtPhi", &mhtPhi, &b_mhtPhi);
   fChain->SetBranchAddress("mhtSumEt", &mhtSumEt, &b_mhtSumEt);
   fChain->SetBranchAddress("mhtSignif", &mhtSignif, &b_mhtSignif);
   fChain->SetBranchAddress("nZCandidates", &nZCandidates, &b_nZCandidates);
   fChain->SetBranchAddress("ZCandidates", ZCandidates, &b_ZCandidates);
   fChain->SetBranchAddress("nTriggerBits", &nTriggerBits, &b_nTriggerBits);
   fChain->SetBranchAddress("TriggerBits", TriggerBits, &b_TriggerBits);
   fChain->SetBranchAddress("HLT_Ele15_LW_L1R", &HLT_Ele15_LW_L1R, &b_HLT_Ele15_LW_L1R);
   fChain->SetBranchAddress("flavorhistory", &flavorhistory, &b_flavorhistory);
   fChain->SetBranchAddress("myProcess", &myProcess, &b_myProcess);
   fChain->SetBranchAddress("genParEta", genParEta, &b_genParEta);
   fChain->SetBranchAddress("genParPhi", genParPhi, &b_genParPhi);
   fChain->SetBranchAddress("genParTheta", genParTheta, &b_genParTheta);
   fChain->SetBranchAddress("genParE", genParE, &b_genParE);
   fChain->SetBranchAddress("genParEt", genParEt, &b_genParEt);
   fChain->SetBranchAddress("genParP", genParP, &b_genParP);
   fChain->SetBranchAddress("genParPt", genParPt, &b_genParPt);
   fChain->SetBranchAddress("genParId", genParId, &b_genParId);
   fChain->SetBranchAddress("genParStat", genParStat, &b_genParStat);
   fChain->SetBranchAddress("ndaughters", ndaughters, &b_ndaughters);
   fChain->SetBranchAddress("genParCharge", genParCharge, &b_genParCharge);
   fChain->SetBranchAddress("genEta", genEta, &b_genEta);
   fChain->SetBranchAddress("genPhi", genPhi, &b_genPhi);
   fChain->SetBranchAddress("genTheta", genTheta, &b_genTheta);
   fChain->SetBranchAddress("genE", genE, &b_genE);
   fChain->SetBranchAddress("genEt", genEt, &b_genEt);
   fChain->SetBranchAddress("genP", genP, &b_genP);
   fChain->SetBranchAddress("genPt", genPt, &b_genPt);
   fChain->SetBranchAddress("genId", genId, &b_genId);
   fChain->SetBranchAddress("genStat", genStat, &b_genStat);
   fChain->SetBranchAddress("genCharge", genCharge, &b_genCharge);
   fChain->SetBranchAddress("numTau", &numTau, &b_numTau);
   fChain->SetBranchAddress("tauE", tauE, &b_tauE);
   fChain->SetBranchAddress("tauPt", tauPt, &b_tauPt);
   fChain->SetBranchAddress("tauPhi", tauPhi, &b_tauPhi);
   fChain->SetBranchAddress("tauEta", tauEta, &b_tauEta);
   fChain->SetBranchAddress("numPhoton", &numPhoton, &b_numPhoton);
   fChain->SetBranchAddress("photonE", photonE, &b_photonE);
   fChain->SetBranchAddress("photonPt", photonPt, &b_photonPt);
   fChain->SetBranchAddress("photonPhi", photonPhi, &b_photonPhi);
   fChain->SetBranchAddress("photonEta", photonEta, &b_photonEta);
   fChain->SetBranchAddress("eventRun", &eventRun, &b_eventRun);
   fChain->SetBranchAddress("eventNum", &eventNum, &b_eventNum);
   fChain->SetBranchAddress("eventLumiblock", &eventLumiblock, &b_eventLumiblock);
   Notify();
}

Bool_t eventLoop::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void eventLoop::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t eventLoop::Cut()
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.

  if (nGoodJets() < minGoodJets_) return -1;
  //number of primary vertices is missing....
  //i'm going to remove all jet pt cuts for the moment
  //  if (jetPt[0] <=100) return -1; //jets are already sorted by Pt
  
  if ( nGoodElectrons() != 0 ) return -1;
  if ( nGoodMuons() != 0 ) return -1;

  return 1;
}

UInt_t eventLoop::nGoodJets() {
  UInt_t ngood=0;
  for (int i=0; i< numJet; i++) {
    if ( isJetGood(i) ) ngood++;
  }
 
  return ngood;
}

Bool_t eventLoop::isJetGood(UInt_t index) {
  return  ( jetPt[index] > minJetPt_Good_ && fabs(jetEta[index])< maxJetEta_Good_ );
}

UInt_t eventLoop::nGoodElectrons() {

  UInt_t ngood = 0;
  for ( int i=0; i< numEle; i++) {
    
    if (eleTrackPt[i] < 15 ) continue;
    if (fabs(eleEta[i]) > 2.5) continue;

    //FIXME...i'm taking this directly from my old code
    //but why would it not be fabs() of this?
    if ( eleBeamSpotCorrectedTrackD0[i] > 0.2 ) continue;
    //&& (fabs(ele->eta)>1.567 || fabs(ele->eta)<1.47) );

    //this may not be exactly the same iso as i was using before
    if ( eleComRelIso[i] >0.5 ) continue;

    //from when i was using dmp's code. how to duplicate?
    //    if ( ele->IDLoose != 1) continue;

    //if we get to here it must be good
    ngood++;
  }
  return ngood;
}

UInt_t eventLoop::nGoodMuons() {
  UInt_t ngood=0;
  for (int i=0; i<numMuo; i++) {
    //it doesn't look like we have the muon ID
    if (muoPt[i] < 10) continue;
    if ( fabs( muoEta[i]) >2.4) continue;
    //probably not the same variable that i was using before
    if (muonComRelIso[i] >0.1) continue;

    if (muonNDOF[i] == 0) {
      std::cout<<"muon ndof is 0!"<<std::endl;
      continue;
    }
    //i am assuming this is chi^2 and not chi^2/dof
    if ( muonChi2[i]/muonNDOF[i] > 10) continue;
    //FIXME ... again, shouldn't this be fabs()?
    if (muonBeamSpotCorrectedD0[i] > 0.2) continue;
    if (muonTrackNHits[i] < 11) continue;

    ngood++;
  }
  return ngood;
}


double eventLoop::getDeltaPhi(Double_t phi1, Double_t phi2) {
  return acos(cos(phi1-phi2));
}


#endif // #ifdef eventLoop_cxx
