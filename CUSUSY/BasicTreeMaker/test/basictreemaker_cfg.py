import FWCore.ParameterSet.Config as cms

#load configuration for Don's code
from CUSUSY.Selection.susybjetsSelection_cfi import *

process = cms.Process("BasicTreeMaker")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'INPUTFILE'
#    '/store/user/wteo/LM3/lm3_spring10_v1/838332be00fca3dd34042237dc5fc2e8/LM3_SPRING10_PAT_1_4.root'
    )
)

process.BasicTreeMaker = cms.EDAnalyzer('BasicTreeMaker',
                                        #b taggers to store in ntuple
                                        btagAlgorithms = cms.vstring("trackCountingHighPurBJetTags","trackCountingHighEffBJetTags",
                                                                     "simpleSecondaryVertexHighEffBJetTags",#"simpleSecondaryVertexNegativeBJetTags",
                                                                     "simpleSecondaryVertexHighPurBJetTags",
                                                                     "simpleSecondaryVertexBJetTags"), #for older samples
                                        susyBJetsSelection = susybjetsSelection,
#extract individual tags from don's python
#first the names of various collections
                                        jetTag = susybjetsSelection.jetSrc,
                                        caloMetTag = susybjetsSelection.metSrc, #we assume don's py is asking for calomet
                                        tcMetTag = cms.InputTag('patMETsTC'),
                                        eleTag = susybjetsSelection.electronSrc,
                                        muoTag = susybjetsSelection.muonSrc,
#configuration of the cut flow and what cuts to use when storing to the ntuple
#the 'SUSY trigger' used in the cut flow is the first one listed here
                                        triggersOfInterest = cms.vstring("HLT_HT100U","HLT_HT200",
                                                                         "HLT_MET45","HLT_MET60","HLT_MET100","HLT_Jet15U",
                                                                         "HLT_BTagIP_Jet50U"),
            
#                                        susyTrigger = susybjetsSelection.susyTrig,

                                        pvSelector = susybjetsSelection.pvSelector,
                                        
                                        jetIdLoose = susybjetsSelection.jetIdLoose,

                                        muonId = susybjetsSelection.muonId,
                                        electronId = susybjetsSelection.electronId,

                                        minJets = susybjetsSelection.minJets,
                                        
                                        jetPtMin = susybjetsSelection.jetPtMin,
                                        jetEtaMax = susybjetsSelection.jetEtaMax,
                                        loosejetPtMin = susybjetsSelection.loosejetPtMin,
                                        loosejetEtaMax = susybjetsSelection.loosejetEtaMax,
                                        muPtMin = susybjetsSelection.muPtMin,
                                        muEtaMax = susybjetsSelection.muEtaMax,
                                        eleEtMin = susybjetsSelection.eleEtMin,
                                        eleEtaMax = susybjetsSelection.eleEtaMax,

                                        mhtMin = susybjetsSelection.mhtMin,
                                        metMin = susybjetsSelection.metMin
)

process.TFileService = cms.Service("TFileService", fileName = cms.string('BasicNtuple.root') )

process.p = cms.Path(process.BasicTreeMaker)
