import FWCore.ParameterSet.Config as cms

isMC = False

process = cms.Process("BasicTreeMaker")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'INPUT'
    #'file:/cu1/joshmt/DonPAT/MoreMSSM_PAT_10_1_osg.root'
    #'file:/cu1/joshmt/DonPAT/TTbarJets_SPRING10_PAT_99_7.root'
#    '/store/user/wteo/MinimumBias/commissioning_jetmettau_jun14thskim/c4203bae27064fcbcff856a81ef504ae/PAT_362_commissioningJune14th_jetmettau_49_3_MI6.root'
#    'file:/cu1/joshmt/data/Run2010B/9E73F38F-ECC9-DF11-8BFD-00304879FC6C.root'
#                                'file:/cu1/joshmt/RECO/FAB95B26-0769-DF11-BF53-003048CB87DE.root'
    #    '/store/user/wteo/LM3/lm3_spring10_v1/838332be00fca3dd34042237dc5fc2e8/LM3_SPRING10_PAT_1_4.root'
    )
)
############################# START ntuple specifics ####################################

#b-tagging efficiency parameters:
process.load ("RecoBTag.PerformanceDB.PoolBTagPerformanceDBOctEx")
process.load ("RecoBTag.PerformanceDB.BTagPerformanceDBOctEx")
process.load ("RecoBTag.PerformanceDB.PoolBTagPerformanceDBMC36X")
process.load ("RecoBTag.PerformanceDB.BTagPerformanceDBMC36X")

from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector as pvSel
from PhysicsTools.SelectorUtils.jetIDSelector_cfi import jetIDSelector

process.BasicTreeMaker = cms.EDAnalyzer('BasicTreeMaker',
                                        MCflag = cms.bool(isMC),
                                        #b taggers to store in ntuple
                                        btagAlgorithms = cms.vstring("trackCountingHighPurBJetTags","trackCountingHighEffBJetTags",
                                                                     "simpleSecondaryVertexHighEffBJetTags",#"simpleSecondaryVertexNegativeBJetTags",
                                                                     "simpleSecondaryVertexHighPurBJetTags",#),
                                                                     "simpleSecondaryVertexBJetTags"), #for older samples
#extract individual tags from don's python
#first the names of various collections
                                        jetTag = cms.InputTag('selectedPatJets'), #cleanPatJetsAK5Calo
                                        caloMetTag = cms.InputTag('patMETs'), #patMETsAK5Calo 
                                        tcMetTag = cms.InputTag(''), #patMETsTC
                                        eleTag = cms.InputTag('selectedPatElectrons'), #cleanPatElectrons
                                        muoTag = cms.InputTag('selectedPatMuons'), #cleanPatMuons
#configuration of the cut flow and what cuts to use when storing to the ntuple
#the 'SUSY trigger' used in the cut flow is the first one listed here
                                        triggersOfInterest = cms.vstring("HLT_HT100U","HLT_HT120U","HLT_HT140U","HLT_HT200",
                                                                         "HLT_MET45","HLT_MET60","HLT_MET65","HLT_MET100",
                                                                         "HLT_Jet15U","HLT_Jet30U","HLT_Jet50U","HLT_Jet70U","HLT_Jet100U",
                                                                         "HLT_DiJetAve15U","HLT_DiJetAve30U","HLT_DiJetAve50U","HLT_DiJetAve70U",
                                                                         "HLT_BTagIP_Jet50U"),
            
                                        pvSelector = cms.PSet( pvSel.clone() ),
                                        
                                        jetIdLoose = jetIDSelector.clone(),

                                        ## muons
                                        muonId = cms.PSet(
                                            version = cms.string('FIRSTDATA'),
                                            Chi2 = cms.double(10.0),
                                            D0 = cms.double(0.20),
                                            ED0 = cms.double(999.0),
                                            SD0 = cms.double(3.0),
                                            NHits = cms.int32(11),
                                            NValMuHits = cms.int32(0),
                                            ECalVeto = cms.double(4.0),
                                            HCalVeto = cms.double(6.0),
                                            RelIso = cms.double(0.1),
                                            cutsToIgnore = cms.vstring('ED0', 'SD0', 'ECalVeto', 'HCalVeto'),
                                            RecalcFromBeamSpot = cms.bool(False),
                                            beamLineSrc = cms.InputTag("offlineBeamSpot")
                                            ),

                                        ## electrons
                                        electronId = cms.PSet(
                                            version = cms.string('FIRSTDATA'),
                                            D0 = cms.double(0.20),
                                            ED0 = cms.double(999.0),
                                            SD0 = cms.double(3.0),
                                            RelIso = cms.double( 0.1 ), #SUSYb selection
                                            #RelIso = cms.double( 0.5 ), #RA2 Reference Selection
                                            cutsToIgnore = cms.vstring('ED0', 'SD0')
                                            ),

                                        minJets = cms.int32( 3 ),
                                        
                                        jetPtMin = cms.double( 50.0 ),
                                        jetEtaMax = cms.double( 2.5 ),  #RA2 Reference Selection
                                        loosejetPtMin = cms.double( 30.0 ),
                                        loosejetEtaMax = cms.double( 5.0 ),
                                        muPtMin = cms.double( 10.0 ),
                                        muEtaMax = cms.double( 2.4 ),
                                        eleEtMin = cms.double( 15.0 ),
                                        eleEtaMax = cms.double( 2.4 ),

                                        mhtMin = cms.double( 150.0 ), #RA2 Reference Selection
                                        metMin = cms.double( 150.0 )
)


process.TFileService = cms.Service("TFileService", fileName = cms.string('BasicNtuple.root') )

process.p = cms.Path(process.BasicTreeMaker)
