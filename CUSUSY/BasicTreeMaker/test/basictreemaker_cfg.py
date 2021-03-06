import FWCore.ParameterSet.Config as cms

isMC = True

process = cms.Process("BasicTreeMaker")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#   'file:/cu1/joshmt/AOD/Fall10-QCD_TuneD6T_HT-1000_7TeV-madgraph-C473857E-99DA-DF11-8531-00163691DC86.SUSYPAT.root'
#    'file:/cu1/joshmt/AOD/387/QCD_Pt_80to120_TuneZ2_pythia6_PU_Fall10_387_PAT_70_1_RCp.root'
#    'file:/cu1/joshmt/Validation/mSugra_tanBeta50_SUSYPAT.root'
  'INPUT'
 #   'file:/afs/cern.ch/user/s/ssekmen/public/Sezen_PAT.root'
#    'file:/cu1/joshmt/AOD/387/TTJets_TuneD6T_Fall10_387_PAT_9_1_rY1.root'
#   'file:/cu1/joshmt/AOD/387/PAT_387_run2010B_multijet_nov4rereco_9_1_9wD.root'
#    'file:/cu1/joshmt/DonPAT/data/PAT_38X_jetmettau_sep17rereco_sep11_9_1_wJJ.root'
#    'file:/cu1/joshmt/DonPAT/MoreMSSM_PAT_10_1_osg.root'
#    'file:/cu1/joshmt/DonPAT/WJets-SUSYPAT-FEEFD640-9277-DF11-948C-001731EB1E20.root'
#    'file:/cu1/joshmt/DonPAT/TTbarJets_SPRING10_PAT_99_7.root'
#    '/store/user/wteo/MinimumBias/commissioning_jetmettau_jun14thskim/c4203bae27064fcbcff856a81ef504ae/PAT_362_commissioningJune14th_jetmettau_49_3_MI6.root'
#    'file:/cu1/joshmt/data/Run2010B/9E73F38F-ECC9-DF11-8BFD-00304879FC6C.root'
    #    '/store/user/wteo/LM3/lm3_spring10_v1/838332be00fca3dd34042237dc5fc2e8/LM3_SPRING10_PAT_1_4.root'
    )
)

pvString = 'offlinePrimaryVertices'
eleCollection = 'cleanPatElectrons'
muonCollection = 'cleanPatMuons'

eleCollectionPF = 'selectedPatElectronsPF'
muonCollectionPF = 'selectedPatMuonsPF'

############################# START RA2 lepton specifics ####################################
from SandBox.Skims.electronSelector_cfi import *

process.RA2electronSelector = electronSelector.clone()
process.RA2electronSelector.ElectronSource = cms.InputTag(eleCollection)
process.RA2electronSelector.VertexSource = cms.InputTag(pvString)

process.RA2electronSelectorPF = electronSelector.clone()
process.RA2electronSelectorPF.ElectronSource = cms.InputTag(eleCollectionPF)
process.RA2electronSelectorPF.VertexSource = cms.InputTag(pvString)

from SandBox.Skims.muonSelector_cfi import *

process.RA2muonSelector = muonSelector.clone()
process.RA2muonSelector.MuonSource = cms.InputTag(muonCollection)
process.RA2muonSelector.VertexSource = cms.InputTag(pvString)

process.RA2muonSelectorPF = muonSelector.clone()
process.RA2muonSelectorPF.MuonSource = cms.InputTag(muonCollectionPF)
process.RA2muonSelectorPF.VertexSource = cms.InputTag(pvString)


############################# START ntuple specifics ####################################

#in MC the global tag is needed for JEC uncertainies
#in data it is needed for JEC uncertainties and

if isMC:
    print "is MC!"
    process.load('Configuration.StandardSequences.Services_cff')
    process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
    process.GlobalTag.globaltag = 'START38_V14::All' #this Global Tag is for the latest 387
else:
    print "is Data!"
    process.load('Configuration.StandardSequences.Services_cff')
    process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
    process.GlobalTag.globaltag = 'GR_R_38X_V15::All' #this Global Tag is for the latest 387


#flavor history tool
process.load("PhysicsTools.HepMCCandAlgos.flavorHistoryPaths_cfi")

#b-tagging efficiency parameters:
process.load ("RecoBTag.PerformanceDB.PoolBTagPerformanceDBOctEx")
process.load ("RecoBTag.PerformanceDB.BTagPerformanceDBOctEx")
process.load ("RecoBTag.PerformanceDB.PoolBTagPerformanceDBMC36X")
process.load ("RecoBTag.PerformanceDB.BTagPerformanceDBMC36X")

from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
from PhysicsTools.SelectorUtils.jetIDSelector_cfi import jetIDSelector

process.BasicTreeMaker = cms.EDAnalyzer('BasicTreeMaker',
                                        MCflag = cms.bool(isMC),
                                        #b taggers to store in ntuple
                                        btagAlgorithms = cms.vstring("trackCountingHighPurBJetTags","trackCountingHighEffBJetTags",
                                                                     "simpleSecondaryVertexHighEffBJetTags",#"simpleSecondaryVertexNegativeBJetTags",
                                                                     "simpleSecondaryVertexHighPurBJetTags",#),
                                                                     "simpleSecondaryVertexBJetTags"), #for older samples
                                        tauidAlgorithms = cms.vstring("againstElectron", "againstMuon", "byIsolation",
                                                                      "byTaNCfrOnePercent","byTaNCfrHalfPercent", "byTaNCfrQuarterPercent"),
#first the names of various collections
                                        pvTag = cms.InputTag(pvString),
                                        jetAlgorithms = cms.vstring( "cleanPatJetsAK5Calo","selectedPatJetsPF" ), #real collection names
                                        jetNames      = cms.vstring( "calo","PF" ), #for TTree use
                                        metAlgorithms = cms.vstring( "patMETsAK5Calo", "patMETsTC","patMETsPF","patMETsTypeIPF"),
                                        metNames = cms.vstring( "calo", "tc","pf","pf1"),

#careful...the lepton code is pretty sensitive to the input here
                                        eleAlgorithms = cms.vstring(eleCollection,eleCollectionPF),
                                        muonAlgorithms = cms.vstring(muonCollection,muonCollectionPF),
                                        tauAlgorithms = cms.vstring("cleanPatTaus","selectedPatTausPF"),

#careful...these must "line" up with the ele and muon Algorithms given above
                                        eleSelected = cms.vstring("RA2electronSelector","RA2electronSelectorPF"),
                                        muonSelected = cms.vstring("RA2muonSelector","RA2muonSelectorPF"),

#trying to add all particle flow candidates for RA2 filter specifications
                                        PFCandSource = cms.InputTag('particleFlow'),

#the 'SUSY trigger' used in the cut flow is the first (valid) one listed here
                                        triggersOfInterest = cms.vstring("HLT_HT100U","HLT_HT120U","HLT_HT140U","HLT_HT150U","HLT_HT150U_v3","HLT_HT100U_v3","HLT_HT140_J30U_Eta3_v3",
                                                                         "HLT_MET45","HLT_MET60","HLT_MET65","HLT_MET100",
                                                                         "HLT_Jet15U","HLT_Jet30U","HLT_Jet50U","HLT_Jet70U","HLT_Jet100U"),
            
                                        pvSelector = cms.PSet(
                                            #!isFake requirement enforced in PhysicsTools/SelectorUtils/interface/PVSelector.h 
                                            pvSrc = cms.InputTag(pvString),
                                            minNdof = cms.double(4), #5 was the default in 36X but it changed to 4 in 38X; don says use 4
                                            maxZ = cms.double(24), #per Don
                                            maxRho = cms.double(2.0)
                                        ),
                                        jetIdLoose = cms.PSet( #for calo and jpt
                                            version = cms.string('PURE09'), #these are the 384 defaults
                                            quality = cms.string('LOOSE')
                                        ),
                                        jetIdTight = cms.PSet( #for calo and jpt
                                            version = cms.string('PURE09'), #these are the 384 defaults
                                            quality = cms.string('TIGHT')
                                        ),
                                        pfjetIdLoose = cms.PSet( #for pf jets
                                            version = cms.string('FIRSTDATA'), #these are the 384 defaults
                                            quality = cms.string('LOOSE')
                                        ),
                                        pfjetIdTight = cms.PSet( #for pf jets
                                            version = cms.string('FIRSTDATA'), #these are the 384 defaults
                                            quality = cms.string('TIGHT')
                                        ),
                                        
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
                                            RelIso = cms.double(0.1), #uses R=0.3 (at least in 384)

                                            #requires PhysicsTools/SelectorUtils V00-02-24
                                            LepZ = cms.double(1),
                                            nPixelHits = cms.int32(1),
                                            nMatchedStations = cms.int32(1),
                                            
                                            cutsToIgnore = cms.vstring('ED0', 'SD0', 'ECalVeto', 'HCalVeto','LepZ','nPixelHits','minNMatches'),
                                            RecalcFromBeamSpot = cms.bool(True),
                                            beamLineSrc = cms.InputTag("offlineBeamSpot"),
                                            pvSrc = cms.InputTag(pvString)

                                            ),
 
                                        ## electrons
                                        electronId = cms.PSet(
                                            version = cms.string('FIRSTDATA'),
                                            D0 = cms.double(0.20),
                                            ED0 = cms.double(999.0),
                                            SD0 = cms.double(3.0),
                                            RelIso = cms.double( 0.1 ), #uses R=0.3 (at least in 384)
                                            cutsToIgnore = cms.vstring('ED0', 'SD0')
                                            ),

                                        minJets = cms.int32( 3 ),
                                        
                                        jetPtMin = cms.double( 50.0 ),
                                        jetEtaMax = cms.double( 2.4 ),
                                        loosejetPtMin = cms.double( 0.0 ), #no cut
                                        loosejetEtaMax = cms.double( 10.0 ) #no cut

)

process.TFileService = cms.Service("TFileService", fileName = cms.string('BasicNtuple.root') )
#process.TFileService = cms.Service("TFileService", fileName = cms.string('/home/joshmt/BasicNtuple.root') )

if isMC:
    process.p = cms.Path(process.RA2electronSelector*process.RA2electronSelectorPF*process.RA2muonSelector*process.RA2muonSelectorPF
                         *process.flavorHistorySeq
                         *process.BasicTreeMaker)
else:
    process.p = cms.Path(process.RA2electronSelector*process.RA2electronSelectorPF*process.RA2muonSelector*process.RA2muonSelectorPF
                         *process.BasicTreeMaker)

