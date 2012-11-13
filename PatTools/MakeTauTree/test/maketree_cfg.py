import FWCore.ParameterSet.Config as cms

process = cms.Process("MakeTauTree")


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#    'file:/eos/uscms/store/user/rossin/PAT/Summer12_DR53X-PU_S10_START53_V7A-v1/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/v0/patuple_455_1_btE.root'
    'root://cmseos:1094//eos/uscms/store/user/nmccoll/PAT/Summer12_DR53X-PU_S10_START53_V7A-v2/TT_CT10_TuneZ2star_8TeV-powheg-tauola/lite_v0/patuple_648_1_0XS.root'  
#    'file:/uscms/home/joshmt/nobackup/ttbar-8TeV-Skim_Mu.root'
#    'file:/afs/cern.ch/work/j/joshmt/private/stop7TeV/CMSSW_4_2_5/src/SandBox/Skims/test/TTJets_Fall11_PAT.root'
#    'file:/afs/cern.ch/user/j/joshmt/work/public/stop7TeV/CMSSW_4_2_5/src//SandBox/Skims/test/PAToutput.root'
#    'file:/afs/cern.ch/user/j/joshmt/work/public/stop7TeV-v1/CMSSW_4_2_5/src/SandBox/Skims/test/PAToutput.root'
#    'file:/afs/cern.ch/user/j/joshmt/work/public/stop8TeV-v1/CMSSW_5_2_5/src/SandBox/Skims/test/susypat.root'
#    'file:/afs/cern.ch/user/j/joshmt/work/public/stop8TeV/CMSSW_5_2_5/src/SandBox/Skims/test/susypat.root'
#    'file:/cu3/joshmt/stop/PAT/TTJets_Fall11_PAT_v2.root'
    )
                            )

#need to run the directional isolation e/mu
process.VetoElectrons = cms.EDFilter("SAKLoosePATElectronSelector",
                                     ElectronSource = cms.InputTag("patElectronsPFchs"),
                                     VertexSource = cms.InputTag("goodVertices"),
                                     PFSource = cms.InputTag("pfNoPileUpPFchs"), #correct?
                                     DoID = cms.bool(True),
                                     DoIso = cms.bool(True),
                                     DoVeto = cms.bool(False)
                                     )
process.VetoMuons = cms.EDFilter("SAKLoosePATMuonSelector",
                                 MuonSource = cms.InputTag("patMuonsPFchs"),
                                 VertexSource = cms.InputTag("goodVertices"),
                                 PFSource = cms.InputTag("pfNoPileUpPFchs"), #correct?
                                 DoID = cms.bool(True),
                                 DoIso = cms.bool(True),
                                 DoVeto = cms.bool(False)
                                 )


#need to run the indirect tau veto code too
#prepare a jet collection with pT and eta cuts
process.tauJetCands = cms.EDFilter("PATJetSelector",
                                   src = cms.InputTag("patJetsAK5PFchs"),
                                   cut = cms.string("pt > 15 & abs( eta ) < 2.4")
                                   )
process.VetoTaus = cms.EDFilter("IndirectTauSelector",
                                JetSource = cms.InputTag("tauJetCands"),
                                METSource = cms.InputTag("pfMet"),
                                DoVeto = cms.bool(False)
                                )


#finally the isolated track finder
process.trackIsolationSelector = cms.EDFilter("IsolatedTrackSelector",
                                             PFSource     = cms.InputTag("pfNoPileUpPFchs"),
                                             VertexSource      = cms.InputTag("goodVertices"),
                                             DR         = cms.double(0.3),
                                             DzMax         = cms.double(0.05),
                                             IsoMax = cms.double(0.1),
                                             PtMin = cms.double(10),
                                             DoVeto = cms.bool(False)
                                             )

#ben hooberman's version here
process.trackIsolationMaker =  cms.EDProducer("TrackIsolationMaker",
                                            pfCandidatesTag     = cms.InputTag("pfNoPileUpPFchs"),
                                            vertexInputTag      = cms.InputTag("goodVertices"),
                                            dR_ConeSize         = cms.double(0.3),
                                            dz_CutValue         = cms.double(0.05),
                                            minPt_PFCandidate   = cms.double(1)
                                            )

process.ttbarDecayProducer = cms.EDProducer('TTbarDecayCoder',
                                            GenParticleSource = cms.InputTag('genParticles')
)


process.MakeTauTree = cms.EDAnalyzer('MakeTauTree',
                                     TauSource = cms.InputTag("patTausPFchs"), #selectedTaus (these are POG taus)
                                     JetSource = cms.InputTag("patJetsAK5PFchs"),
                                     MetSource = cms.InputTag("pfMet"), # reco::PFMET
                                     ElectronVetoSource = cms.InputTag("VetoElectrons"),
                                     MuonVetoSource = cms.InputTag("VetoMuons"),
                                     TauVetoSource = cms.InputTag("VetoTaus"),
#                                     GenParticleSource = cms.InputTag("genParticles"), #prunedGenParticles
                                     VertexSource = cms.InputTag("goodVertices"),
                                     StoreAllJetInfo = cms.bool(False)
                                     )

#process.TFileService = cms.Service("TFileService", fileName = cms.string('/cu3/joshmt/stop/reducedTrees/reducedTree.TTbarJets.root') )
process.TFileService = cms.Service("TFileService", fileName = cms.string('reducedTree.root') )

#stupid that i'm basically running the track isolation code twice. but let's leave well enough alone
process.p = cms.Path(process.VetoElectrons*process.VetoMuons*process.tauJetCands*process.VetoTaus*process.trackIsolationSelector*process.trackIsolationMaker*process.ttbarDecayProducer
                     *process.MakeTauTree)
