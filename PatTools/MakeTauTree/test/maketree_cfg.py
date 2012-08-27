import FWCore.ParameterSet.Config as cms

process = cms.Process("MakeTauTree")


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#    'file:/afs/cern.ch/work/j/joshmt/private/stop7TeV/CMSSW_4_2_5/src/SandBox/Skims/test/TTJets_Fall11_PAT.root'
#    'file:/afs/cern.ch/user/j/joshmt/work/public/stop7TeV/CMSSW_4_2_5/src//SandBox/Skims/test/PAToutput.root'
#    'file:/afs/cern.ch/user/j/joshmt/work/public/stop7TeV-v1/CMSSW_4_2_5/src/SandBox/Skims/test/PAToutput.root'
    'file:/afs/cern.ch/user/j/joshmt/work/public/stop8TeV-v1/CMSSW_5_2_5/src/SandBox/Skims/test/susypat.root'
#    'file:/afs/cern.ch/user/j/joshmt/work/public/stop8TeV/CMSSW_5_2_5/src/SandBox/Skims/test/susypat.root'
#    'file:/cu3/joshmt/stop/PAT/TTJets_Fall11_PAT_v2.root'
    )
                            )


process.MakeTauTree = cms.EDAnalyzer('MakeTauTree',
                                     TauSource = cms.InputTag("selectedTaus"),
                                     JetSource = cms.InputTag("selectedJets"),
                                     MetSource = cms.InputTag("patMETsPF"), # chs is or is not required
                                     ElectronVetoSource = cms.InputTag("VetoElectrons"),
                                     MuonVetoSource = cms.InputTag("VetoMuons"),
                                     TauVetoSource = cms.InputTag("VetoTaus"),
                                     GenParticleSource = cms.InputTag("prunedGenParticles")
                                     )

#process.TFileService = cms.Service("TFileService", fileName = cms.string('/cu3/joshmt/stop/reducedTrees/reducedTree.TTbarJets.root') )
process.TFileService = cms.Service("TFileService", fileName = cms.string('reducedTree.root') )

process.p = cms.Path(process.MakeTauTree)
