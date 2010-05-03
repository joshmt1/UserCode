import FWCore.ParameterSet.Config as cms

SupplementaryNtuple = cms.EDAnalyzer('SupplementaryNtuple',
                                     genParticles = cms.InputTag("genParticles")
                                     )
