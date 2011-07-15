import FWCore.ParameterSet.Config as cms

trackMPT = cms.EDProducer('MPT',
                       minpt = cms.double(5)
                          )
