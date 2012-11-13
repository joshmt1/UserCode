import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'root://cmseos:1094//eos/uscms/store/user/nmccoll/PAT/Summer12_DR53X-PU_S10_START53_V7A-v2/TT_CT10_TuneZ2star_8TeV-powheg-tauola/lite_v0/patuple_648_1_0XS.root'
#        'file:/cu1/joshmt/Validation/store__mc__Summer12__TTJets_MassiveBinDECAY__AODSIM.root' #52X ttbar
    )
)

process.myProducerLabel = cms.EDProducer('TTbarDecayCoder',
                                         GenParticleSource = cms.InputTag('genParticles')

)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)

  
process.p = cms.Path(process.myProducerLabel)

process.e = cms.EndPath(process.out)
