import FWCore.ParameterSet.Config as cms

#what = "embedded"
what = "MC"

if what=="MC":
#    inputfile = 'file:/uscms_data/d2/joshmt/ttbar-8TeV-Skim_Mu.root'
#    inputfile = 'ttbar-mu-skim.txt'
#    inputfile = 'ttbar-powheg.txt'
    inputfile = 'localttbar.txt'
    jetcollection = "patJetsAK5PFchs"
    metcollection = "pfMet"
    pfcandcollection = "pfNoPileUpPFchs"
    ttbardecsrc = "ttbarCoder"
    mdrsrc = ""
    originalmuon = ""
    outfile = "embedTree_MC.root"
elif what=="embedded":
#    inputfile = 'file:/uscms/home/joshmt/nobackup/smallSkims/embedTests/embedFlipReco.justright.root'
    inputfile = 'ttbar-embedFlipReco.txt'
    jetcollection = "ak5PFJetsPFchsMerged"
    metcollection = "pfMetMerged"
    pfcandcollection = "embedmerger"
    ttbardecsrc = ""
    mdrsrc = "minDRmuonjet"
    originalmuon = "embedsplitter"
    outfile = "embedTree_MC_embedFlipReco.root"

#same for MC and embedded
genparticlesource = "genParticles"
vtxsource = "goodVertices"

process = cms.Process("MakeTree")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport  = cms.untracked.PSet(
  reportEvery = cms.untracked.int32(10000)
#  limit = cms.untracked.int32(1)
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#input, to be filled in my the SAK magic
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring() #add 'inputfile' here to do tests
)

#use steven's tool for reading many files
from StevenCode.Skimming.sourceFromText import *
sourceFromText(process, inputfile)


process.ttbarCoder = cms.EDProducer("TTbarDecayCoder",
                                    GenParticleSource = cms.InputTag(genparticlesource)
                                    )

process.trackIsolationMaker =  cms.EDProducer("TrackIsolationMaker",
                                              pfCandidatesTag     = cms.InputTag(pfcandcollection),
                                              vertexInputTag      = cms.InputTag(vtxsource),
                                              dR_ConeSize         = cms.double(0.3),
                                              dz_CutValue         = cms.double(0.05),
                                              minPt_PFCandidate   = cms.double(5)
                                              )

#directional iso
process.vetoMuons = cms.EDFilter("SAKLooseMuonSelector",
                                 MuonSource = cms.InputTag("muons"),
                                 VertexSource = cms.InputTag(vtxsource), #to be checked
                                 PFSource = cms.InputTag(pfcandcollection),
                                 DoMuonVeto = cms.bool(False)
                                 )
process.vetoElectrons = cms.EDFilter("SAKLooseElectronSelector",
                                 ElectronSource = cms.InputTag("gsfElectrons"),
                                 VertexSource = cms.InputTag(vtxsource), #to be checked
                                 PFSource = cms.InputTag(pfcandcollection),
                                 DoElectronVeto = cms.bool(False)
                                 )


process.maketree = cms.EDAnalyzer('EmbedTree',
                                  JetSource = cms.InputTag(jetcollection),
                                  MetSource = cms.InputTag(metcollection),
                                  TtbarDecaySource = cms.InputTag(ttbardecsrc),
                                  MuonSource = cms.InputTag("vetoMuons"),
                                  ElectronSource = cms.InputTag("vetoElectrons"),
                                  MinDRSource = cms.InputTag(mdrsrc),
                                  OriginalMuonSource = cms.InputTag(originalmuon)
                                  )

process.TFileService = cms.Service("TFileService", fileName = cms.string(outfile) )

if what=="MC":
    process.p = cms.Path(process.ttbarCoder
                         *process.trackIsolationMaker
                         *process.vetoMuons
                         *process.vetoElectrons
                         *process.maketree)
elif what=="embedded":
    process.p = cms.Path(process.trackIsolationMaker
                         *process.vetoMuons
                         *process.vetoElectrons
                         *process.maketree)
    
#magic from SAK
def InputConfiguration(postfix, process):
  import os
  allNames                = vars(process).keys()
  for name in allNames:
    object                = getattr(process, name)
    try:
      if object.type_() == 'PoolOutputModule':
        (name, extension) = os.path.splitext(object.fileName.value())
        object.fileName   = "%s%s%s" % (name, postfix, extension)
    except AttributeError:
      pass
    try:
      if object.type_() == 'TFileService':
        (name, extension) = os.path.splitext(object.fileName.value())
        object.fileName   = "%s%s%s" % (name, postfix, extension)
    except AttributeError:
      pass
        

InputConfiguration("", process)
