import FWCore.ParameterSet.Config as cms

#what = "embedded"
what = "MC"

isTest = False

if what=="MC":
#    inputfile = 'file:/uscms_data/d2/joshmt/ttbar-8TeV-Skim_Mu.root'
#    inputfile = 'ttbar-mu-skim.txt'
#    inputfile = 'ttbar-powheg.txt'
    jetcollection = "patJetsAK5PFchs"
    metcollection = "pfMet"
    pfcandcollection = "pfNoPileUpPFchs"
    originalmuon = "promptMuonSelector:mu"
    outfile = "embedTree_MC.root"
elif what=="embedded":
#    inputfile = 'file:/uscms/home/joshmt/nobackup/smallSkims/embedTests/embedFlipReco.justright.root'
#    inputfile = 'ttbar-embed-TauGenTauHad-2012-11-28.txt'
    inputfile = 'ttbar-embedTauGenAllTauHad-2012-12-04.txt'
#    inputfile = 'test-input.txt'
    jetcollection = "ak5PFJetsPFchsMerged"
    metcollection = "pfMetMerged"
    pfcandcollection = "embedmerger:pf"
    originalmuon = "embedsplitter:pfLep"
    outfile = "/uscmst1b_scratch/lpc1/3DayLifetime/joshmt/embedTreeMC_embedTauGenAllTauHad-2012-12-04.root"
#    outfile = "/uscmst1b_scratch/lpc1/3DayLifetime/joshmt/embedTreeMC_embedTauGenTauHad-test.root"
#    outfile = "test.root"

#same for MC and embedded
ttbardecsrc = "ttbarCoder"
genparticlesource = "genParticles::PAT"
vtxsource = "goodVertices"
mdrsrc = "minDRmuonjet"

process = cms.Process("MakeTree")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport  = cms.untracked.PSet(
  reportEvery = cms.untracked.int32(10000)
#  limit = cms.untracked.int32(1)
)

if isTest:
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
else:
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


#input, to be filled in my the SAK magic
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring() #add 'inputfile' here to do tests
)

#load a test file
if isTest and what=="MC":
    process.source.fileNames.extend(cms.vstring('dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/resilient/sakoay/PAT/Summer12_DR53X-PU_S10_START53_V7A-v2/TT_CT10_TuneZ2star_8TeV-powheg-tauola/genLeptonic_lite_v0/patuple_1000_1_utx.root'))

#use steven's tool for reading many files
if what=="embedded":
    from StevenCode.Skimming.sourceFromText import *
    sourceFromText(process, inputfile)

####only for raw MC
if what=="MC":
    # load steven's skim tools
    process.load('StevenCode.Skimming.hadronicSkim_cff')
    # this is how we select the prompt muon for the Control Sample
    process.load('StevenCode.LeptonSelection.promptMuonSelector_cfi')
    process.promptMuonSelector.MuonSource = cms.InputTag("patMuonsPFchs")
    process.promptMuonSelector.ParticleFlowSource = cms.InputTag("pfNoPileUpPFchs")
    process.promptMuonSelector.DoFilter = cms.bool(False) #don't cut events!
    #then form the W
    process.WtoMuNu = cms.EDProducer(
        "CandViewShallowCloneCombiner",
        decay = cms.string('promptMuonSelector:mu pfMet'),
        cut = cms.string(''),
        checkCharge = cms.bool(False)
        )
    process.minDRmuonjet = cms.EDProducer( "MinDeltaRToJet",
                                           pfCandSource = cms.InputTag(""),
                                           GenCandSource = cms.InputTag("ttbarCoder:genLeptons"),
                                           JetSource = cms.InputTag("goodJets")
                                           )
 

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
                                 VertexSource = cms.InputTag(vtxsource),
                                 PFSource = cms.InputTag(pfcandcollection),
                                 DoMuonVeto = cms.bool(False)
                                 )
process.vetoElectrons = cms.EDFilter("SAKLooseElectronSelector",
                                 ElectronSource = cms.InputTag("gsfElectrons"),
                                 VertexSource = cms.InputTag(vtxsource),
                                 PFSource = cms.InputTag(pfcandcollection),
                                 DoElectronVeto = cms.bool(False)
                                 )


process.maketree = cms.EDAnalyzer('EmbedTree',
                                  VtxSource = cms.InputTag(vtxsource),
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
    process.p = cms.Path(
        process.goodVertices #calculate good vertices
        * process.promptMuonSelector #find control sample muons
        * process.WtoMuNu #calculate W pt
        * process.goodJets #make a goodJets collection
        * process.ttbarCoder
        * process.minDRmuonjet
        *process.trackIsolationMaker*process.vetoMuons*process.vetoElectrons*process.maketree
        )
elif what=="embedded":
    process.p = cms.Path(process.ttbarCoder*process.trackIsolationMaker*process.vetoMuons*process.vetoElectrons*process.maketree)
    
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

