import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntupler")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
process.load("Geometry.CaloEventSetup.CaloTopology_cfi");

# Define input data to read
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
inputFilesMiniAOD = cms.untracked.vstring(
    'file:/wk_cms/pchen/work/HWAnalysis/data/eos/cms/store/data/Run2016B/SingleElectron/MINIAOD/23Sep2016-v2/80000/08A02DC3-608C-E611-ADA5-0025905B85B6.root'
)
process.source = cms.Source ("PoolSource", fileNames = inputFilesMiniAOD )

#process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')

# Set up electron ID (VID framework)

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
                 #'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',
                 #'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

#process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('calibratedPatElectrons')
#process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('selectedElectrons')

# Configure the ntupler module
process.selectedElectrons = cms.EDFilter("PATElectronSelector", 
    src = cms.InputTag("slimmedElectrons"), 
    cut = cms.string("pt > 9 && abs(eta)<2.5")
)

# should have no effect on slimmedElectrons
#process.calibratedPatElectrons.isMC = cms.bool(True)

process.ntupler = cms.EDAnalyzer('runEffAnaMiniAOD',
    isMC     = cms.untracked.bool(False),
    xsec     = cms.untracked.double(1.),
    dtag     = cms.untracked.string("SingleElectron"),

    trigger  = cms.InputTag("TriggerResults::HLT"),
    prescale = cms.InputTag("patTrigger"),
    # pileup   = cms.InputTag("slimmedAddPileupInfo"),
    rho      = cms.InputTag("fixedGridRhoFastjetAll"),
    beamSpot = cms.InputTag('offlineBeamSpot'),
    eventWeight   = cms.InputTag("generator"),

    # Objects specific to MiniAOD format
    muonsMiniAOD    = cms.InputTag("slimmedMuons"),
    electronsMiniAOD    = cms.InputTag("slimmedElectrons"),
    photonsMiniAOD      = cms.InputTag("slimmedPhotons"),
    metsMiniAOD         = cms.InputTag("slimmedMETs"),
    genParticlesMiniAOD = cms.InputTag("prunedGenParticles"),
    verticesMiniAOD     = cms.InputTag("offlineSlimmedPrimaryVertices"),
    conversionsMiniAOD  = cms.InputTag('reducedEgamma:reducedConversions'),

    # ID decisions (common to all formats)
    eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
    eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
    eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
    eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
    elenontrigMVAlooseIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90"),
    elenontrigMVAtightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80"),
    phoLooseIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-loose"),
    phoMediumIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-medium"),
    phoTightIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-tight"),
    # ValueMaps with MVA results
    mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
    mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories"),

    # This is a fairly verbose mode if switched on, with full cut flow 
    # diagnostics for each candidate. Use it in a low event count test job.
    eleMediumIdFullInfoMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium",),
    eleIdVerbose = cms.bool(False),
    objects = cms.InputTag('selectedPatTrigger')
)

process.primaryVertexFilter  = cms.EDFilter("VertexSelector",
      src = cms.InputTag('offlineSlimmedPrimaryVertices'),
      cut = cms.string('!isFake && ndof > 4.0 && position.Rho < 2.0 && abs(z) < 24'),
      filter = cms.bool(True)  ## otherwise it won't filter the events, just produce an empty vertex collection.
      )

process.hcalDDDRecConstants = cms.ESProducer( "HcalDDDRecConstantsESModule",
  appendToDataLabel = cms.string( "" )
)
process.hcalDDDSimConstants = cms.ESProducer( "HcalDDDSimConstantsESModule",
  appendToDataLabel = cms.string( "" )
)

#from PhysicsTools.PatAlgos.preselection_eltau_cff import *
#process.leadingElectron = process.slimmedElectrons.clone(
#    src = 'slimmedElectrons', cut = 'pt>20'
#    )
#process.load("PhysicsTools.PatAlgos.selectionLayer1.electronCountFilter_cfi")
#process.leadingElectronRequirement = process.countPatElectrons.clone(minNumber = 1, src = 'leadingElectron')
#process.leadingleptonsequence = cms.Sequence(process.leadingElectron+process.leadingElectronRequirement)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('MC.root')
)

#process.p = cms.Path(process.selectedElectrons + process.calibratedPatElectrons + process.egmGsfElectronIDSequence + process.primaryVertexFilter + process.ntupler) 
#process.p = cms.Path(process.selectedElectrons * process.egmGsfElectronIDSequence * process.primaryVertexFilter * process.ntupler)
process.p = cms.Path(process.egmGsfElectronIDSequence * process.primaryVertexFilter * process.ntupler)
