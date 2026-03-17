# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: -s L1,L1P2GT --conditions auto:phase2_realistic_T35 --geometry ExtendedRun4D121 --era Phase2C17I13M9 --eventcontent FEVTDEBUGHLT --datatier GEN-SIM-DIGI-RAW-MINIAOD --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000,Configuration/DataProcessing/Utils.addMonitoring,L1Trigger/Configuration/customisePhase2.addHcalTriggerPrimitives --filein /store/mc/Phase2Spring24DIGIRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_140X_mcRun4_realistic_v4-v1/2810000/5bed6cad-6cdb-4a5e-87af-08cc8c8e0ff8.root --fileout file:output_Phase2_L1T.root --python_filename rerunL1_cfg.py '--inputCommands=keep *, drop l1tPFJets_*_*_*, drop l1tTrackerMuons_l1tTkMuonsGmt*_*_HLT' --mc -n 10 --nThreads 1
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

process = cms.Process('L1P2GT2',Phase2C17I13M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtendedRun4D121Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimPhase2L1GlobalTriggerEmulator_cff')
process.load('L1Trigger.Configuration.Phase2GTMenus.SeedDefinitions.step1_2024.l1tGTMenu_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('L1Trigger.Phase2L1GMT.gmtKMTFMuons_cfi')
process.load('L1Trigger.Phase2L1GMT.gmtStubs_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(50),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    fileNames = cms.untracked.vstring('file:output_Phase2_L1T.root'),
    inputCommands = cms.untracked.vstring(
        'keep *'
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    TryToContinue = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToCallForTryToContinue = cms.untracked.vstring(),
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('-s nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
process.FEVTDEBUGHLTEventContent.outputCommands=cms.untracked.vstring(
    'drop *',
    'keep *_simDtTriggerPrimitiveDigis_*_*',
    "keep *_dtTriggerPhase2PrimitivePairDigis_*_*",
    'keep *_genParticles_*_*',
    'keep *_dtTriggerPhase2PrimitiveDigis_*_*',
    "keep *_simMuonDTDigis_*_*",
    "keep *_gmtStubs_*_*",
    "keep *_gmtKMTFMuons_*_*",
    "keep *_l1tSAMuonsGmt_*_*",
    "keep *_l1tKMTFMuonsGmt_*_*",
    "keep *_l1tStubsGmt_*_*",
    "keep *_simBmtfDigis_*_*",
    "keep *_dtTriggerPhase2PrimitivePairDigis_*_*",
    "keep *_simCscTriggerPrimitiveDigis_*_*",
    "keep *_simMuonRPCDigis_*_*")
process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW-MINIAOD'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:output_DY_Phase2_L1T_KMTF.root'),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T35', '')

# Path and EndPath definitions
process.GMTStubProducer_offline = cms.Path(process.gmtStubs)
process.GMTKMTFMuons_offline = cms.Path(process.gmtKMTFMuons)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
process.schedule = cms.Schedule(process.GMTStubProducer_offline, process.GMTKMTFMuons_offline, process.endjob_step, process.FEVTDEBUGHLToutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.aging
from SLHCUpgradeSimulations.Configuration.aging import customise_aging_1000 

#call to customisation function customise_aging_1000 imported from SLHCUpgradeSimulations.Configuration.aging
process = customise_aging_1000(process)

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# Automatic addition of the customisation function from L1Trigger.Configuration.customisePhase2
from L1Trigger.Configuration.customisePhase2 import addHcalTriggerPrimitives 

#call to customisation function addHcalTriggerPrimitives imported from L1Trigger.Configuration.customisePhase2
process = addHcalTriggerPrimitives(process)

# End of customisation functions


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
