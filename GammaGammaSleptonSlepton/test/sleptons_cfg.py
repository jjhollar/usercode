import FWCore.ParameterSet.Config as cms

process = cms.Process("gamgam2leplepanalysis")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load("Configuration.EventContent.EventContent_cff")

process.source = cms.Source("PoolSource",
    debugFlag = cms.untracked.bool(False),
    fileNames = cms.untracked.vstring('file:/scratch/jjhollar/Madgraph_gamgamsmusmu_RECO.root')
)

process.smuanal = cms.EDFilter("GammaGammaSleptonSlepton",
    outfilename = cms.untracked.string('gamgamsmusmufp420.anal.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False)
)
process.p = cms.Path(process.smuanal)


