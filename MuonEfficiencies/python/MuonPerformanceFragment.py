import FWCore.ParameterSet.Config as cms

MuonPerformanceESProducer_1 = cms.ESProducer("MuonPerformanceESProducer",
# this is what it makes available
    ComponentName = cms.string('StandaloneMuon'),
# this is where it gets the payload from                                                
    PayloadName = cms.string('SAMU_T'),
    WorkingPointName = cms.string('SAMU_WP')
)


    





