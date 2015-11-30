import FWCore.ParameterSet.Config as cms

#from FWCore.ParameterSet.VarParsing import VarParsing
###############################
####### Parameters ############
###############################
#options = VarParsing ('python')
#options.parseArguments()

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
	                    fileNames = cms.untracked.vstring('file:AOD_data.root')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.TFileService=cms.Service(
    "TFileService",
    fileName=cms.string('testoutput_thrust.root')
)

process.ThrustAKPu4PF = cms.EDAnalyzer(
    'ThrustAnalyzer',
    JetType = cms.untracked.string('pf'),
    UEAlgo = cms.untracked.string('Pu'),
    radius = cms.int(4),
    src = cms.InputTag("akVs4PFJets"),
    #centralitycollection = cms.InputTag("hiCentrality"),
    #centralitybincollection = cms.InputTag("centralityBin","HFtowers"),
    #JetCorrections = cms.string(""),
    #PFcands = cms.InputTag("particleFlowTmp"),
    RThreshold = cms.double(0.3),  
    mRecoJetPtThreshold = cms.double(10),        
    mReco_SubJetPtThreshold = cms.double(10)        
)


process.ThrustAKVs4PF = process.ThrustAKVs4PF.clone(
    UEAlgo = cms.untracked.string('Pu')
)


process.p = cms.Path(
    process.ThrustAKPu4PF *
    process.ThrustAKVs4PF
)

process.MessageLogger.cerr.FwkReport.reportEvery = 10000
