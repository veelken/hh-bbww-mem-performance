import FWCore.ParameterSet.Config as cms
import os

#--------------------------------------------------------------------------------
# CV: imports needed by analyzeConfig.py base-class
from tthAnalysis.HiggsToTauTau.configs.recommendedMEtFilters_cfi import *
from tthAnalysis.HiggsToTauTau.configs.EvtYieldHistManager_cfi import *
from tthAnalysis.HiggsToTauTau.configs.hhWeight_cfi import hhWeight
#--------------------------------------------------------------------------------

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames = cms.vstring(),
    maxEvents = cms.int32(-1),
    outputEvery = cms.uint32(100)
)

process.fwliteOutput = cms.PSet(
    fileName = cms.string('')
)

process.analyze_hh_bbwwMEM_singlelepton = cms.PSet(
    treeName = cms.string('Events'),

    skipSelEvents = cms.int32(0),
    maxSelEvents = cms.int32(1000),

    process = cms.string(''),
    histogramDir = cms.string(''),
    era = cms.string('2017'),

    apply_jetSmearing = cms.bool(True),
    jetSmearing_coeff = cms.double(1.00),
    apply_metSmearing = cms.bool(True),
    metSmearing_sigmaX = cms.double(10.),
    metSmearing_sigmaY = cms.double(10.),

    apply_genWeight = cms.bool(True),
    hasLHE = cms.bool(True),

    branchName_genLeptons = cms.string('GenLep'),
    branchName_genNeutrinos = cms.string('GenNu'),
    branchName_genJets = cms.string('GenJet'),

    # branches specific to HH signal
    branchName_genParticlesFromHiggs = cms.string('GenHiggsDaughters'),
    branchName_genWBosons = cms.string('GenVbosons'),
    branchName_genWJets = cms.string('GenWZQuark'),

    # branches specific to ttbar background
    branchName_genLeptonsFromTop = cms.string('GenLepFromTop'),
    branchName_genNeutrinosFromTop = cms.string('GenNuFromTop'),
    branchName_genBQuarksFromTop = cms.string('GenBQuarkFromTop'),
    branchName_genWJetsFromTop = cms.string('GenQuarkFromTop'),

    selEventsFileName_input = cms.string(''),
    selEventsFileName_output = cms.string(''),
  
    # general configuration parameters, required by our analysis framework
    leptonFakeRateWeight = cms.PSet(),
    hhWeight_cfg = hhWeight,
    gen_mHH = cms.vdouble(),
    nonRes_BMs = cms.vstring(),

    isDEBUG = cms.bool(False)
)