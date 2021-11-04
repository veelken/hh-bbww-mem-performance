#include "FWCore/ParameterSet/interface/ParameterSet.h" // edm::ParameterSet
#include "FWCore/Utilities/interface/Exception.h" // cms::Exception
#include "PhysicsTools/FWLite/interface/TFileService.h" // fwlite::TFileService
#include "DataFormats/FWLite/interface/InputSource.h" // fwlite::InputSource
#include "DataFormats/FWLite/interface/OutputFiles.h" // fwlite::OutputFiles
#include "DataFormats/Math/interface/LorentzVector.h" // math::PtEtaPhiMLorentzVector
#include "DataFormats/Math/interface/deltaR.h" // deltaR
#include "DataFormats/Math/interface/deltaPhi.h" // deltaPhi

#if __has_include (<FWCore/ParameterSetReader/interface/ParameterSetReader.h>)
#  include <FWCore/ParameterSetReader/interface/ParameterSetReader.h> // edm::readPSetsFrom()
#else
#  include <FWCore/PythonParameterSet/interface/MakeParameterSets.h> // edm::readPSetsFrom()
#endif

#include <TBenchmark.h> // TBenchmark
#include <TString.h> // TString, Form
#include <TError.h> // gErrorAbortLevel, kError
#include <TRandom3.h> // TRandom3
#include <TLorentzVector.h> // TLorentzVector 
#include <TMatrixD.h> // TMatrixD

#include "tthAnalysis/HiggsToTauTau/interface/GenLepton.h" // GenLepton
#include "tthAnalysis/HiggsToTauTau/interface/GenJet.h" // GenJet
#include "tthAnalysis/HiggsToTauTau/interface/EventInfo.h" // EventInfo
#include "tthAnalysis/HiggsToTauTau/interface/GenLeptonReader.h" // GenLeptonReader
#include "tthAnalysis/HiggsToTauTau/interface/GenLeptonCollectionSelector.h" // GenLeptonCollectionSelector
#include "tthAnalysis/HiggsToTauTau/interface/GenParticleReader.h" // GenParticleReader
#include "tthAnalysis/HiggsToTauTau/interface/GenJetReader.h" // GenJetReader
#include "tthAnalysis/HiggsToTauTau/interface/GenJetCollectionSelector.h" // GenJetCollectionSelector
#include "tthAnalysis/HiggsToTauTau/interface/LHEInfoReader.h" // LHEInfoReader
#include "tthAnalysis/HiggsToTauTau/interface/EventInfoReader.h" // EventInfoReader
#include "tthAnalysis/HiggsToTauTau/interface/convert_to_ptrs.h" // convert_to_ptrs
#include "tthAnalysis/HiggsToTauTau/interface/ParticleCollectionCleaner.h" // GenJetCollectionCleaner
#include "tthAnalysis/HiggsToTauTau/interface/RunLumiEventSelector.h" // RunLumiEventSelector
#include "tthAnalysis/HiggsToTauTau/interface/CutFlowTableHistManager.h" // CutFlowTableHistManager
#include "tthAnalysis/HiggsToTauTau/interface/WeightHistManager.h" // WeightHistManager
#include "tthAnalysis/HiggsToTauTau/interface/GenEvtHistManager.h" // GenEvtHistManager
#include "tthAnalysis/HiggsToTauTau/interface/LHEInfoHistManager.h" // LHEInfoHistManager
#include "tthAnalysis/HiggsToTauTau/interface/leptonTypes.h" // getLeptonType, kElectron, kMuon
#include "tthAnalysis/HiggsToTauTau/interface/analysisAuxFunctions.h" // isHigherPt, isMatched, contains, findFile
#include "tthAnalysis/HiggsToTauTau/interface/generalAuxFunctions.h" // format_vstring
#include "tthAnalysis/HiggsToTauTau/interface/cutFlowTable.h" // cutFlowTableType
#include "tthAnalysis/HiggsToTauTau/interface/TTreeWrapper.h" // TTreeWrapper
#include "tthAnalysis/HiggsToTauTau/interface/hltFilter.h" // hltFilter()

#include <boost/math/special_functions/sign.hpp> // boost::math::sign()
#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/MEMbbwwHistManager.h" // MEMbbwwHistManager
#include "tthAnalysis/HiggsToTauTau/interface/LocalFileInPath.h" // LocalFileInPath
#include "hhAnalysis/bbwwMEM/interface/MEMbbwwAlgoSingleLepton.h"
#include "hhAnalysis/bbwwMEM/interface/MeasuredParticle.h" // MeasuredParticle
#include "hhAnalysis/bbwwMEM/interface/memAuxFunctions.h"
#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/BJetTF_toy.h" // BJetTF_toy
#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/HadWJetTF_toy.h" // HadWJetTF_toy
#include "hhAnalysis/bbww/interface/genMatchingAuxFunctions.h" // findGenLepton_and_NeutrinoFromWBoson
#include "tthAnalysis/HiggsToTauTau/interface/histogramAuxFunctions.h" // fillWithOverFlow()
#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/GenJetSmearer.h" // GenJetSmearer
#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/GenMEtSmearer.h" // GenMEtSmearer
#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/MEMEvent_singlelepton.h" // MEMEvent_singlelepton
#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/MEMbbwwNtupleManager_singlelepton.h" // MEMbbwwNtupleManager_singlelepton
#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/memNtupleAuxFunctions_singlelepton.h" // addGenMatches_singlelepton
#include "hhAnalysis/multilepton/interface/AnalysisConfig_hh.h" // AnalysisConfig_hh

#include <iostream> // std::cerr, std::fixed
#include <iomanip> // std::setprecision(), std::setw()
#include <string> // std::string
#include <vector> // std::vector<>
#include <cstdlib> // EXIT_SUCCESS, EXIT_FAILURE
#include <fstream> // std::ofstream
#include <assert.h> // assert

typedef math::PtEtaPhiMLorentzVector LV;
typedef std::vector<std::string> vstring;
typedef std::vector<double> vdouble;

//std::string
//format_vsize(const std::set<size_t> & vs)
//{
//  std::ostringstream os;
//  const unsigned numEntries = vs.size();
//
//  os << "{ ";
//  unsigned iEntry = 0;
//  for ( std::set<size_t>::const_iterator vs_iter = vs.begin(); vs_iter != vs.end(); ++vs_iter )
//  {
//    os << (*vs_iter);
//    if(iEntry < numEntries - 1)
//    {
//      os << ", ";
//    }
//    ++iEntry;
//  }
//  os << " }";
//
//  return os.str();
//}

/**
 * @brief Produce MEM performance plots for single-lepton channel of the HH->bbWW analysis.
 */
int main(int argc, char* argv[])
{
//--- throw an exception in case ROOT encounters an error
  gErrorAbortLevel = kError;

//--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "<analyze_hh_bbwwMEM_singlelepton>:" << std::endl;

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("analyze_hh_bbwwMEM_singlelepton");

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") )
    throw cms::Exception("analyze_hh_bbwwMEM_singlelepton")
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfg_analyze = cfg.getParameter<edm::ParameterSet>("analyze_hh_bbwwMEM_singlelepton");
  AnalysisConfig_hh analysisConfig("HH->bbWW", cfg_analyze);

  std::string treeName = cfg_analyze.getParameter<std::string>("treeName");

  std::string process_string = cfg_analyze.getParameter<std::string>("process");
  bool isSignal = ( process_string.find("signal") != std::string::npos ) ? true : false;

  std::string histogramDir = cfg_analyze.getParameter<std::string>("histogramDir");
  
  std::string era_string = cfg_analyze.getParameter<std::string>("era");
  const Era era = get_era(era_string);

  std::string branchName_genLeptons = cfg_analyze.getParameter<std::string>("branchName_genLeptons");
  std::string branchName_genNeutrinos = cfg_analyze.getParameter<std::string>("branchName_genNeutrinos");
  std::string branchName_genJets = cfg_analyze.getParameter<std::string>("branchName_genJets");

  // branches specific to HH signal
  std::string branchName_genParticlesFromHiggs = cfg_analyze.getParameter<std::string>("branchName_genParticlesFromHiggs");
  std::string branchName_genWBosons = cfg_analyze.getParameter<std::string>("branchName_genWBosons");
  std::string branchName_genWJets = cfg_analyze.getParameter<std::string>("branchName_genWJets");

  // branches specific to ttbar background
  std::string branchName_genLeptonsFromTop = cfg_analyze.getParameter<std::string>("branchName_genLeptonsFromTop");
  std::string branchName_genNeutrinosFromTop = cfg_analyze.getParameter<std::string>("branchName_genNeutrinosFromTop");
  std::string branchName_genBQuarksFromTop = cfg_analyze.getParameter<std::string>("branchName_genBQuarksFromTop");
  std::string branchName_genWJetsFromTop = cfg_analyze.getParameter<std::string>("branchName_genWJetsFromTop");

  bool apply_jetSmearing = cfg_analyze.getParameter<bool>("apply_jetSmearing");
  double jetSmearing_coeff = cfg_analyze.getParameter<double>("jetSmearing_coeff");
  GenJetSmearer genJetSmearer;
  genJetSmearer.set_coeff(jetSmearing_coeff);
  mem::BJetTF_toy bjetTF;
  bjetTF.set_coeff(jetSmearing_coeff);
  mem::HadWJetTF_toy hadWJetTF;
  hadWJetTF.set_coeff(jetSmearing_coeff);

  bool apply_metSmearing = cfg_analyze.getParameter<bool>("apply_metSmearing");
  double metSmearing_sigmaX = cfg_analyze.getParameter<double>("metSmearing_sigmaX");
  double metSmearing_sigmaY = cfg_analyze.getParameter<double>("metSmearing_sigmaY");
  GenMEtSmearer genMEtSmearer;
  genMEtSmearer.set_sigmaX(metSmearing_sigmaX);
  genMEtSmearer.set_sigmaY(metSmearing_sigmaY);
  TMatrixD metCov(2,2);
  metCov[0][0] = mem::square(metSmearing_sigmaX);
  metCov[1][0] = 0.;
  metCov[0][1] = 0.;
  metCov[1][1] = mem::square(metSmearing_sigmaY);

  // random number generator for choosing fake b-jets
  TRandom3 rnd;
  rnd.SetSeed(12345);

  const std::string central_or_shift = "central";
  bool hasLHE = cfg_analyze.getParameter<bool>("hasLHE");
  bool apply_genWeight = cfg_analyze.getParameter<bool>("apply_genWeight");

  bool isDEBUG = cfg_analyze.getParameter<bool>("isDEBUG");

  std::string selEventsFileName_input = cfg_analyze.getParameter<std::string>("selEventsFileName_input");
  std::cout << "selEventsFileName_input = " << selEventsFileName_input << std::endl;
  RunLumiEventSelector* run_lumi_eventSelector = 0;
  if ( selEventsFileName_input != "" ) {
    edm::ParameterSet cfg_runLumiEventSelector;
    cfg_runLumiEventSelector.addParameter<std::string>("inputFileName", selEventsFileName_input);
    cfg_runLumiEventSelector.addParameter<std::string>("separator", ":");
    run_lumi_eventSelector = new RunLumiEventSelector(cfg_runLumiEventSelector);
  }

  std::string selEventsFileName_output = cfg_analyze.getParameter<std::string>("selEventsFileName_output");
  std::cout << "selEventsFileName_output = " << selEventsFileName_output << std::endl;

  fwlite::InputSource inputFiles(cfg);
  int maxEvents = inputFiles.maxEvents();
  std::cout << " maxEvents = " << maxEvents << std::endl;
  int skipSelEvents = cfg_analyze.getParameter<int>("skipSelEvents");
  std::cout << " skipSelEvents = " << skipSelEvents << std::endl;
  int maxSelEvents = cfg_analyze.getParameter<int>("maxSelEvents");
  std::cout << " maxSelEvents = " << maxSelEvents << std::endl;
  unsigned reportEvery = inputFiles.reportAfter();

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());

  TTreeWrapper* inputTree = new TTreeWrapper(treeName.data(), inputFiles.files(), maxEvents);
  std::cout << "Loaded " << inputTree->getFileCount() << " file(s)." << std::endl;

//--- declare event-level variables
  EventInfo eventInfo(analysisConfig);
  EventInfoReader eventInfoReader(&eventInfo);
  inputTree->registerReader(&eventInfoReader);

//--- declare collections of generator-level particles
  GenLeptonReader* genLeptonReader = nullptr;
  if ( branchName_genLeptons != "" ) {
    genLeptonReader = new GenLeptonReader(branchName_genLeptons);
    inputTree->registerReader(genLeptonReader);
  }
  GenParticleReader* genNeutrinoReader = nullptr;
  if ( branchName_genNeutrinos != "" ) {
    genNeutrinoReader = new GenParticleReader(branchName_genNeutrinos);
    inputTree->registerReader(genNeutrinoReader);
  }
  GenJetReader* genJetReader = nullptr;
  if ( branchName_genJets != "" ) {
    genJetReader = new GenJetReader(branchName_genJets);
    inputTree->registerReader(genJetReader);
  }
    
  // collections specific to HH signal
  GenParticleReader* genParticleFromHiggsReader = nullptr;
  if ( branchName_genParticlesFromHiggs != "" ) {
    genParticleFromHiggsReader = new GenParticleReader(branchName_genParticlesFromHiggs);
    inputTree->registerReader(genParticleFromHiggsReader);
  }
  GenParticleReader* genWBosonReader = nullptr;
  if ( branchName_genWBosons != "" ) {
    genWBosonReader = new GenParticleReader(branchName_genWBosons);
    inputTree->registerReader(genWBosonReader);
  }
  GenParticleReader* genWJetReader = nullptr;
  if ( branchName_genWJets != "" ) {
    genWJetReader = new GenParticleReader(branchName_genWJets);
    inputTree->registerReader(genWJetReader);
  }

  // collections specific to ttbar background
  GenLeptonReader* genLeptonFromTopReader = nullptr;  
  if ( branchName_genLeptonsFromTop != "" ) {
    genLeptonFromTopReader = new GenLeptonReader(branchName_genLeptonsFromTop);
    inputTree->registerReader(genLeptonFromTopReader);
  }
  GenParticleReader* genNeutrinoFromTopReader = nullptr;
  if ( branchName_genNeutrinosFromTop != "" ) {
    genNeutrinoFromTopReader = new GenParticleReader(branchName_genNeutrinosFromTop);
    inputTree->registerReader(genNeutrinoFromTopReader);
  }
  GenParticleReader* genBQuarksFromTopReader = nullptr;
  if ( branchName_genBQuarksFromTop != "" ) {
    genBQuarksFromTopReader = new GenParticleReader(branchName_genBQuarksFromTop);
    inputTree->registerReader(genBQuarksFromTopReader);
  }
  GenParticleReader* genWJetsFromTopReader = nullptr;
  if ( branchName_genWJetsFromTop != "" ) {
    genWJetsFromTopReader = new GenParticleReader(branchName_genWJetsFromTop);
    inputTree->registerReader(genWJetsFromTopReader);
  }

  GenLeptonCollectionSelector genLeptonSelector(era, -1, isDEBUG);
  genLeptonSelector.getSelector().set_max_absEta_muon(2.4);
  genLeptonSelector.getSelector().set_max_absEta_electron(2.4);

  GenJetCollectionSelector genJetSelector(era, -1, isDEBUG);
  genJetSelector.getSelector().set_min_pt(20.);
  genJetSelector.getSelector().set_max_absEta(2.4);
  GenJetCollectionCleaner genJetCleaner(0.4, isDEBUG);

//--- declare other generator level information
  LHEInfoReader* lheInfoReader = new LHEInfoReader(hasLHE);
  inputTree->registerReader(lheInfoReader);

//--- open output file containing run:lumi:event numbers of events passing final event selection criteria
  std::ostream* selEventsFile = ( selEventsFileName_output != "" ) ? new std::ofstream(selEventsFileName_output.data(), std::ios::out) : 0;
  std::cout << "selEventsFileName_output = " << selEventsFileName_output << std::endl;

//--- declare histograms
  GenEvtHistManager* genEvtHistManager_beforeCuts = new GenEvtHistManager(makeHistManager_cfg(process_string,
    Form("%s/unbiased/genEvt", histogramDir.data()), era_string, central_or_shift));
  genEvtHistManager_beforeCuts->bookHistograms(fs);
  LHEInfoHistManager* lheInfoHistManager_beforeCuts = new LHEInfoHistManager(makeHistManager_cfg(process_string,
    Form("%s/unbiased/lheInfo", histogramDir.data()), era_string, central_or_shift));
  lheInfoHistManager_beforeCuts->bookHistograms(fs);
  
  struct selHistManagerType
  {
    MEMbbwwHistManagerSingleLepton* mem_2genuineBJets_2genuineWJets_;
    MEMbbwwHistManagerSingleLepton* mem_1genuineBJet_2genuineWJets_;
    MEMbbwwHistManagerSingleLepton* mem_2genuineBJets_1genuineWJet_;
    MEMbbwwHistManagerSingleLepton* mem_1genuineBJet_1genuineWJet_;
    MEMbbwwHistManagerSingleLepton* mem_missingBJet_genuineBJet_2genuineWJets_;
    MEMbbwwHistManagerSingleLepton* mem_missingBJet_fakeBJet_2genuineWJets_;
    MEMbbwwHistManagerSingleLepton* mem_missingWJet_2genuineBJets_genuineWJet_;
    MEMbbwwHistManagerSingleLepton* mem_missingWJet_2genuineBJets_fakeWJet_;
    MEMbbwwHistManagerSingleLepton* mem_missingBnWJet_genuineBJet_genuineWJet_;
    MEMbbwwHistManagerSingleLepton* mem_missingBnWJet_fakeBJet_genuineWJet_;
    MEMbbwwHistManagerSingleLepton* mem_missingBnWJet_genuineBJet_fakeWJet_;
    MEMbbwwHistManagerSingleLepton* mem_missingBnWJet_fakeBJet_fakeWJet_;
    GenEvtHistManager* genEvtHistManager_afterCuts_;
    LHEInfoHistManager* lheInfoHistManager_afterCuts_;
    WeightHistManager* weights_;
  };
  selHistManagerType* selHistManager = new selHistManagerType();
  selHistManager->mem_2genuineBJets_2genuineWJets_ = new MEMbbwwHistManagerSingleLepton(makeHistManager_cfg(process_string,
    Form("%s/sel/mem_2genuineBJets_2genuineWJets", histogramDir.data()), era_string, central_or_shift));
  selHistManager->mem_2genuineBJets_2genuineWJets_->bookHistograms(fs);
  selHistManager->mem_1genuineBJet_2genuineWJets_ = new MEMbbwwHistManagerSingleLepton(makeHistManager_cfg(process_string,
    Form("%s/sel/mem_1genuineBJet_2genuineWJets", histogramDir.data()), era_string, central_or_shift));
  selHistManager->mem_1genuineBJet_2genuineWJets_->bookHistograms(fs);
  selHistManager->mem_2genuineBJets_1genuineWJet_ = new MEMbbwwHistManagerSingleLepton(makeHistManager_cfg(process_string,
    Form("%s/sel/mem_2genuineBJets_1genuineWJet", histogramDir.data()), era_string, central_or_shift));
  selHistManager->mem_2genuineBJets_1genuineWJet_->bookHistograms(fs);
  selHistManager->mem_1genuineBJet_1genuineWJet_ = new MEMbbwwHistManagerSingleLepton(makeHistManager_cfg(process_string,
    Form("%s/sel/mem_1genuineBJet_1genuineWJet", histogramDir.data()), era_string, central_or_shift));
  selHistManager->mem_1genuineBJet_1genuineWJet_->bookHistograms(fs);
  selHistManager->mem_missingBJet_genuineBJet_2genuineWJets_ = new MEMbbwwHistManagerSingleLepton(makeHistManager_cfg(process_string,
    Form("%s/sel/mem_missingBJet_genuineBJet_2genuineWJets", histogramDir.data()), era_string, central_or_shift));
  selHistManager->mem_missingBJet_genuineBJet_2genuineWJets_->bookHistograms(fs);
  selHistManager->mem_missingBJet_fakeBJet_2genuineWJets_ = new MEMbbwwHistManagerSingleLepton(makeHistManager_cfg(process_string,
    Form("%s/sel/mem_missingBJet_fakeBJet_2genuineWJets", histogramDir.data()), era_string, central_or_shift));
  selHistManager->mem_missingBJet_fakeBJet_2genuineWJets_->bookHistograms(fs);
  selHistManager->mem_missingWJet_2genuineBJets_genuineWJet_ = new MEMbbwwHistManagerSingleLepton(makeHistManager_cfg(process_string,
    Form("%s/sel/mem_missingWJet_2genuineBJets_genuineWJet", histogramDir.data()), era_string, central_or_shift));
  selHistManager->mem_missingWJet_2genuineBJets_genuineWJet_->bookHistograms(fs);
  selHistManager->mem_missingWJet_2genuineBJets_fakeWJet_ = new MEMbbwwHistManagerSingleLepton(makeHistManager_cfg(process_string,
    Form("%s/sel/mem_missingWJet_2genuineBJets_fakeWJet", histogramDir.data()), era_string, central_or_shift));
  selHistManager->mem_missingWJet_2genuineBJets_fakeWJet_->bookHistograms(fs);
  selHistManager->mem_missingBnWJet_genuineBJet_genuineWJet_ = new MEMbbwwHistManagerSingleLepton(makeHistManager_cfg(process_string,
    Form("%s/sel/mem_missingBnWJet_genuineBJet_genuineWJet", histogramDir.data()), era_string, central_or_shift));
  selHistManager->mem_missingBnWJet_genuineBJet_genuineWJet_->bookHistograms(fs);
  selHistManager->mem_missingBnWJet_fakeBJet_genuineWJet_ = new MEMbbwwHistManagerSingleLepton(makeHistManager_cfg(process_string,
    Form("%s/sel/mem_missingBnWJet_fakeBJet_genuineWJet", histogramDir.data()), era_string, central_or_shift));
  selHistManager->mem_missingBnWJet_fakeBJet_genuineWJet_->bookHistograms(fs);
  selHistManager->mem_missingBnWJet_genuineBJet_fakeWJet_ = new MEMbbwwHistManagerSingleLepton(makeHistManager_cfg(process_string,
    Form("%s/sel/mem_missingBnWJet_genuineBJet_fakeWJet", histogramDir.data()), era_string, central_or_shift));
  selHistManager->mem_missingBnWJet_genuineBJet_fakeWJet_->bookHistograms(fs);
  selHistManager->mem_missingBnWJet_fakeBJet_fakeWJet_ = new MEMbbwwHistManagerSingleLepton(makeHistManager_cfg(process_string,
    Form("%s/sel/mem_missingBnWJet_fakeBJet_fakeWJet", histogramDir.data()), era_string, central_or_shift));
  selHistManager->mem_missingBnWJet_fakeBJet_fakeWJet_->bookHistograms(fs);
  selHistManager->genEvtHistManager_afterCuts_ = new GenEvtHistManager(makeHistManager_cfg(process_string,
    Form("%s/sel/genEvt", histogramDir.data()), era_string, central_or_shift));
  selHistManager->genEvtHistManager_afterCuts_->bookHistograms(fs);
  selHistManager->lheInfoHistManager_afterCuts_ = new LHEInfoHistManager(makeHistManager_cfg(process_string,
    Form("%s/sel/lheInfo", histogramDir.data()), era_string, central_or_shift));
  selHistManager->lheInfoHistManager_afterCuts_->bookHistograms(fs);
  selHistManager->weights_ = new WeightHistManager(makeHistManager_cfg(process_string,
    Form("%s/sel/weights", histogramDir.data()), era_string, central_or_shift));
  selHistManager->weights_->bookHistograms(fs, { "genWeight", "pileupWeight" });

  std::string ntupleDir = Form("%s/ntuples/%s", histogramDir.data(), process_string.data());
  MEMbbwwNtupleManager_singlelepton* mem_ntuple = new MEMbbwwNtupleManager_singlelepton(ntupleDir, "mem");
  mem_ntuple->makeTree(fs);
  mem_ntuple->initializeBranches();
  MEMbbwwNtupleManager_singlelepton* mem_ntuple_missingBJet = new MEMbbwwNtupleManager_singlelepton(ntupleDir, "mem_missingBJet");
  mem_ntuple_missingBJet->makeTree(fs);
  mem_ntuple_missingBJet->initializeBranches();
  MEMbbwwNtupleManager_singlelepton* mem_ntuple_missingWJet = new MEMbbwwNtupleManager_singlelepton(ntupleDir, "mem_missingWJet");
  mem_ntuple_missingWJet->makeTree(fs);
  mem_ntuple_missingWJet->initializeBranches();
  MEMbbwwNtupleManager_singlelepton* mem_ntuple_missingBnWJet = new MEMbbwwNtupleManager_singlelepton(ntupleDir, "mem_missingBnWJet");
  mem_ntuple_missingBnWJet->makeTree(fs);
  mem_ntuple_missingBnWJet->initializeBranches();

  int analyzedEntries = 0;
  int skippedEntries = 0;
  int selectedEntries = 0;
  double selectedEntries_weighted = 0.;
  TH1* histogram_analyzedEntries = fs.make<TH1D>("analyzedEntries", "analyzedEntries", 1, -0.5, +0.5);
  TH1* histogram_selectedEntries = fs.make<TH1D>("selectedEntries", "selectedEntries", 1, -0.5, +0.5);
  cutFlowTableType cutFlowTable;
  const edm::ParameterSet cutFlowTableCfg = makeHistManager_cfg(
    process_string, Form("%s/sel/cutFlow", histogramDir.data()), era_string, central_or_shift
  );
  const std::vector<std::string> cuts = {
    "run:ls:event selection",
    ">= 1 gen lepton",
    //"gen electron (muon) pT > 32 (25) GeV",
    "gen lepton pT > 25 GeV",
    ">= 2 gen b-jets",
    ">= 2 gen jets from W->jj"
  };
  CutFlowTableHistManager * cutFlowHistManager = new CutFlowTableHistManager(cutFlowTableCfg, cuts);
  cutFlowHistManager->bookHistograms(fs);
  while ( inputTree->hasNextEvent() && (! run_lumi_eventSelector || (run_lumi_eventSelector && ! run_lumi_eventSelector -> areWeDone())) && selectedEntries < maxSelEvents ) {
    if ( inputTree -> canReport(reportEvery) ) {
      std::cout << "processing Entry " << inputTree -> getCurrentMaxEventIdx()
                << " or " << inputTree -> getCurrentEventIdx() << " entry in #"
                << (inputTree -> getProcessedFileCount() - 1)
                << " (" << eventInfo
                << ") file (" << selectedEntries << " Entries selected)\n";
    }
    ++analyzedEntries;
    histogram_analyzedEntries->Fill(0.);

    if ( isDEBUG ) {
      std::cout << "event #" << inputTree -> getCurrentMaxEventIdx() << ' ' << eventInfo << '\n';
    }

    double evtWeight = 1.;
    if ( apply_genWeight ) evtWeight *= boost::math::sign(eventInfo.genWeight);
    lheInfoReader->read();
    evtWeight *= lheInfoReader->getWeight_scale(kLHE_scale_central);
    evtWeight *= eventInfo.pileupWeight;
    
    if ( run_lumi_eventSelector && !(*run_lumi_eventSelector)(eventInfo) ) continue;
    cutFlowTable.update("run:ls:event selection");
    cutFlowHistManager->fillHistograms("run:ls:event selection", evtWeight);

    if ( run_lumi_eventSelector ) {
      std::cout << "processing Entry #" << inputTree->getCumulativeMaxEventCount() << ": " << eventInfo << std::endl;
      if ( inputTree -> isOpen() ) {
        std::cout << "input File = " << inputTree->getCurrentFileName() << std::endl;
      }
    }

//--- build collections of generator level particles (before any cuts are applied, to check distributions in unbiased event samples)
    std::vector<GenLepton> genLeptons;
    std::vector<GenLepton> genElectrons;
    std::vector<GenLepton> genMuons;
    if ( genLeptonReader ) {
      genLeptons = genLeptonReader->read();
      for ( std::vector<GenLepton>::const_iterator genLepton = genLeptons.begin();
	    genLepton != genLeptons.end(); ++genLepton ) {
	int abs_pdgId = std::abs(genLepton->pdgId());
	if      ( abs_pdgId == 11 ) genElectrons.push_back(*genLepton);
	else if ( abs_pdgId == 13 ) genMuons.push_back(*genLepton);
      }
    }
    std::vector<GenParticle> genNeutrinos;
    if ( genNeutrinoReader ) {
      genNeutrinos = genNeutrinoReader->read();
    }
    std::vector<GenJet> genJets;
    if ( genJetReader ) {
      genJets = genJetReader->read();
    }
    if ( isDEBUG ) {
      printCollection("genLeptons", genLeptons);
      printCollection("genNeutrinos", genNeutrinos);
      printCollection("genJets", genJets);
    }

    genEvtHistManager_beforeCuts->fillHistograms(genElectrons, genMuons, {}, {}, genJets, evtWeight);

    std::vector<GenParticle> genParticlesFromHiggs;
    std::vector<GenParticle> genWBosons;
    std::vector<GenParticle> genWJets;
    if ( isSignal ) {
      genParticlesFromHiggs = genParticleFromHiggsReader->read();
      genWBosons = genWBosonReader->read();
      genWJets = genWJetReader->read();
      if ( isDEBUG ) {
	printCollection("genParticlesFromHiggs", genParticlesFromHiggs);
        printCollection("genWBosons", genWBosons);
        printCollection("genWJets", genWJets);
      }
      if ( !(genParticlesFromHiggs.size() == 4 && genWBosons.size() == 2 && genWJets.size() == 2) ) {
	if ( run_lumi_eventSelector ) {
	  std::cout << "event " << eventInfo.str() << " FAILS generator-level selection." << std::endl;
	  std::cout << "#genParticlesFromHiggs = " << genParticlesFromHiggs.size() << std::endl;
          std::cout << "#genWBosons = " << genWBosons.size() << std::endl;
          std::cout << "#genWJets = " << genWJets.size() << std::endl;
	}
	continue;
      }
    } 
    std::vector<GenLepton> genLeptonsFromTop;
    std::vector<GenParticle> genNeutrinosFromTop;
    std::vector<GenParticle> genBQuarksFromTop;
    std::vector<GenParticle> genWJetsFromTop;
    if ( !isSignal ) {
      genLeptonsFromTop = genLeptonFromTopReader->read();
      genNeutrinosFromTop = genNeutrinoFromTopReader->read();
      genBQuarksFromTop = genBQuarksFromTopReader->read();
      genWJetsFromTop = genWJetsFromTopReader->read();
      if ( isDEBUG ) {
	printCollection("genLeptonsFromTop", genLeptonsFromTop);
	printCollection("genNeutrinosFromTop", genNeutrinosFromTop);
	printCollection("genBQuarksFromTop", genBQuarksFromTop);
        printCollection("genWJetsFromTop", genWJetsFromTop);
      }
      if ( !(genLeptonsFromTop.size() == 1 && genNeutrinosFromTop.size() == 1 && genBQuarksFromTop.size() == 2 && genWJetsFromTop.size() == 2) ) {
	if ( run_lumi_eventSelector ) {
	  std::cout << "event " << eventInfo.str() << " FAILS generator-level selection." << std::endl;
	  std::cout << "#genLeptonsFromTop = " << genLeptonsFromTop.size() << std::endl;
	  std::cout << "#genNeutrinosFromTop = " << genNeutrinosFromTop.size() << std::endl;
	  std::cout << "#genBQuarksFromTop = " << genBQuarksFromTop.size() << std::endl;
          std::cout << "#genWJetsFromTop = " << genWJetsFromTop.size() << std::endl;
	}
	continue;
      }
    }
    cutFlowTable.update("generator-level selection (1)", evtWeight);

//--- select lepton, light-quark jets, and b-jets from from H->WW->lnuqq and H->bb decays (signal)
//    and from tt->bWbW->blnu bqq decays (background)
    std::vector<GenLepton> genLeptonsForMatching;
    std::vector<GenJet> genWJetsForMatching;
    std::vector<GenJet> genBJetsForMatching;
    double genMEtPx = 0.;
    double genMEtPy = 0.;
    if ( isSignal ) {
      const GenParticle* genBQuark      = nullptr;
      const GenParticle* genAntiBQuark  = nullptr;
      const GenParticle* genWBosonPlus  = nullptr;
      const GenParticle* genWBosonMinus = nullptr;
      for ( std::vector<GenParticle>::const_iterator genParticle = genParticlesFromHiggs.begin();
	    genParticle != genParticlesFromHiggs.end(); ++genParticle ) {
	if      ( genParticle->pdgId() ==  +5 ) genBQuark      = &(*genParticle);
	else if ( genParticle->pdgId() ==  -5 ) genAntiBQuark  = &(*genParticle);
	else if ( genParticle->pdgId() == +24 ) genWBosonPlus  = &(*genParticle);
	else if ( genParticle->pdgId() == -24 ) genWBosonMinus = &(*genParticle);
      }
      if ( genBQuark && genAntiBQuark ) {
	genBJetsForMatching.push_back(GenJet(
          genBQuark->pt(), genBQuark->eta(), genBQuark->phi(), mem::bottomQuarkMass, genBQuark->pdgId()));
        genBJetsForMatching.push_back(GenJet(
          genAntiBQuark->pt(), genAntiBQuark->eta(), genAntiBQuark->phi(), mem::bottomQuarkMass, genAntiBQuark->pdgId()));	
      }
      if ( genWBosonPlus && genWBosonMinus ) {
        const GenParticle* genHadWBoson = nullptr;
	std::pair<const GenLepton*, const GenParticle*> genLepton_and_NeutrinoFromWBosonPlus =
          findGenLepton_and_NeutrinoFromWBoson(*genWBosonPlus, genLeptons, genNeutrinos);
        if ( genLepton_and_NeutrinoFromWBosonPlus.first && genLepton_and_NeutrinoFromWBosonPlus.second ) {
          genLeptonsForMatching.push_back(*genLepton_and_NeutrinoFromWBosonPlus.first);
          genMEtPx += genLepton_and_NeutrinoFromWBosonPlus.second->p4().px();
          genMEtPy += genLepton_and_NeutrinoFromWBosonPlus.second->p4().py();
          genHadWBoson = genWBosonMinus;
        }
	std::pair<const GenLepton*, const GenParticle*> genLepton_and_NeutrinoFromWBosonMinus =
          findGenLepton_and_NeutrinoFromWBoson(*genWBosonMinus, genLeptons, genNeutrinos);
        if ( genLepton_and_NeutrinoFromWBosonMinus.first && genLepton_and_NeutrinoFromWBosonMinus.second ) {
          genLeptonsForMatching.push_back(*genLepton_and_NeutrinoFromWBosonMinus.first);
          genMEtPx += genLepton_and_NeutrinoFromWBosonMinus.second->p4().px();
          genMEtPy += genLepton_and_NeutrinoFromWBosonMinus.second->p4().py();
          genHadWBoson = genWBosonPlus;
        }
        // CV: skip events for which matching of generator-level lepton+neutrino to W boson is ambiguous
        if ( genLeptonsForMatching.size() != 1 ) continue;
        std::vector<GenJet> genWJets_tmp;
        for ( std::vector<GenParticle>::const_iterator genWJet = genWJets.begin();
              genWJet != genWJets.end(); ++genWJet ) {
          genWJets_tmp.push_back(GenJet(
            genWJet->pt(), genWJet->eta(), genWJet->phi(), genWJet->mass(), genWJet->pdgId()));
        }
        assert(genHadWBoson);
        std::vector<const GenJet*> genWJetsForMatching_tmp = findGenJetsFromWBoson(*genHadWBoson, genWJets_tmp);        
        std::sort(genWJetsForMatching_tmp.begin(), genWJetsForMatching_tmp.end(), isHigherPt);
        for ( std::vector<const GenJet*>::const_iterator genWJet = genWJetsForMatching_tmp.begin();
              genWJet != genWJetsForMatching_tmp.end(); ++genWJet ) {
          genWJetsForMatching.push_back(**genWJet);
        }
      }
    } else {
      genLeptonsForMatching = genLeptonsFromTop;
      assert(genNeutrinosFromTop.size() == 1);
      genMEtPx = genNeutrinosFromTop[0].p4().px();
      genMEtPy = genNeutrinosFromTop[0].p4().py();
      for ( std::vector<GenParticle>::const_iterator genBQuark = genBQuarksFromTop.begin();
	    genBQuark != genBQuarksFromTop.end(); ++genBQuark ) {
	genBJetsForMatching.push_back(GenJet(
          genBQuark->pt(), genBQuark->eta(), genBQuark->phi(), mem::bottomQuarkMass, genBQuark->pdgId()));
      }
      for ( std::vector<GenParticle>::const_iterator genWJet = genWJetsFromTop.begin();
            genWJet != genWJetsFromTop.end(); ++genWJet ) {
        genWJetsForMatching.push_back(GenJet(
          genWJet->pt(), genWJet->eta(), genWJet->phi(), genWJet->mass(), genWJet->pdgId()));
      }
      std::sort(genWJetsForMatching.begin(), genWJetsForMatching.end(), isHigherPtT<GenJet>);
    }
    if ( !(genLeptonsForMatching.size() == 1 && genWJetsForMatching.size() == 2 && genBJetsForMatching.size() == 2) ) {
      if ( run_lumi_eventSelector ) {
	std::cout << "event " << eventInfo.str() << " FAILS generator-level selection." << std::endl;
	std::cout << "#genLeptonsForMatching = " << genLeptonsForMatching.size() << std::endl;
        std::cout << "#genWJetsForMatching = " << genWJetsForMatching.size() << std::endl;
	std::cout << "#genBJetsForMatching = " << genBJetsForMatching.size() << std::endl;
      }
      continue;
    }
    cutFlowTable.update("generator-level selection (2)", evtWeight);

//--- apply pT and eta cuts to generator-level lepton, light-quark jets, and b-jets,
//    clean collection of generator-level b-jets with respect to leptons,
//    and collection of generator-level light-quark jets with respect to leptons and b-jets
    std::vector<const GenLepton*> genLeptonsForMatching_ptrs = convert_to_ptrs(genLeptonsForMatching);
    std::vector<const GenLepton*> selGenLeptons = genLeptonSelector(genLeptonsForMatching_ptrs, isHigherPt);

    std::vector<const GenJet*> genBJetsForMatching_ptrs = convert_to_ptrs(genBJetsForMatching);
    std::vector<const GenJet*> cleanedGenBJets = genJetCleaner(genBJetsForMatching_ptrs, genLeptonsForMatching_ptrs);
    std::vector<const GenJet*> selGenBJets = genJetSelector(cleanedGenBJets, isHigherPt);

    std::vector<const GenJet*> genWJetsForMatching_ptrs = convert_to_ptrs(genWJetsForMatching);
    std::vector<const GenJet*> cleanedGenWJets = genJetCleaner(genWJetsForMatching_ptrs, genLeptonsForMatching_ptrs, genBJetsForMatching_ptrs);
    std::vector<const GenJet*> selGenWJets = genJetSelector(cleanedGenWJets, isHigherPt);
    std::set<size_t> usedGenWJets;

    std::vector<const GenJet*> genJets_ptrs = convert_to_ptrs(genJets);
    std::vector<const GenJet*> cleanedGenJets = genJetCleaner(genJets_ptrs, genLeptonsForMatching_ptrs, genBJetsForMatching_ptrs, genWJetsForMatching_ptrs);
    std::vector<const GenJet*> selGenJets = genJetSelector(cleanedGenJets, isHigherPt);
    std::set<size_t> usedGenJets;

    GenMEt genMEt(genMEtPx, genMEtPy);

    //std::cout << "#selGenBJets = " << selGenBJets.size() << std::endl;
    //std::cout << "#selGenWJets = " << selGenWJets.size() << std::endl;
    //std::cout << "#selGenJets = " << selGenJets.size() << std::endl;
    
//--- apply pT smearing to generator-level b-jets (and other jets)
    const double genBJet_pFake = 0.10;
    std::vector<GenJet> selGenBJets_smeared;
    bool selGenBJet_lead_isFake = false;
    bool selGenBJet_sublead_isFake = false;
    for ( size_t idxGenBJet = 0; idxGenBJet < selGenBJets.size(); ++idxGenBJet ) {
      const GenJet* selGenBJet = selGenBJets[idxGenBJet];
      const GenJet* genJet = nullptr;
      bool genJet_isFake;
      double u = rnd.Uniform();
      assert(u >= 0. && u <= 1.);
      if ( u > genBJet_pFake ) {
	genJet = selGenBJet;
        genJet_isFake = false;
      } else if ( (selGenWJets.size() + selGenJets.size()) > (usedGenWJets.size() + usedGenJets.size()) ) {
        int idxGenWJet = -1;
        int idxGenJet = -1;
        while ( idxGenWJet == -1 && idxGenJet == -1 ) {
          int idxGenJet_tmp = TMath::Nint(rnd.Uniform(-0.5, selGenWJets.size() + selGenJets.size() - 0.5));
          if ( idxGenJet_tmp < (int)selGenWJets.size() ) {
            if ( usedGenWJets.find(idxGenJet_tmp) == usedGenWJets.end() ) {
              idxGenWJet = idxGenJet_tmp;
              usedGenWJets.insert(idxGenWJet); 
            }
          } else {
            idxGenJet_tmp -= selGenWJets.size();
            if ( usedGenJets.find(idxGenJet_tmp) == usedGenJets.end() ) {
              idxGenJet = idxGenJet_tmp;
              usedGenJets.insert(idxGenJet); 
            }
          }
        }
        if ( idxGenWJet >= 0 && idxGenWJet < (int)selGenWJets.size() ) {
          genJet = selGenWJets[idxGenWJet];
	  genJet_isFake = true;
        } else if ( idxGenJet >= 0 && idxGenJet < (int)selGenJets.size() ) {
          genJet = selGenJets[idxGenJet];
	  genJet_isFake = true;
        } else assert(0);
      }
      if ( genJet ) {
	double genJetPt_smeared;
	if ( apply_jetSmearing ) genJetPt_smeared = genJetSmearer(*genJet).pt();
	else genJetPt_smeared = genJet->pt();
	if ( genJetPt_smeared > genJetSelector.getSelector().get_min_pt() ) {
  	  selGenBJets_smeared.push_back(GenJet(
	    genJetPt_smeared, genJet->eta(), genJet->phi(), genJet->mass(), genJet->pdgId()));
	}
	if      ( idxGenBJet == 0 ) selGenBJet_lead_isFake = genJet_isFake;
	else if ( idxGenBJet == 1 ) selGenBJet_sublead_isFake = genJet_isFake;
	else assert(0);
      }
    }

//--- apply pT smearing to generator-level light-quark jets (and other jets)
    const double genWJet_lead_pFake = 0.10;
    const double genWJet_sublead_pFake = 0.30;
    std::vector<GenJet> selGenWJets_smeared;
    bool selGenWJet_lead_isFake = false;
    bool selGenWJet_sublead_isFake = false;
    for ( size_t idxGenWJet = 0; idxGenWJet < selGenWJets.size(); ++idxGenWJet ) {
      const GenJet* selGenWJet = selGenWJets[idxGenWJet];
      const GenJet* genJet = nullptr;
      bool genJet_isFake;
      double u = rnd.Uniform();
      assert(u >= 0. && u <= 1.);
      double genWJet_pFake = ( idxGenWJet == 0 ) ? genWJet_lead_pFake : genWJet_sublead_pFake;
      if ( u > genWJet_pFake && usedGenWJets.find(idxGenWJet) == usedGenWJets.end() ) {
	genJet = selGenWJet;
        genJet_isFake = false;
      } else if ( selGenJets.size() > usedGenJets.size() ) {
        int idxGenJet = -1;
        while ( idxGenJet == -1 ) {
          int idxGenJet_tmp = TMath::Nint(rnd.Uniform(-0.5, selGenJets.size() - 0.5));
          if ( usedGenJets.find(idxGenJet_tmp) == usedGenJets.end() ) {
            idxGenJet = idxGenJet_tmp;
            usedGenJets.insert(idxGenJet); 
          }
        }
        assert(idxGenJet >= 0 && idxGenJet < (int)selGenJets.size());
        genJet = selGenJets[idxGenJet];
	genJet_isFake = true;
      }
      if ( genJet ) {
	double genJetPt_smeared;
	if ( apply_jetSmearing ) genJetPt_smeared = genJetSmearer(*genJet).pt();
	else genJetPt_smeared = genJet->pt();
	if ( genJetPt_smeared > genJetSelector.getSelector().get_min_pt() ) {
  	  selGenWJets_smeared.push_back(GenJet(
	    genJetPt_smeared, genJet->eta(), genJet->phi(), genJet->mass(), genJet->pdgId()));
	}
	if      ( idxGenWJet == 0 ) selGenWJet_lead_isFake = genJet_isFake;
	else if ( idxGenWJet == 1 ) selGenWJet_sublead_isFake = genJet_isFake;
	else assert(0);
      }
    }    
    
//--- apply pX, pY smearing to generator-level missing transverse momentum (MET)
    GenMEt genMEt_smeared;
    if ( apply_metSmearing ) {
      genMEt_smeared = genMEtSmearer(genMEt);
    } else {
      genMEt_smeared = genMEt;
    }

    // require one or more generator-level leptons 
    if ( !(selGenLeptons.size() >= 1) ) {
      if ( run_lumi_eventSelector ) {
        std::cout << "event " << eventInfo.str() << " FAILS selGenLeptons selection." << std::endl;
        printCollection("selGenLeptons", selGenLeptons);
      }
      continue;
    }
    cutFlowTable.update(">= 1 gen lepton", evtWeight);
    cutFlowHistManager->fillHistograms(">= 1 gen lepton", evtWeight);
    const GenLepton* selGenLepton = selGenLeptons[0];
    
    //const double minPt_lepton = ( TMath::Abs(selGenLepton->pdgId()) == 11 ) ? 32. : 25.;
    const double minPt_lepton = 25.;
    const double maxAbsEta_lepton = 2.4;
    if ( !(selGenLepton->pt() > minPt_lepton && selGenLepton->absEta() < maxAbsEta_lepton) ) {
      if ( run_lumi_eventSelector ) {
        std::cout << "event " << eventInfo.str() << " FAILS lepton pT selection." << std::endl;
        std::cout << " selGenLepton: pT = " << selGenLepton->pt() << ", eta = " << selGenLepton->eta() << " "
		  << "(minPt = " << minPt_lepton <<  ", maxAbsEta = " << maxAbsEta_lepton << ")" << std::endl;
      }
      continue;
    }
    //cutFlowTable.update("gen electron (muon) pT > 32 (25) GeV", evtWeight);
    //cutFlowHistManager->fillHistograms("gen electron (muon) pT > 32 (25) GeV", evtWeight);
    cutFlowTable.update("gen lepton pT > 25 GeV", evtWeight);
    cutFlowHistManager->fillHistograms("gen lepton pT > 25 GeV", evtWeight);

    if ( !(selGenBJets_smeared.size() == 2) ) {
      if ( run_lumi_eventSelector ) {
	std::cout << "event " << eventInfo.str() << " FAILS gen smeared b-jets selection." << std::endl;
	std::cout << "#cleanedGenBJets = " << cleanedGenBJets.size() << std::endl;
	std::cout << "#selGenBJets = " << selGenBJets.size() << std::endl;
        std::cout << "#cleanedGenWJets = " << cleanedGenWJets.size() << std::endl;
	std::cout << "#selGenWJets = " << selGenWJets.size() << std::endl;
	std::cout << "#cleanedGenJets = " << cleanedGenJets.size() << std::endl;
	std::cout << "#selGenJets = " << selGenJets.size() << std::endl;
	std::cout << "#selGenBJets_smeared = " << selGenBJets_smeared.size() << std::endl;
      }
      continue;
    }
    cutFlowTable.update(">= 2 gen b-jets", evtWeight);
    cutFlowHistManager->fillHistograms(">= 2 gen b-jets", evtWeight);
    const GenJet* selGenBJet_lead = &selGenBJets_smeared[0];
    const GenJet* selGenBJet_sublead = &selGenBJets_smeared[1];

    if ( !(selGenWJets_smeared.size() == 2) ) {
      if ( run_lumi_eventSelector ) {
	std::cout << "event " << eventInfo.str() << " FAILS gen smeared jets from W->jj selection." << std::endl;
        std::cout << "#cleanedGenWJets = " << cleanedGenWJets.size() << std::endl;
	std::cout << "#selGenWJets = " << selGenWJets.size() << std::endl;
	std::cout << "#cleanedGenJets = " << cleanedGenJets.size() << std::endl;        
	std::cout << "#selGenJets = " << selGenJets.size() << std::endl;
	std::cout << "#selGenWJets_smeared = " << selGenWJets_smeared.size() << std::endl;
      }
      continue;
    }
    cutFlowTable.update(">= 2 gen jets from W->jj", evtWeight);
    cutFlowHistManager->fillHistograms(">= 2 gen jets from W->jj", evtWeight);
    const GenJet* selGenWJet_lead = &selGenWJets_smeared[0];
    const GenJet* selGenWJet_sublead = &selGenWJets_smeared[1];
    
    //---------------------------------------------------------------------------
    // CV: Skip running matrix element method (MEM) computation for the first 'skipSelEvents' events.
    //     This feature allows to process the HH signal samples in chunks of 'maxSelEvents' events per job.
    ++skippedEntries;
    if ( skippedEntries < skipSelEvents ) continue;
    //---------------------------------------------------------------------------

    //---------------------------------------------------------------------------
    // CV: Compute MEM likelihood ratio of HH signal and ttbar background hypotheses

    if ( isDEBUG ) {
      std::cout << "selGenLepton: pT = " << selGenLepton->pt() << "," 
                << " eta = " << selGenLepton->eta() << ", phi = " << selGenLepton->phi() << std::endl;
      std::cout << "selGenBJet_lead: pT = " << selGenBJet_lead->pt() << "," 
                << " eta = " << selGenBJet_lead->eta() << ", phi = " << selGenBJet_lead->phi() 
	        << " (isFake = " << selGenBJet_lead_isFake << ")" << std::endl;
      std::cout << "selGenBJet_sublead: pT = " << selGenBJet_sublead->pt() << "," 
                << " eta = " << selGenBJet_sublead->eta() << ", phi = " << selGenBJet_sublead->phi() 
	        << " (isFake = " << selGenBJet_sublead_isFake << ")" << std::endl;
      std::cout << "selGenWJet_lead: pT = " << selGenWJet_lead->pt() << "," 
                << " eta = " << selGenWJet_lead->eta() << ", phi = " << selGenWJet_lead->phi() 
	        << " (isFake = " << selGenWJet_lead_isFake << ")" << std::endl;
      std::cout << "selGenWJet_sublead: pT = " << selGenWJet_sublead->pt() << "," 
                << " eta = " << selGenWJet_sublead->eta() << ", phi = " << selGenWJet_sublead->phi() 
	        << " (isFake = " << selGenWJet_sublead_isFake << ")" << std::endl;
    }

    int memLeptonType;
    double memLeptonMass;
    if ( selGenLepton->is_electron() ) {
      memLeptonType = mem::MeasuredParticle::kElectron;
      memLeptonMass = mem::electronMass;
    } else if ( selGenLepton->is_muon() ) {
      memLeptonType = mem::MeasuredParticle::kMuon;
      memLeptonMass = mem::muonMass;
    } else assert(0);

    mem::MeasuredParticle memMeasuredLepton(memLeptonType, 
      selGenLepton->pt(), selGenLepton->eta(), selGenLepton->phi(), 
      memLeptonMass, selGenLepton->charge());
    mem::MeasuredParticle memMeasuredBJet_lead(mem::MeasuredParticle::kBJet,
      selGenBJet_lead->pt(), selGenBJet_lead->eta(), selGenBJet_lead->phi(), 
      mem::bottomQuarkMass);
    mem::MeasuredParticle memMeasuredBJet_sublead(mem::MeasuredParticle::kBJet,
      selGenBJet_sublead->pt(), selGenBJet_sublead->eta(), selGenBJet_sublead->phi(), 
      mem::bottomQuarkMass);
    mem::MeasuredParticle memMeasuredWJet_lead(mem::MeasuredParticle::kHadWJet,
      selGenWJet_lead->pt(), selGenWJet_lead->eta(), selGenWJet_lead->phi(), 
      selGenWJet_lead->mass());
    mem::MeasuredParticle memMeasuredWJet_sublead(mem::MeasuredParticle::kHadWJet,
      selGenWJet_sublead->pt(), selGenWJet_sublead->eta(), selGenWJet_sublead->phi(), 
      selGenWJet_sublead->mass());

    std::vector<mem::MeasuredParticle> memMeasuredParticles;
    memMeasuredParticles.push_back(memMeasuredLepton);
    memMeasuredParticles.push_back(memMeasuredBJet_lead);
    memMeasuredParticles.push_back(memMeasuredBJet_sublead);
    memMeasuredParticles.push_back(memMeasuredWJet_lead);
    memMeasuredParticles.push_back(memMeasuredWJet_sublead);
    
    MEMEvent_singlelepton memEvent(
      eventInfo, isSignal, 
      &memMeasuredBJet_lead, &memMeasuredBJet_sublead,
      &memMeasuredWJet_lead, &memMeasuredWJet_sublead,
      &memMeasuredLepton,
      genMEt_smeared.px(), genMEt_smeared.py(), metCov);
    addGenMatches_singlelepton(memEvent, genBJetsForMatching_ptrs, genWJetsForMatching_ptrs, genLeptonsForMatching_ptrs, genMEtPx, genMEtPy);

    std::vector<mem::MeasuredParticle> memMeasuredParticles_missingBJet;
    memMeasuredParticles_missingBJet.push_back(memMeasuredLepton);
    const mem::MeasuredParticle* memMeasuredBJet_missingBJet = nullptr;
    bool selGenBJet_isFake_missingBJet;
    double u1 = rnd.Uniform();
    assert(u1 >= 0. && u1 <= 1.);
    if ( u1 > 0.50 ) {
      memMeasuredParticles_missingBJet.push_back(memMeasuredBJet_lead);
      memMeasuredBJet_missingBJet = &memMeasuredBJet_lead;
      selGenBJet_isFake_missingBJet = selGenBJet_lead_isFake;
    } else {
      memMeasuredParticles_missingBJet.push_back(memMeasuredBJet_sublead);
      memMeasuredBJet_missingBJet = &memMeasuredBJet_sublead;
      selGenBJet_isFake_missingBJet = selGenBJet_sublead_isFake;
    }
    memMeasuredParticles_missingBJet.push_back(memMeasuredWJet_lead);
    memMeasuredParticles_missingBJet.push_back(memMeasuredWJet_sublead);

    MEMEvent_singlelepton memEvent_missingBJet(
      eventInfo, isSignal, 
      memMeasuredBJet_missingBJet, nullptr,
      &memMeasuredWJet_lead, &memMeasuredWJet_sublead,
      &memMeasuredLepton,
      genMEt_smeared.px(), genMEt_smeared.py(), metCov);
    addGenMatches_singlelepton(memEvent_missingBJet, genBJetsForMatching_ptrs, genWJetsForMatching_ptrs, genLeptonsForMatching_ptrs, genMEtPx, genMEtPy);

    std::vector<mem::MeasuredParticle> memMeasuredParticles_missingWJet;
    memMeasuredParticles_missingWJet.push_back(memMeasuredLepton);
    memMeasuredParticles_missingWJet.push_back(memMeasuredBJet_lead);
    memMeasuredParticles_missingWJet.push_back(memMeasuredBJet_sublead);
    const mem::MeasuredParticle* memMeasuredWJet_missingWJet = nullptr;
    bool selGenWJet_isFake_missingWJet;
    double u2 = rnd.Uniform();
    assert(u2 >= 0. && u2 <= 1.);
    if ( u2 > 0.50 ) {
      memMeasuredParticles_missingWJet.push_back(memMeasuredWJet_lead);
      memMeasuredWJet_missingWJet = &memMeasuredWJet_lead;
      selGenWJet_isFake_missingWJet = selGenWJet_lead_isFake;
    } else {
      memMeasuredParticles_missingWJet.push_back(memMeasuredWJet_sublead);
      memMeasuredWJet_missingWJet = &memMeasuredWJet_sublead;
      selGenWJet_isFake_missingWJet = selGenWJet_sublead_isFake;
    }

    MEMEvent_singlelepton memEvent_missingWJet(
      eventInfo, isSignal, 
      &memMeasuredBJet_lead, &memMeasuredBJet_sublead,
      memMeasuredWJet_missingWJet, nullptr,
      &memMeasuredLepton,
      genMEt_smeared.px(), genMEt_smeared.py(), metCov);
    addGenMatches_singlelepton(memEvent_missingWJet, genBJetsForMatching_ptrs, genWJetsForMatching_ptrs, genLeptonsForMatching_ptrs, genMEtPx, genMEtPy);

    std::vector<mem::MeasuredParticle> memMeasuredParticles_missingBnWJet;
    memMeasuredParticles_missingBnWJet.push_back(memMeasuredLepton);
    const mem::MeasuredParticle* memMeasuredBJet_missingBnWJet = nullptr;
    bool selGenBJet_isFake_missingBnWJet;
    double u3 = rnd.Uniform();
    assert(u3 >= 0. && u3 <= 1.);
    if ( u3 > 0.50 ) {
      memMeasuredParticles_missingBnWJet.push_back(memMeasuredBJet_lead);
      memMeasuredBJet_missingBnWJet = &memMeasuredBJet_lead;
      selGenBJet_isFake_missingBnWJet = selGenBJet_lead_isFake;
    } else {
      memMeasuredParticles_missingBnWJet.push_back(memMeasuredBJet_sublead);
      memMeasuredBJet_missingBnWJet = &memMeasuredBJet_sublead;
      selGenBJet_isFake_missingBnWJet = selGenBJet_sublead_isFake;
    }
    const mem::MeasuredParticle* memMeasuredWJet_missingBnWJet = nullptr;
    bool selGenWJet_isFake_missingBnWJet;
    double u4 = rnd.Uniform();
    assert(u4 >= 0. && u4 <= 1.);
    if ( u4 > 0.50 ) {
      memMeasuredParticles_missingBnWJet.push_back(memMeasuredWJet_lead);
      memMeasuredWJet_missingBnWJet = &memMeasuredWJet_lead;
      selGenWJet_isFake_missingBnWJet = selGenWJet_lead_isFake;
    } else {
      memMeasuredParticles_missingBnWJet.push_back(memMeasuredWJet_sublead);
      memMeasuredWJet_missingBnWJet = &memMeasuredWJet_sublead;
      selGenWJet_isFake_missingBnWJet = selGenWJet_sublead_isFake;
    }

    MEMEvent_singlelepton memEvent_missingBnWJet(
      eventInfo, isSignal, 
      memMeasuredBJet_missingBnWJet, nullptr,
      memMeasuredWJet_missingBnWJet, nullptr,
      &memMeasuredLepton,
      genMEt_smeared.px(), genMEt_smeared.py(), metCov);
    addGenMatches_singlelepton(memEvent_missingBnWJet, genBJetsForMatching_ptrs, genWJetsForMatching_ptrs, genLeptonsForMatching_ptrs, genMEtPx, genMEtPy);

    const double sqrtS = 13.e+3;
    const std::string pdfName = "MSTW2008lo68cl";
    const std::string madgraphFileName_signal     = "hhAnalysis/bbwwMEM/data/param_hh_SM.dat";
    const std::string madgraphFileName_background = "hhAnalysis/bbwwMEM/data/param_ttbar.dat";
    const bool applyOnshellWmassConstraint_signal = false;
    const int memAlgo_verbosity = 0;
    //const int maxObjFunctionCalls_signal = 2500;
    //const int maxObjFunctionCalls_background = 25000;
    const int maxObjFunctionCalls_signal = 1000;
    const int maxObjFunctionCalls_background = 10000;

    clock.Reset();
    clock.Start("memAlgo");
    MEMbbwwAlgoSingleLepton memAlgo(sqrtS, pdfName, findFile(madgraphFileName_signal), findFile(madgraphFileName_background), memAlgo_verbosity);
    memAlgo.setBJet1TF(&bjetTF);
    memAlgo.setBJet2TF(&bjetTF);
    memAlgo.setHadWJet1TF(&hadWJetTF);
    memAlgo.setHadWJet2TF(&hadWJetTF);
    memAlgo.applyOnshellWmassConstraint_signal(applyOnshellWmassConstraint_signal);
    memAlgo.setIntMode(MEMbbwwAlgoSingleLepton::kVAMP);
    memAlgo.setMaxObjFunctionCalls_signal(maxObjFunctionCalls_signal);
    memAlgo.setMaxObjFunctionCalls_background(maxObjFunctionCalls_background);
    memAlgo.integrate(memMeasuredParticles, genMEt_smeared.px(), genMEt_smeared.py(), metCov);
    MEMbbwwResultSingleLepton memResult = memAlgo.getResult();
    clock.Stop("memAlgo");

    double memCpuTime = clock.GetCpuTime("memAlgo");
    if ( isDEBUG ) {
      std::cout << "MEM:"
	        << " probability for signal hypothesis = " << memResult.getProb_signal() 
                << " +/- " << memResult.getProbErr_signal() << ","
	        << " probability for background hypothesis = " << memResult.getProb_background() 
                << " +/- " << memResult.getProbErr_background() << " " 
	        << "--> likelihood ratio = " << memResult.getLikelihoodRatio() 
                << " +/- " << memResult.getLikelihoodRatioErr() 
	        << " (CPU time = " << memCpuTime << ")" << std::endl;
    }

    (const_cast<MEMEvent_singlelepton*>(&memEvent))->set_memResult(memResult);
    (const_cast<MEMEvent_singlelepton*>(&memEvent))->set_memCpuTime(memCpuTime);

    mem_ntuple->read(memEvent);
    mem_ntuple->fill();

    clock.Reset();
    clock.Start("memAlgo_missingBJet");
    MEMbbwwAlgoSingleLepton memAlgo_missingBJet(sqrtS, pdfName, findFile(madgraphFileName_signal), findFile(madgraphFileName_background), memAlgo_verbosity);
    memAlgo_missingBJet.setBJet1TF(&bjetTF);
    memAlgo_missingBJet.setBJet2TF(&bjetTF);
    memAlgo_missingBJet.setHadWJet1TF(&hadWJetTF);
    memAlgo_missingBJet.setHadWJet2TF(&hadWJetTF);
    memAlgo_missingBJet.applyOnshellWmassConstraint_signal(applyOnshellWmassConstraint_signal);
    memAlgo_missingBJet.setIntMode(MEMbbwwAlgoSingleLepton::kVAMP);
    memAlgo_missingBJet.setMaxObjFunctionCalls_signal(maxObjFunctionCalls_signal);
    memAlgo_missingBJet.setMaxObjFunctionCalls_background(maxObjFunctionCalls_background);
    memAlgo_missingBJet.integrate(memMeasuredParticles_missingBJet, genMEt_smeared.px(), genMEt_smeared.py(), metCov);
    MEMbbwwResultSingleLepton memResult_missingBJet = memAlgo_missingBJet.getResult();
    clock.Stop("memAlgo_missingBJet");
    
    double memCpuTime_missingBJet = clock.GetCpuTime("memAlgo_missingBJet");
    if ( isDEBUG ) {
      std::cout << "MEM (missing b-jet case):" 
	        << " probability for signal hypothesis = " << memResult_missingBJet.getProb_signal() 
                << " +/- " << memResult_missingBJet.getProbErr_signal() << ","
	        << " probability for background hypothesis = " << memResult_missingBJet.getProb_background() 
                << " +/- " << memResult_missingBJet.getProbErr_background() << " " 
	        << "--> likelihood ratio = " << memResult_missingBJet.getLikelihoodRatio() 
                << " +/- " << memResult_missingBJet.getLikelihoodRatioErr() 
	        << " (CPU time = " << memCpuTime_missingBJet << ")" << std::endl;
    }

    (const_cast<MEMEvent_singlelepton*>(&memEvent_missingBJet))->set_memResult(memResult_missingBJet);
    (const_cast<MEMEvent_singlelepton*>(&memEvent_missingBJet))->set_memCpuTime(memCpuTime_missingBJet);

    mem_ntuple_missingBJet->read(memEvent_missingBJet);
    mem_ntuple_missingBJet->fill();

    clock.Reset();
    clock.Start("memAlgo_missingWJet");
    MEMbbwwAlgoSingleLepton memAlgo_missingWJet(sqrtS, pdfName, findFile(madgraphFileName_signal), findFile(madgraphFileName_background), memAlgo_verbosity);
    memAlgo_missingWJet.setBJet1TF(&bjetTF);
    memAlgo_missingWJet.setBJet2TF(&bjetTF);
    memAlgo_missingWJet.setHadWJet1TF(&hadWJetTF);
    memAlgo_missingWJet.setHadWJet2TF(&hadWJetTF);
    memAlgo_missingWJet.applyOnshellWmassConstraint_signal(applyOnshellWmassConstraint_signal);
    memAlgo_missingWJet.setIntMode(MEMbbwwAlgoSingleLepton::kVAMP);
    memAlgo_missingWJet.setMaxObjFunctionCalls_signal(maxObjFunctionCalls_signal);
    memAlgo_missingWJet.setMaxObjFunctionCalls_background(maxObjFunctionCalls_background);
    memAlgo_missingWJet.integrate(memMeasuredParticles_missingWJet, genMEt_smeared.px(), genMEt_smeared.py(), metCov);
    MEMbbwwResultSingleLepton memResult_missingWJet = memAlgo_missingWJet.getResult();
    clock.Stop("memAlgo_missingWJet");
    
    double memCpuTime_missingWJet = clock.GetCpuTime("memAlgo_missingWJet");
    if ( isDEBUG ) {
      std::cout << "MEM (missing jet from W->jj case):" 
	        << " probability for signal hypothesis = " << memResult_missingWJet.getProb_signal() 
                << " +/- " << memResult_missingWJet.getProbErr_signal() << ","
	        << " probability for background hypothesis = " << memResult_missingWJet.getProb_background() 
                << " +/- " << memResult_missingWJet.getProbErr_background() << " " 
	        << "--> likelihood ratio = " << memResult_missingWJet.getLikelihoodRatio() 
                << " +/- " << memResult_missingWJet.getLikelihoodRatioErr() 
	        << " (CPU time = " << memCpuTime_missingWJet << ")" << std::endl;
    }

    (const_cast<MEMEvent_singlelepton*>(&memEvent_missingWJet))->set_memResult(memResult_missingWJet);
    (const_cast<MEMEvent_singlelepton*>(&memEvent_missingWJet))->set_memCpuTime(memCpuTime_missingWJet);

    mem_ntuple_missingWJet->read(memEvent_missingWJet);
    mem_ntuple_missingWJet->fill();

    clock.Reset();
    clock.Start("memAlgo_missingBnWJet");
    MEMbbwwAlgoSingleLepton memAlgo_missingBnWJet(sqrtS, pdfName, findFile(madgraphFileName_signal), findFile(madgraphFileName_background), memAlgo_verbosity);
    memAlgo_missingBnWJet.setBJet1TF(&bjetTF);
    memAlgo_missingBnWJet.setBJet2TF(&bjetTF);
    memAlgo_missingBnWJet.setHadWJet1TF(&hadWJetTF);
    memAlgo_missingBnWJet.setHadWJet2TF(&hadWJetTF);
    memAlgo_missingBnWJet.applyOnshellWmassConstraint_signal(applyOnshellWmassConstraint_signal);
    memAlgo_missingBnWJet.setIntMode(MEMbbwwAlgoSingleLepton::kVAMP);
    memAlgo_missingBnWJet.setMaxObjFunctionCalls_signal(maxObjFunctionCalls_signal);
    memAlgo_missingBnWJet.setMaxObjFunctionCalls_background(maxObjFunctionCalls_background);
    memAlgo_missingBnWJet.integrate(memMeasuredParticles_missingBnWJet, genMEt_smeared.px(), genMEt_smeared.py(), metCov);
    MEMbbwwResultSingleLepton memResult_missingBnWJet = memAlgo_missingBnWJet.getResult();
    clock.Stop("memAlgo_missingBnWJet");
    
    double memCpuTime_missingBnWJet = clock.GetCpuTime("memAlgo_missingBnWJet");
    if ( isDEBUG ) {
      std::cout << "MEM (missing b-jet && jet from W->jj case):" 
	        << " probability for signal hypothesis = " << memResult_missingBnWJet.getProb_signal() 
                << " +/- " << memResult_missingBnWJet.getProbErr_signal() << ","
	        << " probability for background hypothesis = " << memResult_missingBnWJet.getProb_background() 
                << " +/- " << memResult_missingBnWJet.getProbErr_background() << " " 
	        << "--> likelihood ratio = " << memResult_missingBnWJet.getLikelihoodRatio() 
                << " +/- " << memResult_missingBnWJet.getLikelihoodRatioErr() 
	        << " (CPU time = " << memCpuTime_missingBnWJet << ")" << std::endl;
    }

    (const_cast<MEMEvent_singlelepton*>(&memEvent_missingBnWJet))->set_memResult(memResult_missingBnWJet);
    (const_cast<MEMEvent_singlelepton*>(&memEvent_missingBnWJet))->set_memCpuTime(memCpuTime_missingBnWJet);

    mem_ntuple_missingBnWJet->read(memEvent_missingBnWJet);
    mem_ntuple_missingBnWJet->fill();
    //---------------------------------------------------------------------------

    int numGenuineBJets = 0;
    if ( !selGenBJet_lead_isFake    ) ++numGenuineBJets;
    if ( !selGenBJet_sublead_isFake ) ++numGenuineBJets;
    int numGenuineWJets = 0;
    if ( !selGenWJet_lead_isFake    ) ++numGenuineWJets;
    if ( !selGenWJet_sublead_isFake ) ++numGenuineWJets;
    if ( numGenuineBJets == 2 && numGenuineWJets == 2 ) {
      selHistManager->mem_2genuineBJets_2genuineWJets_->fillHistograms(memResult, memCpuTime, evtWeight);
    } else if ( numGenuineBJets == 1 && numGenuineWJets == 2 ) {
      selHistManager->mem_1genuineBJet_2genuineWJets_->fillHistograms(memResult, memCpuTime, evtWeight);
    } else if ( numGenuineBJets == 2 && numGenuineWJets == 1 ) {
      selHistManager->mem_2genuineBJets_1genuineWJet_->fillHistograms(memResult, memCpuTime, evtWeight);
    } else if ( numGenuineBJets == 1 && numGenuineWJets == 1 ) {
      selHistManager->mem_1genuineBJet_1genuineWJet_->fillHistograms(memResult, memCpuTime, evtWeight);
    } 
    int numGenuineBJets_missingBJet = ( !selGenBJet_isFake_missingBJet ) ? 1 : 0;
    if ( numGenuineBJets_missingBJet == 1 && numGenuineWJets == 2 ) {
      selHistManager->mem_missingBJet_genuineBJet_2genuineWJets_->fillHistograms(memResult_missingBJet, memCpuTime_missingBJet, evtWeight);
    } else if ( numGenuineBJets_missingBJet == 0 && numGenuineWJets == 2 ) {
      selHistManager->mem_missingBJet_fakeBJet_2genuineWJets_->fillHistograms(memResult_missingBJet, memCpuTime_missingBJet, evtWeight);
    }
    int numGenuineWJets_missingWJet = ( !selGenWJet_isFake_missingWJet ) ? 1 : 0;
    if ( numGenuineBJets == 2 && numGenuineWJets_missingWJet == 1 ) {
      selHistManager->mem_missingWJet_2genuineBJets_genuineWJet_->fillHistograms(memResult_missingWJet, memCpuTime_missingWJet, evtWeight);
    } else if ( numGenuineBJets == 2 && numGenuineWJets_missingWJet == 0 ) {
      selHistManager->mem_missingWJet_2genuineBJets_fakeWJet_->fillHistograms(memResult_missingWJet, memCpuTime_missingWJet, evtWeight);
    }
    int numGenuineBJets_missingBnWJet = ( !selGenBJet_isFake_missingBnWJet ) ? 1 : 0;
    int numGenuineWJets_missingBnWJet = ( !selGenWJet_isFake_missingBnWJet ) ? 1 : 0;
    if ( numGenuineBJets_missingBnWJet == 1 && numGenuineWJets_missingBnWJet == 1 ) {
      selHistManager->mem_missingBnWJet_genuineBJet_genuineWJet_->fillHistograms(memResult_missingBnWJet, memCpuTime_missingBnWJet, evtWeight);
    } else if ( numGenuineBJets_missingBnWJet == 0 && numGenuineWJets_missingBnWJet == 1 ) {
      selHistManager->mem_missingBnWJet_fakeBJet_genuineWJet_->fillHistograms(memResult_missingBnWJet, memCpuTime_missingBnWJet, evtWeight);
    } else if ( numGenuineBJets_missingBnWJet == 1 && numGenuineWJets_missingBnWJet == 0 ) {
      selHistManager->mem_missingBnWJet_genuineBJet_fakeWJet_->fillHistograms(memResult_missingBnWJet, memCpuTime_missingBnWJet, evtWeight);
    } else if ( numGenuineBJets_missingBnWJet == 0 && numGenuineWJets_missingBnWJet == 0 ) {
      selHistManager->mem_missingBnWJet_fakeBJet_fakeWJet_->fillHistograms(memResult_missingBnWJet, memCpuTime_missingBnWJet, evtWeight);
    }
    selHistManager->genEvtHistManager_afterCuts_->fillHistograms(genElectrons, genMuons, {}, {}, genJets, evtWeight);
    selHistManager->lheInfoHistManager_afterCuts_->fillHistograms(*lheInfoReader, evtWeight);
    selHistManager->weights_->fillHistograms("genWeight", eventInfo.genWeight);
    selHistManager->weights_->fillHistograms("pileupWeight", eventInfo.pileupWeight);

    if ( selEventsFile ) {
      (*selEventsFile) << eventInfo.run << ':' << eventInfo.lumi << ':' << eventInfo.event << '\n';
    }

    ++selectedEntries;
    selectedEntries_weighted += evtWeight;
    histogram_selectedEntries->Fill(0.);
  }

  std::cout << "max num. Entries = " << inputTree -> getCumulativeMaxEventCount()
            << " (limited by " << maxEvents << ") processed in "
            << inputTree -> getProcessedFileCount() << " file(s) (out of "
            << inputTree -> getFileCount() << ")\n"
            << " analyzed = " << analyzedEntries << '\n'
            << " selected = " << selectedEntries << " (weighted = " << selectedEntries_weighted << ")\n\n"
            << "cut-flow table" << std::endl;
  cutFlowTable.print(std::cout);
  std::cout << std::endl;

  delete run_lumi_eventSelector;

  delete selEventsFile;

  delete genLeptonReader;
  delete genNeutrinoReader;
  delete genJetReader;
  delete lheInfoReader;

  delete genEvtHistManager_beforeCuts;

  delete mem_ntuple;
  delete mem_ntuple_missingBJet;
  delete mem_ntuple_missingWJet;
  delete mem_ntuple_missingBnWJet;

  delete inputTree;

  clock.Show("analyze_hh_bbwwMEM_singlelepton");

  return EXIT_SUCCESS;
}
