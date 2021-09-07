#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/MEMbbwwNtupleManager.h"

#include "tthAnalysis/HiggsToTauTau/interface/histogramAuxFunctions.h" // createSubdirectory_recursively

MEMbbwwNtupleManager::MEMbbwwNtupleManager(const std::string & outputTreeName)
  : outputTreeName_(outputTreeName)
  , tree_(nullptr)
  , bjet1_("bjet1")
  , bjet2_("bjet2")
  , gen_bjet1_("gen_bjet1")
  , gen_bjet2_("gen_bjet2")
  , met_("met")
  , gen_met_("gen_met")
{}

MEMbbwwNtupleManager::~MEMbbwwNtupleManager()
{}

void 
MEMbbwwNtupleManager::makeTree(TFileDirectory & dir)
{
  TDirectory * subDir = createSubdirectory_recursively(dir, "ntuples");
  subDir->cd();
  tree_ = new TTree(outputTreeName_.c_str(), outputTreeName_.c_str());
  dir.cd();
}

void 
MEMbbwwNtupleManager::initializeBranches()
{
  assert(tree_);

  tree_->Branch("run",         &run_,         Form("run/%s",         Traits<UInt_t>::TYPE_NAME));
  tree_->Branch("ls",          &ls_,          Form("run/%s",         Traits<UInt_t>::TYPE_NAME));
  tree_->Branch("event",       &event_,       Form("run/%s",         Traits<ULong64_t>::TYPE_NAME));
  
  tree_->Branch("memProbS",    &memProbS_,    Form("memProbS/%s",    Traits<Double_t>::TYPE_NAME));
  tree_->Branch("memProbSerr", &memProbSerr_, Form("memProbSerr/%s", Traits<Double_t>::TYPE_NAME));
  tree_->Branch("memProbB",    &memProbB_,    Form("memProbB/%s",    Traits<Double_t>::TYPE_NAME));
  tree_->Branch("memProbBerr", &memProbBerr_, Form("memProbBerr/%s", Traits<Double_t>::TYPE_NAME));
  tree_->Branch("memLR",       &memLR_,       Form("memLR/%s",       Traits<Double_t>::TYPE_NAME));
  tree_->Branch("memLRerr",    &memLRerr_,    Form("memLRerr/%s",    Traits<Double_t>::TYPE_NAME));

  bjet1_.initializeBranches(tree_);
  bjet2_.initializeBranches(tree_);
  tree_->Branch("nbjets",      &nbjets_,      Form("nbjets/%s",      Traits<Int_t>::TYPE_NAME));
  gen_bjet1_.initializeBranches(tree_);
  gen_bjet2_.initializeBranches(tree_);
  tree_->Branch("gen_nbjets",  &gen_nbjets_,  Form("gen_nbjets/%s",  Traits<Int_t>::TYPE_NAME));
  
  tree_->Branch("ptbb",        &ptbb_,        Form("ptbb/%s",        Traits<Float_t>::TYPE_NAME));
  tree_->Branch("drbb",        &drbb_,        Form("drbb/%s",        Traits<Float_t>::TYPE_NAME));
  tree_->Branch("mbb",         &mbb_,         Form("mbb/%s",         Traits<Float_t>::TYPE_NAME));
}

void 
MEMbbwwNtupleManager::read(const EventInfo & eventInfo)
{
  run_         = eventInfo.run;
  ls_          = eventInfo.lumi;
  event_       = eventInfo.event;
}

void 
MEMbbwwNtupleManager::read(const MEMResultBase & memResult)
{
  memProbS_    = memResult.getProb_signal();
  memProbSerr_ = memResult.getProbErr_signal();
  memProbB_    = memResult.getProb_background();
  memProbBerr_ = memResult.getProbErr_background();
  memLR_       = memResult.getLikelihoodRatio();
  memLRerr_    = memResult.getLikelihoodRatioErr();
}

void 
MEMbbwwNtupleManager::resetBranches()
{
  run_         = 0;
  ls_          = 0;
  event_       = 0;

  memProbS_    = 0.;
  memProbSerr_ = 0.;
  memProbB_    = 0.;
  memProbBerr_ = 0.;
  memLR_       = 0.;
  memLRerr_    = 0.;

  bjet1_.resetBranches();
  bjet2_.resetBranches();
  nbjets_      = 0;
  gen_bjet1_.resetBranches();
  gen_bjet2_.resetBranches();
  gen_nbjets_  = 0;

  ptbb_        = 0.;
  drbb_        = 0.;
  mbb_         = 0.;
}
