#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/MEMbbwwNtupleManager.h"

#include "tthAnalysis/HiggsToTauTau/interface/histogramAuxFunctions.h" // createSubdirectory_recursively

MEMbbwwNtupleManager::MEMbbwwNtupleManager(const std::string & outputTreeName)
  : outputTreeName_(outputTreeName)
  , tree_(nullptr)
  , run_(0)
  , ls_(0)
  , event_(0)
  , genWeight_(0.)
  , isSignal_(false)
  , memProbS_(0.)
  , memProbSerr_(0.)
  , memProbB_(0.)
  , memProbBerr_(0.)
  , memLR_(0.)
  , memLRerr_(0.)
  , memCpuTime_(0.)
  , bjet1_("bjet1")
  , bjet2_("bjet2")
  , nbjets_(0)
  , gen_bjet1_("gen_bjet1")
  , gen_bjet2_("gen_bjet2")
  , gen_nbjets_(0)
  , met_("met")
  , gen_met_("gen_met")
  , barcode_(-1)
  , ptbb_(0.)
  , drbb_(0.)
  , mbb_(0.)
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
  tree_->Branch("ls",          &ls_,          Form("ls/%s",          Traits<UInt_t>::TYPE_NAME));
  tree_->Branch("event",       &event_,       Form("event/%s",       Traits<ULong64_t>::TYPE_NAME));
  
  tree_->Branch("memProbS",    &memProbS_,    Form("memProbS/%s",    Traits<Double_t>::TYPE_NAME));
  tree_->Branch("memProbSerr", &memProbSerr_, Form("memProbSerr/%s", Traits<Double_t>::TYPE_NAME));
  tree_->Branch("memProbB",    &memProbB_,    Form("memProbB/%s",    Traits<Double_t>::TYPE_NAME));
  tree_->Branch("memProbBerr", &memProbBerr_, Form("memProbBerr/%s", Traits<Double_t>::TYPE_NAME));
  tree_->Branch("memLR",       &memLR_,       Form("memLR/%s",       Traits<Double_t>::TYPE_NAME));
  tree_->Branch("memLRerr",    &memLRerr_,    Form("memLRerr/%s",    Traits<Double_t>::TYPE_NAME));
  tree_->Branch("memCpuTime",  &memCpuTime_,  Form("memCpuTime/%s",  Traits<Float_t>::TYPE_NAME));

  tree_->Branch("genWeight",   &genWeight_,   Form("genWeight/%s",   Traits<Float_t>::TYPE_NAME));

  tree_->Branch("isSignal",    &isSignal_,    Form("isSignal/%s",    Traits<Bool_t>::TYPE_NAME));

  bjet1_.initializeBranches(tree_);
  bjet2_.initializeBranches(tree_);
  tree_->Branch("nbjets",      &nbjets_,      Form("nbjets/%s",      Traits<Int_t>::TYPE_NAME));
  gen_bjet1_.initializeBranches(tree_);
  gen_bjet2_.initializeBranches(tree_);
  tree_->Branch("gen_nbjets",  &gen_nbjets_,  Form("gen_nbjets/%s",  Traits<Int_t>::TYPE_NAME));
  
  met_.initializeBranches(tree_);
  gen_met_.initializeBranches(tree_);

  tree_->Branch("barcode",     &barcode_,     Form("barcode/%s",     Traits<Int_t>::TYPE_NAME));

  tree_->Branch("ptbb",        &ptbb_,        Form("ptbb/%s",        Traits<Float_t>::TYPE_NAME));
  tree_->Branch("drbb",        &drbb_,        Form("drbb/%s",        Traits<Float_t>::TYPE_NAME));
  tree_->Branch("mbb",         &mbb_,         Form("mbb/%s",         Traits<Float_t>::TYPE_NAME));
}

void 
MEMbbwwNtupleManager::read(const MEMEvent & memEvent)
{
  run_         = memEvent.eventInfo().run;
  ls_          = memEvent.eventInfo().lumi;
  event_       = memEvent.eventInfo().event;

  genWeight_   = memEvent.eventInfo().genWeight;

  isSignal_    = memEvent.isSignal();

  memProbS_    = memEvent.memResult().getProb_signal();
  memProbSerr_ = memEvent.memResult().getProbErr_signal();
  memProbB_    = memEvent.memResult().getProb_background();
  memProbBerr_ = memEvent.memResult().getProbErr_background();
  memLR_       = memEvent.memResult().getLikelihoodRatio();
  memLRerr_    = memEvent.memResult().getLikelihoodRatioErr();
  memCpuTime_  = memEvent.memCpuTime();

  bjet1_.read(memEvent.measuredBJet1(), memEvent.genBJet1() != nullptr);
  bjet2_.read(memEvent.measuredBJet2(), memEvent.genBJet2() != nullptr);
  nbjets_      = memEvent.numMeasuredBJets();
  gen_bjet1_.read(memEvent.genBJet1());
  gen_bjet2_.read(memEvent.genBJet2());
  gen_nbjets_  = memEvent.numGenBJets();

  met_.read(memEvent.measuredMEtPx(), memEvent.measuredMEtPy(), &memEvent.measuredMEtCov());
  gen_met_.read(memEvent.genMEtPx(), memEvent.genMEtPy());

  barcode_     = memEvent.barcode();

  if ( memEvent.measuredBJet1() && memEvent.measuredBJet2() )
  {
    const mem::MeasuredParticle * bjet1 = memEvent.measuredBJet1();
    math::PtEtaPhiMLorentzVector bjet1P4(bjet1->pt(), bjet1->eta(), bjet1->phi(), bjet1->mass());
    const mem::MeasuredParticle * bjet2 = memEvent.measuredBJet2();
    math::PtEtaPhiMLorentzVector bjet2P4(bjet2->pt(), bjet2->eta(), bjet2->phi(), bjet2->mass());
    math::PtEtaPhiMLorentzVector hbbP4 = bjet1P4 + bjet2P4;
    ptbb_      = hbbP4.pt();
    drbb_      = deltaR(bjet1P4, bjet2P4); 
    mbb_       = hbbP4.mass();
  }
}

void MEMbbwwNtupleManager::fill()
{
  tree_->Fill();
  resetBranches();
}

void 
MEMbbwwNtupleManager::resetBranches()
{
  run_         = 0;
  ls_          = 0;
  event_       = 0;

  genWeight_   = 0.;

  isSignal_    = false;

  memProbS_    = 0.;
  memProbSerr_ = 0.;
  memProbB_    = 0.;
  memProbBerr_ = 0.;
  memLR_       = 0.;
  memLRerr_    = 0.;
  memCpuTime_  = -1.;

  bjet1_.resetBranches();
  bjet2_.resetBranches();
  nbjets_      = 0;
  gen_bjet1_.resetBranches();
  gen_bjet2_.resetBranches();
  gen_nbjets_  = 0;

  met_.resetBranches();
  gen_met_.resetBranches();

  barcode_     = -1;

  ptbb_        = 0.;
  drbb_        = 0.;
  mbb_         = 0.;
}
