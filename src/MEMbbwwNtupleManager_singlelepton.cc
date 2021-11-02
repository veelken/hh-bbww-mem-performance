#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/MEMbbwwNtupleManager_singlelepton.h"

#include "tthAnalysis/HiggsToTauTau/interface/TypeTraits.h"           // Traits<T>::TYPE_NAME
#include "tthAnalysis/HiggsToTauTau/interface/mvaInputVariables.h"    // comp_MT_met

#include <DataFormats/Math/interface/LorentzVector.h> // math::PtEtaPhiMLorentzVector

#include <algorithm> // std::sort
#include <math.h>    // sqrt, atan2

MEMbbwwNtupleManager_singlelepton::MEMbbwwNtupleManager_singlelepton(const std::string & outputTreeName)
  : MEMbbwwNtupleManager(outputTreeName)
  , wjet1_("wjet1")
  , wjet2_("wjet2")
  , nwjets_(0)
  , gen_wjet1_("gen_wjet1")
  , gen_wjet2_("gen_wjet2")
  , gen_nwjets_(0)
  , lepton_("lepton")
  , nleptons_(0)
  , gen_lepton_("gen_lepton")
  , gen_nleptons_(0)
  , ptjj_(0.)
  , drjj_(0.)
  , mjj_(0.)
  , ptww_(0.)
  , mww_(0.)
  , mt_(0.)
  , ptmiss_(0.)
{}

MEMbbwwNtupleManager_singlelepton::~MEMbbwwNtupleManager_singlelepton()
{}

void 
MEMbbwwNtupleManager_singlelepton::initializeBranches()
{
  MEMbbwwNtupleManager::initializeBranches();

  assert(tree_);

  wjet1_.initializeBranches(tree_);
  wjet2_.initializeBranches(tree_);
  tree_->Branch("nwjets",       &nwjets_,       Form("nwjets/%s",       Traits<Int_t>::TYPE_NAME));
  gen_wjet1_.initializeBranches(tree_);
  gen_wjet2_.initializeBranches(tree_);
  tree_->Branch("gen_nwjets",   &gen_nwjets_,   Form("gen_nwjets/%s",   Traits<Int_t>::TYPE_NAME));
  
  lepton_.initializeBranches(tree_);
  tree_->Branch("nleptons",     &nleptons_,     Form("nleptons/%s",     Traits<Int_t>::TYPE_NAME));
  gen_lepton_.initializeBranches(tree_);
  tree_->Branch("gen_nleptons", &gen_nleptons_, Form("gen_nleptons/%s", Traits<Int_t>::TYPE_NAME));

  tree_->Branch("ptjj",       &ptjj_,       Form("ptjj/%s",       Traits<Float_t>::TYPE_NAME));
  tree_->Branch("drjj",       &drjj_,       Form("drjj/%s",       Traits<Float_t>::TYPE_NAME));
  tree_->Branch("mjj",        &mjj_,        Form("mjj/%s",        Traits<Float_t>::TYPE_NAME));
  tree_->Branch("ptww",       &ptww_,       Form("ptww/%s",       Traits<Float_t>::TYPE_NAME));
  tree_->Branch("mww",        &mww_,        Form("mww/%s",        Traits<Float_t>::TYPE_NAME));
  tree_->Branch("mt",         &mt_,         Form("mt/%s",         Traits<Float_t>::TYPE_NAME));
  tree_->Branch("ptmiss",     &ptmiss_,     Form("ptmiss/%s",     Traits<Float_t>::TYPE_NAME));
}

void 
MEMbbwwNtupleManager_singlelepton::resetBranches()
{
  MEMbbwwNtupleManager::resetBranches();

  wjet1_.resetBranches();
  wjet2_.resetBranches();
  nwjets_     = 0;
  gen_wjet1_.resetBranches();
  gen_wjet2_.resetBranches();
  gen_nwjets_ = 0;

  lepton_.resetBranches();
  gen_lepton_.resetBranches();

  ptjj_       = 0.;
  drjj_       = 0.;
  mjj_        = 0.;
  ptww_       = 0.;
  mww_        = 0.;
  mt_         = 0.;
  ptmiss_     = 0.;
}

void 
MEMbbwwNtupleManager_singlelepton::read(const MEMEvent_singlelepton & memEvent)
{
  wjet1_.read(memEvent.measuredWJet1(), memEvent.genWJet1() != nullptr);
  wjet2_.read(memEvent.measuredWJet2(), memEvent.genWJet2() != nullptr);
  nwjets_     = memEvent.numMeasuredWJets();
  gen_wjet1_.read(memEvent.genWJet1());
  gen_wjet2_.read(memEvent.genWJet2());
  gen_nwjets_ = memEvent.numGenWJets();

  if ( memEvent.measuredWJet1() && memEvent.measuredWJet2() )
  {
    const mem::MeasuredParticle * wjet1 = memEvent.measuredWJet1();
    math::PtEtaPhiMLorentzVector wjet1P4(wjet1->pt(), wjet1->eta(), wjet1->phi(), wjet1->mass());
    const mem::MeasuredParticle * wjet2 = memEvent.measuredWJet2();
    math::PtEtaPhiMLorentzVector wjet2P4(wjet2->pt(), wjet2->eta(), wjet2->phi(), wjet2->mass());
    math::PtEtaPhiMLorentzVector whadP4 = wjet1P4 + wjet2P4;
    ptjj_     = whadP4.pt();
    drjj_     = deltaR(wjet1P4, wjet2P4); 
    mjj_      = whadP4.mass();
    if ( memEvent.measuredLepton() )
    {
      const mem::MeasuredParticle * lepton = memEvent.measuredLepton();
      math::PtEtaPhiMLorentzVector leptonP4(lepton->pt(), lepton->eta(), lepton->phi(), lepton->mass());
      double metPx = memEvent.measuredMEtPx();
      double metPy = memEvent.measuredMEtPy();
      double metPt = sqrt(metPx*metPx + metPy*metPy);
      double metPhi = atan2(metPy, metPx);
      math::PtEtaPhiMLorentzVector metP4(metPt, 0., metPhi, 0.);
      math::PtEtaPhiMLorentzVector hwwP4 = whadP4 + leptonP4 + metP4;
      ptww_   = hwwP4.pt();
      mww_    = hwwP4.mass();
      mt_     = comp_MT_met(leptonP4, metPt, metPhi);
      ptmiss_ = metPt;
    }
  }
}
