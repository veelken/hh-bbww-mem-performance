#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/MEMbbwwNtupleManager_dilepton.h"

#include "tthAnalysis/HiggsToTauTau/interface/TypeTraits.h"           // Traits<T>::TYPE_NAME
#include "tthAnalysis/HiggsToTauTau/interface/mvaInputVariables.h"    // comp_MT_met

#include "DataFormats/Math/interface/LorentzVector.h" // math::PtEtaPhiMLorentzVector
#include "DataFormats/Math/interface/deltaR.h"        // deltaR
#include "DataFormats/Math/interface/deltaPhi.h"      // deltaPhi

#include <algorithm> // std::sort
#include <math.h>    // sqrt, atan2

MEMbbwwNtupleManager_dilepton::MEMbbwwNtupleManager_dilepton(const std::string & outputDirectoryName, const std::string & outputTreeName)
  : MEMbbwwNtupleManager(outputDirectoryName, outputTreeName)
  , lepton1_("lepton1")
  , lepton2_("lepton2")
  , nleptons_(0)
  , gen_lepton1_("gen_lepton1")
  , gen_lepton2_("gen_lepton2")
  , gen_nleptons_(0)
  , ptww_(0.)
  , mww_(0.)
  , ptll_(0.)
  , drll_(0.)
  , dphill_(0.)
  , mll_(0.)
  , ptmiss_(0.)
{}

MEMbbwwNtupleManager_dilepton::~MEMbbwwNtupleManager_dilepton()
{}

void 
MEMbbwwNtupleManager_dilepton::initializeBranches()
{
  MEMbbwwNtupleManager::initializeBranches();

  assert(tree_);

  lepton1_.initializeBranches(tree_);
  lepton2_.initializeBranches(tree_);
  tree_->Branch("nleptons",     &nleptons_,     Form("nleptons/%s",     Traits<Int_t>::TYPE_NAME));
  gen_lepton1_.initializeBranches(tree_);
  gen_lepton2_.initializeBranches(tree_);
  tree_->Branch("gen_nleptons", &gen_nleptons_, Form("gen_nleptons/%s", Traits<Int_t>::TYPE_NAME));
  
  tree_->Branch("ptww",   &ptww_,   Form("ptww/%s",   Traits<Float_t>::TYPE_NAME));
  tree_->Branch("mww",    &mww_,    Form("mww/%s",    Traits<Float_t>::TYPE_NAME));
  tree_->Branch("ptll",   &ptll_,   Form("ptll/%s",   Traits<Float_t>::TYPE_NAME));
  tree_->Branch("drll",   &drll_,   Form("drll/%s",   Traits<Float_t>::TYPE_NAME));
  tree_->Branch("dphill", &dphill_, Form("dphill/%s", Traits<Float_t>::TYPE_NAME));  
  tree_->Branch("mll",    &mll_,    Form("mll/%s",    Traits<Float_t>::TYPE_NAME));
  tree_->Branch("ptmiss", &ptmiss_, Form("ptmiss/%s", Traits<Float_t>::TYPE_NAME));
}

void 
MEMbbwwNtupleManager_dilepton::resetBranches()
{
  MEMbbwwNtupleManager::resetBranches();

  lepton1_.resetBranches();
  lepton2_.resetBranches();
  nleptons_     = 0;
  gen_lepton1_.resetBranches();
  gen_lepton2_.resetBranches();
  gen_nleptons_ = 0;

  ptww_         = 0.;
  mww_          = 0.;
  ptll_         = 0.;
  drll_         = 0.;
  dphill_       = 0.;
  mll_          = 0.;
  ptmiss_       = 0.;
}

void 
MEMbbwwNtupleManager_dilepton::read(const MEMEvent_dilepton & memEvent)
{
  MEMbbwwNtupleManager::read(memEvent);

  lepton1_.read(memEvent.measuredLepton1(), memEvent.genLepton1() != nullptr);
  lepton2_.read(memEvent.measuredLepton2(), memEvent.genLepton2() != nullptr);
  nleptons_     = memEvent.numMeasuredLeptons();
  gen_lepton1_.read(memEvent.genLepton1());
  gen_lepton2_.read(memEvent.genLepton2());
  gen_nleptons_ = memEvent.numGenLeptons();

  if ( memEvent.measuredLepton1() && memEvent.measuredLepton2() )
  {
    const mem::MeasuredParticle * lepton1 = memEvent.measuredLepton1();
    math::PtEtaPhiMLorentzVector lepton1P4(lepton1->pt(), lepton1->eta(), lepton1->phi(), lepton1->mass());
    const mem::MeasuredParticle * lepton2 = memEvent.measuredLepton2();
    math::PtEtaPhiMLorentzVector lepton2P4(lepton2->pt(), lepton2->eta(), lepton2->phi(), lepton2->mass());
    double metPx = memEvent.measuredMEtPx();
    double metPy = memEvent.measuredMEtPy();
    double metPt = sqrt(metPx*metPx + metPy*metPy);
    double metPhi = atan2(metPy, metPx);
    math::PtEtaPhiMLorentzVector metP4(metPt, 0., metPhi, 0.);
    math::PtEtaPhiMLorentzVector hwwP4 = lepton1P4 + lepton2P4 + metP4;
    ptww_       = hwwP4.pt();
    mww_        = hwwP4.mass();
    math::PtEtaPhiMLorentzVector dileptonP4 = lepton1P4 + lepton2P4;
    ptll_       = dileptonP4.pt();
    drll_       = deltaR(lepton1P4, lepton2P4);
    dphill_     = deltaPhi(lepton1P4.phi(), lepton2P4.phi());
    mll_        = dileptonP4.mass();
    ptmiss_     = metPt;
  }
}
