#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/MEMbbwwNtupleManager_dilepton.h"

#include "tthAnalysis/HiggsToTauTau/interface/TypeTraits.h"           // Traits<T>::TYPE_NAME
#include "tthAnalysis/HiggsToTauTau/interface/mvaInputVariables.h"    // comp_MT_met

#include <DataFormats/Math/interface/LorentzVector.h> // math::PtEtaPhiMLorentzVector

#include <algorithm> // std::sort
#include <math.h>    // sqrt, atan2

MEMbbwwNtupleManager_dilepton::MEMbbwwNtupleManager_dilepton(const std::string & outputTreeName)
  : MEMbbwwNtupleManager(outputTreeName)
  , lepton1_("lepton1")
  , lepton2_("lepton2")
  , gen_lepton1_("gen_lepton1")
  , gen_lepton2_("gen_lepton2")
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

  met_.initializeBranches(tree_);
  gen_met_.initializeBranches(tree_);

  tree_->Branch("ptww",   &ptww_,   Form("ptww/%s",   Traits<Float_t>::TYPE_NAME));
  tree_->Branch("mww",    &mww_,    Form("mww/%s",    Traits<Float_t>::TYPE_NAME));
  tree_->Branch("ptll",   &ptll_,   Form("ptll/%s",   Traits<Float_t>::TYPE_NAME));
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
  mll_          = 0.;
  ptmiss_       = 0.;

}

void 
MEMbbwwNtupleManager_dilepton::read(const std::vector<mem::MeasuredParticle> & measuredParticles,
                                    double measuredMEtPx, double measuredMEtPy, const TMatrixD & measuredMEtCov)
{
  std::vector<const mem::MeasuredParticle *> measuredBJets;
  std::vector<const mem::MeasuredParticle *> measuredLeptons;
  for ( std::vector<mem::MeasuredParticle>::const_iterator measuredParticle = measuredParticles.begin();
        measuredParticle != measuredParticles.end(); ++measuredParticle ) {
    if      ( measuredParticle->type() == mem::MeasuredParticle::kElectron ) measuredLeptons.push_back(&(*measuredParticle));
    else if ( measuredParticle->type() == mem::MeasuredParticle::kMuon     ) measuredLeptons.push_back(&(*measuredParticle));
    else if ( measuredParticle->type() == mem::MeasuredParticle::kBJet     ) measuredBJets.push_back(&(*measuredParticle));
  }

  std::sort(measuredBJets.begin(), measuredBJets.end(), mem::isHigherPt);
  if ( measuredBJets.size() >= 1 ) bjet1_.read(measuredBJets[0]);
  if ( measuredBJets.size() >= 2 ) bjet2_.read(measuredBJets[1]);
  nbjets_ = countMeasuredJets(&bjet1_, &bjet2_);
  //std::cout << "nbjets = " << nbjets_ << std::endl;

  std::sort(measuredLeptons.begin(), measuredLeptons.end(), mem::isHigherPt);
  if ( measuredLeptons.size() >= 1 ) lepton1_.read(measuredLeptons[0]);
  if ( measuredLeptons.size() >= 2 ) lepton2_.read(measuredLeptons[1]);
  nleptons_ = countMeasuredLeptons(&lepton1_, &lepton2_);
  //std::cout << "nleptons = " << nleptons_ << std::endl;

  met_.read(measuredMEtPx, measuredMEtPy, &measuredMEtCov);

  if ( bjet1_.measuredJet_ && bjet2_.measuredJet_ )
  {
    const mem::MeasuredParticle * bjet1 = bjet1_.measuredJet_;
    math::PtEtaPhiMLorentzVector bjet1P4(bjet1->pt(), bjet1->eta(), bjet1->phi(), bjet1->mass());
    const mem::MeasuredParticle * bjet2 = bjet2_.measuredJet_;
    math::PtEtaPhiMLorentzVector bjet2P4(bjet2->pt(), bjet2->eta(), bjet2->phi(), bjet2->mass());
    math::PtEtaPhiMLorentzVector hbbP4 = bjet1P4 + bjet2P4;
    ptbb_ = hbbP4.pt();
    drbb_ = deltaR(bjet1P4, bjet2P4); 
    mbb_ = hbbP4.mass();
  }

  if ( lepton1_.measuredLepton_ && lepton2_.measuredLepton_ )
  {
    const mem::MeasuredParticle * lepton1 = lepton1_.measuredLepton_;
    math::PtEtaPhiMLorentzVector lepton1P4(lepton1->cone_pt(), lepton1->eta(), lepton1->phi(), lepton1->mass());
    const mem::MeasuredParticle * lepton2 = lepton2_.measuredLepton_;
    math::PtEtaPhiMLorentzVector lepton2P4(lepton2->cone_pt(), lepton2->eta(), lepton2->phi(), lepton2->mass());
    double metPt = sqrt(met_.px_*met_.px_ + met_.py_*met_.py_);
    double metPhi = atan2(met_.py_, met_.px_);
    math::PtEtaPhiMLorentzVector metP4(metPt, 0., metPhi, 0.);
    math::PtEtaPhiMLorentzVector hwwP4 = lepton1P4 + lepton2P4 + metP4;
    ptww_ = hwwP4.pt();
    mww_ = hwwP4.mass();
    math::PtEtaPhiMLorentzVector dileptonP4 = lepton1P4 + lepton2P4;
    ptll_ = dileptonP4.pt();
    mll_ = dileptonP4.mass();
    ptmiss_ = metPt;
  }
}

void 
MEMbbwwNtupleManager_dilepton::read(const std::vector<const GenJet *> & genBJets,
                                    const std::vector<const GenLepton *> & genLeptons,
                                    double genMEtPx, double genMEtPy)
{
  if ( bjet1_.measuredJet_ ) gen_bjet1_.read(mem::findGenMatch(bjet1_.measuredJet_, genBJets));
  if ( bjet2_.measuredJet_ ) gen_bjet2_.read(mem::findGenMatch(bjet2_.measuredJet_, genBJets));
  gen_nbjets_ = countGenJets(&gen_bjet1_, &gen_bjet2_);

  if ( lepton1_.measuredLepton_ ) gen_lepton1_.read(mem::findGenMatch(lepton1_.measuredLepton_, genLeptons));
  if ( lepton2_.measuredLepton_ ) gen_lepton2_.read(mem::findGenMatch(lepton2_.measuredLepton_, genLeptons));
  gen_nleptons_ = countGenLeptons(&gen_lepton1_, &gen_lepton2_);

  gen_met_.read(genMEtPx, genMEtPy);
}
