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
  , gen_wjet1_("gen_wjet1")
  , gen_wjet2_("gen_wjet2")
  , lepton_("lepton")
  , gen_lepton_("gen_lepton")
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

  met_.initializeBranches(tree_);
  gen_met_.initializeBranches(tree_);

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
  wjet1_.resetBranches();
  wjet2_.resetBranches();
  nwjets_     = 0;
  gen_wjet1_.resetBranches();
  gen_wjet2_.resetBranches();
  gen_nwjets_ = 0;

  lepton_.resetBranches();
  gen_lepton_.resetBranches();

  met_.resetBranches();
  gen_met_.resetBranches();

  ptjj_       = 0.;
  drjj_       = 0.;
  mjj_        = 0.;
  ptww_       = 0.;
  mww_        = 0.;
  mt_         = 0.;
  ptmiss_     = 0.;

}

void 
MEMbbwwNtupleManager_singlelepton::read(const std::vector<mem::MeasuredParticle> & measuredParticles,
                                        double measuredMEtPx, double measuredMEtPy, const TMatrixD & measuredMEtCov)
{
  std::vector<const mem::MeasuredParticle *> measuredBJets;
  std::vector<const mem::MeasuredParticle *> measuredWJets;
  std::vector<const mem::MeasuredParticle *> measuredLeptons;
  for ( std::vector<mem::MeasuredParticle>::const_iterator measuredParticle = measuredParticles.begin();
        measuredParticle != measuredParticles.end(); ++measuredParticle ) {
    if      ( measuredParticle->type() == mem::MeasuredParticle::kElectron ) measuredLeptons.push_back(&(*measuredParticle));
    else if ( measuredParticle->type() == mem::MeasuredParticle::kMuon     ) measuredLeptons.push_back(&(*measuredParticle));
    else if ( measuredParticle->type() == mem::MeasuredParticle::kBJet     ) measuredBJets.push_back(&(*measuredParticle));
    else if ( measuredParticle->type() == mem::MeasuredParticle::kHadWJet  ) measuredWJets.push_back(&(*measuredParticle));
  }

  std::sort(measuredBJets.begin(), measuredBJets.end(), mem::isHigherPt);
  if ( measuredBJets.size() >= 1 ) bjet1_.read(measuredBJets[0]);
  if ( measuredBJets.size() >= 2 ) bjet2_.read(measuredBJets[1]);
  nbjets_ = countMeasuredJets(&bjet1_, &bjet2_);
  
  std::sort(measuredWJets.begin(), measuredWJets.end(), mem::isHigherPt);
  if ( measuredWJets.size() >= 1 ) wjet1_.read(measuredWJets[0]);
  if ( measuredWJets.size() >= 2 ) wjet2_.read(measuredWJets[1]);
  nwjets_ = countMeasuredJets(&wjet1_, &wjet2_);

  std::sort(measuredLeptons.begin(), measuredLeptons.end(), mem::isHigherPt);
  if ( measuredLeptons.size() >= 1 ) lepton_.read(measuredLeptons[0]);
  nleptons_ = countMeasuredLeptons(&lepton_);

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

  if ( wjet1_.measuredJet_ && wjet2_.measuredJet_ )
  {
    const mem::MeasuredParticle * wjet1 = wjet1_.measuredJet_;
    math::PtEtaPhiMLorentzVector wjet1P4(wjet1->pt(), wjet1->eta(), wjet1->phi(), wjet1->mass());
    const mem::MeasuredParticle * wjet2 = wjet2_.measuredJet_;
    math::PtEtaPhiMLorentzVector wjet2P4(wjet2->pt(), wjet2->eta(), wjet2->phi(), wjet2->mass());
    math::PtEtaPhiMLorentzVector whadP4 = wjet1P4 + wjet2P4;
    ptjj_ = whadP4.pt();
    drjj_ = deltaR(wjet1P4, wjet2P4); 
    mjj_ = whadP4.mass();
    if ( lepton_.measuredLepton_ )
    {
      const mem::MeasuredParticle * lepton = lepton_.measuredLepton_;
      math::PtEtaPhiMLorentzVector leptonP4(lepton->pt(), lepton->eta(), lepton->phi(), lepton->mass());
      double metPt = sqrt(met_.px_*met_.px_ + met_.py_*met_.py_);
      double metPhi = atan2(met_.py_, met_.px_);
      math::PtEtaPhiMLorentzVector metP4(metPt, 0., metPhi, 0.);
      math::PtEtaPhiMLorentzVector hwwP4 = whadP4 + leptonP4 + metP4;
      ptww_ = hwwP4.pt();
      mww_ = hwwP4.mass();
      mt_ = comp_MT_met(leptonP4, metPt, metPhi);
      ptmiss_ = metPt;
    }
  }
}

void 
MEMbbwwNtupleManager_singlelepton::read(const std::vector<GenJet *> & genBJets,
                                        const std::vector<GenJet *> & genWJets,
                                        const GenLepton * genLepton,
                                        double genMEtPx, double genMEtPy)
{
  if ( bjet1_.measuredJet_ ) gen_bjet1_.read(findGenMatch(bjet1_.measuredJet_, genBJets));
  if ( bjet2_.measuredJet_ ) gen_bjet2_.read(findGenMatch(bjet2_.measuredJet_, genBJets));
  gen_nbjets_ = countGenJets(&gen_bjet1_, &gen_bjet2_);

  if ( wjet1_.measuredJet_ ) gen_wjet1_.read(findGenMatch(wjet1_.measuredJet_, genWJets));
  if ( wjet2_.measuredJet_ ) gen_wjet2_.read(findGenMatch(wjet2_.measuredJet_, genWJets));
  gen_nwjets_ = countGenJets(&gen_wjet1_, &gen_wjet2_);

  if ( lepton_.measuredLepton_ )
  {
    std::vector<GenLepton *> genLeptons = { const_cast<GenLepton *>(genLepton) };
    gen_lepton_.read(findGenMatch(lepton_.measuredLepton_, genLeptons));
  }
  gen_nleptons_ = countGenLeptons(&gen_lepton_);

  gen_met_.read(genMEtPx, genMEtPy);
}
