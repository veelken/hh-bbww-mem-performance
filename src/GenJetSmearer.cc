#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/GenJetSmearer.h"

#include <TMath.h> // TMath::Sqrt, TMath::Max

GenJetSmearer::GenJetSmearer()
  : coeff_(1.00) // CV: take jet pT resolution to be 100%*sqrt(pT)
{
  rnd_.SetSeed(0);
}

GenJetSmearer::~GenJetSmearer()
{}

void 
GenJetSmearer::set_coeff(double coeff)
{
  coeff_ = coeff;
}

double 
GenJetSmearer::get_coeff() const
{
  return coeff_;
}

GenJet
GenJetSmearer::operator()(const GenJet& jet) const
{
  double sigma = coeff_*TMath::Sqrt(TMath::Max(1., jet.pt()));
  double jetPt_smeared = rnd_.Gaus(jet.pt(), sigma);
  GenJet jet_smeared(jetPt_smeared, jet.eta(), jet.phi(), jet.mass(), jet.pdgId());
  return jet_smeared;
}
