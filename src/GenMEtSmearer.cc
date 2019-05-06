#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/GenMEtSmearer.h"

GenMEtSmearer::GenMEtSmearer()
  : sigmaX_(10.)
  , sigmaY_(10.)
{
  rnd_.SetSeed(0);
}

GenMEtSmearer::~GenMEtSmearer()
{}

void 
GenMEtSmearer::set_sigmaX(double sigmaX)
{
  sigmaX_ = sigmaX;
}

void 
GenMEtSmearer::set_sigmaY(double sigmaY)
{
  sigmaY_ = sigmaY;
}

double 
GenMEtSmearer::get_sigmaX() const
{
  return sigmaX_;
}
 
double 
GenMEtSmearer::get_sigmaY() const
{
  return sigmaY_;
}

GenMEt 
GenMEtSmearer::operator()(const GenMEt& met) const
{
  double metPx_smeared = rnd_.Gaus(met.px(), sigmaX_);
  double metPy_smeared = rnd_.Gaus(met.py(), sigmaY_);
  GenMEt met_smeared(metPx_smeared, metPy_smeared);
  return met_smeared;
}
