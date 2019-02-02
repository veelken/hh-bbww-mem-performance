#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/GenMEt.h"

#include <TMath.h> // TMath::Sqrt, TMath::ATan2

GenMEt::GenMEt()
  : px_(0.)
  , py_(0.)
  , pt_(0.)
  , phi_(0.)
{}

GenMEt::GenMEt(Float_t px, Float_t py)
  : px_(px)
  , py_(py)
{
  pt_ = TMath::Sqrt(px_*px_ + py_*py_);
  phi_ = TMath::ATan2(py_, px_);
}

Double_t 
GenMEt::px() const
{
  return px_;
}

Double_t 
GenMEt::py() const
{
  return py_;
}
 
Double_t 
GenMEt::pt() const
{
  return pt_;
}

Double_t 
GenMEt::phi() const
{
  return phi_;
}

std::ostream& operator<<(std::ostream& stream, const GenMEt& met)
{
  stream << " pX = " << met.px() << ", pY = " << met.py() << " (pT = " << met.pt() << ", phi = " << met.phi() << ")\n";
  return stream;
}
