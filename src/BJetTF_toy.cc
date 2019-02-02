#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/BJetTF_toy.h"

#include "tthAnalysis/tthMEM/interface/JetTransferFunction.h" // gaussianPDF
#include "hhAnalysis/bbwwMEM/interface/memAuxFunctions.h" // mem::bottomQuarkMass2

#include <TMath.h> // TMath::Sqrt, TMath::CosH

using namespace mem;

BJetTF_toy::BJetTF_toy(int verbosity)
  : BJetTF(verbosity)
  , coeff_(1.00) // CV: take jet pT resolution to be 100%/sqrt(pT)
{}

BJetTF_toy::~BJetTF_toy()
{}

void BJetTF_toy::set_coeff(double coeff)
{
  coeff_ = coeff;
}

void BJetTF_toy::setInputs(const mem::LorentzVector& measuredP4)
{
  BJetTF::setInputs(measuredP4);
  measuredPt_ = measuredP4.pt();
}

double BJetTF_toy::Eval(double trueEn) const
{
  // CV: formulae taken from https://en.wikipedia.org/wiki/Pseudorapidity
  double prob = 0.;
  if ( trueEn > bottomQuarkMass ) {
    double trueP = TMath::Sqrt(trueEn*trueEn - bottomQuarkMass2);
    double truePt = trueP/TMath::CosH(measuredEta_);
    double sigma = coeff_/TMath::Sqrt(TMath::Max(1., truePt));
    prob = tthMEM::functions::gaussianPDF(measuredPt_, truePt, sigma);
  }
  return prob;
}
