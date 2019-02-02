#ifndef tthAnalysis_HiggsToTauTau_GenJetSmearer_h
#define tthAnalysis_HiggsToTauTau_GenJetSmearer_h

#include "tthAnalysis/HiggsToTauTau/interface/GenJet.h" // GenJet

#include <TRandom3.h> // TRandom3

class GenJetSmearer
{
 public:
  GenJetSmearer();
  ~GenJetSmearer();

  /**
   * @brief Set resolution
   */
  void set_coeff(double coeff);

  /**
   * @brief Get resolution
   */
  double get_coeff() const;

  /**
   * @brief Smear generator-level jet 
   * @return Smeared jet
   */

  GenJet operator()(const GenJet&) const;

 protected:
  double coeff_;
  mutable TRandom3 rnd_;
};

#endif // hhAnalysis_bbwwMEMPerformanceStudies_GenJetSmearer_h
