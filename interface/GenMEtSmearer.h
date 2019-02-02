#ifndef tthAnalysis_HiggsToTauTau_GenMEtSmearer_h
#define tthAnalysis_HiggsToTauTau_GenMEtSmearer_h

#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/GenMEt.h" // GenMEt

#include <TRandom3.h> // TRandom3

class GenMEtSmearer
{
 public:
  GenMEtSmearer();
  ~GenMEtSmearer();

  /**
   * @brief Set resolutions
   */
  void set_sigmaX(double sigmaX);
  void set_sigmaY(double sigmaY);

  /**
   * @brief Get resolutions
   */
  double get_sigmaX() const;
  double get_sigmaY() const;

  /**
   * @brief Smear generator-level missing transverse momentum (MET)
   * @return Smeared MET
   */
  GenMEt operator()(const GenMEt&) const;

 protected:
  double sigmaX_;
  double sigmaY_;
  mutable TRandom3 rnd_;
};

#endif // hhAnalysis_bbwwMEMPerformanceStudies_GenMEtSmearer_h
