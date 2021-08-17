#ifndef hhAnalysis_bbwwMEMPerformanceStudies_HadWJetTF_toy_h
#define hhAnalysis_bbwwMEMPerformanceStudies_HadWJetTF_toy_h

#include "hhAnalysis/bbwwMEM/interface/HadWJetTF.h" // mem::HadWJetTF
#include "hhAnalysis/bbwwMEM/interface/memAuxFunctions.h" // mem::LorentzVector

namespace mem
{
  
class HadWJetTF_toy : public HadWJetTF
{
 public:
  HadWJetTF_toy(int = 0);
  ~HadWJetTF_toy();

  /// set resolution on jet pT
  void set_coeff(double coeff);
  
  /// set measured energy, pT and pseudo-rapidity of jet from W->jj decay
  void setInputs(const mem::LorentzVector&);

  /// evaluate transfer function (TF)
  double Eval(double) const;

 protected:  
  /// measured pT of jet from W->jj decay
  double measuredPt_;

  /// jet pT resolution parameter
  double coeff_;
};

}

#endif // hhAnalysis_bbwwMEMPerformanceStudies_HadWJetTF_toy_h
