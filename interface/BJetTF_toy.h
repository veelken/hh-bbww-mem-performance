#ifndef hhAnalysis_bbwwMEMPerformanceStudies_BJetTF_toy_h
#define hhAnalysis_bbwwMEMPerformanceStudies_BJetTF_toy_h

#include "hhAnalysis/bbwwMEM/interface/BJetTF.h" // mem::BJetTF
#include "hhAnalysis/bbwwMEM/interface/memAuxFunctions.h" // mem::LorentzVector

namespace mem
{
  
class BJetTF_toy : public BJetTF
{
 public:
  BJetTF_toy(int = 0);
  ~BJetTF_toy();

  /// set resolution on jet pT
  void set_coeff(double coeff);
  
  /// set measured b-jet energy, pT and pseudo-rapidity
  void setInputs(const mem::LorentzVector&);

  /// evaluate transfer function (TF)
  double Eval(double) const;

 protected:  
  /// measured b-jet pT
  double measuredPt_;

  /// jet pT resolution parameter
  double coeff_;
};

}

#endif // hhAnalysis_bbwwMEMPerformanceStudies_BJetTF_toy_h
