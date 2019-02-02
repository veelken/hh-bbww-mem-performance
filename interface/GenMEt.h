#ifndef tthAnalysis_HiggsToTauTau_GenMEt_h
#define tthAnalysis_HiggsToTauTau_GenMEt_h

#include <Rtypes.h> // Float_t

#include <ostream>

class GenMEt
{
 public:
  GenMEt();

  GenMEt(Float_t px, Float_t py);

  /**
   * @brief Funtions to access data-members
   * @return Values of data-members
   */  
  Double_t px() const;
  Double_t py() const;
  Double_t pt() const;
  Double_t phi() const;

 private:
  Float_t px_;
  Float_t py_;
  Float_t pt_;
  Float_t phi_;
};

std::ostream& operator<<(std::ostream& stream, const GenMEt& met);

#endif // hhAnalysis_bbwwMEMPerformanceStudies_GenMEt_h
