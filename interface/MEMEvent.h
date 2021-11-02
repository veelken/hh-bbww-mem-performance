#ifndef hhAnalysis_bbwwMEMPerformanceStudies_MEMEvent_h
#define hhAnalysis_bbwwMEMPerformanceStudies_MEMEvent_h

#include "tthAnalysis/HiggsToTauTau/interface/EventInfo.h" // EventInfo
#include "tthAnalysis/HiggsToTauTau/interface/GenJet.h"    // GenJet

#include "hhAnalysis/bbwwMEM/interface/MeasuredParticle.h" // MeasuredParticle
#include "hhAnalysis/bbwwMEM/interface/MEMResult.h"        // MEMResultBase

#include <TMatrixD.h> // TMatrixD

class MEMEvent
{
 public:
  MEMEvent(const EventInfo & eventInfo, bool isSignal,
           const mem::MeasuredParticle* measuredBJet1, const mem::MeasuredParticle* measuredBJet2, 
           double measuredMEtPx, double measuredMEtPy, const TMatrixD& measuredMEtCov);
  ~MEMEvent();

  void set_genBJet1(const GenJet* genBJet1);
  void set_genBJet2(const GenJet* genBJet2);

  void set_isBoosted_Hbb(bool isBoosted_Hbb);

  void set_genMEtPx(double genMEtPx);
  void set_genMEtPy(double genMEtPy);

  void set_memResult(const MEMResultBase& memResult);
  void set_memCpuTime(double memCpuTime);

  void set_barcode(int barcode);

  const EventInfo & eventInfo() const;
  bool isSignal() const;

  const mem::MeasuredParticle* measuredBJet1() const;
  const GenJet* genBJet1() const;
  const mem::MeasuredParticle* measuredBJet2() const;
  const GenJet* genBJet2() const;
  bool isBoosted_Hbb() const;
  int numMeasuredBJets() const;
  int numGenBJets() const;

  double measuredMEtPx() const;
  double genMEtPx() const;
  double measuredMEtPy() const;
  double genMEtPy() const;
  const TMatrixD & measuredMEtCov() const;

  const MEMResultBase & memResult() const;
  double memCpuTime() const;

  int barcode() const;

 private:
  void countMeasuredBJets();
  void countGenBJets();

  const EventInfo* eventInfo_;
  bool isSignal_;

  const mem::MeasuredParticle* measuredBJet1_;
  const GenJet* genBJet1_;
  const mem::MeasuredParticle* measuredBJet2_;
  const GenJet* genBJet2_;
  bool isBoosted_Hbb_;
  int numMeasuredBJets_;
  int numGenBJets_;

  double measuredMEtPx_;
  double genMEtPx_;
  double measuredMEtPy_;
  double genMEtPy_;
  TMatrixD measuredMEtCov_;

  MEMResultBase memResult_;
  double memCpuTime_;

  mutable int barcode_;
};

#endif // hhAnalysis_bbwwMEMPerformanceStudies_MEMEvent_h

