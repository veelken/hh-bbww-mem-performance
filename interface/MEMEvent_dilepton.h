#ifndef hhAnalysis_bbwwMEMPerformanceStudies_MEMEvent_dilepton_h
#define hhAnalysis_bbwwMEMPerformanceStudies_MEMEvent_dilepton_h

#include "tthAnalysis/HiggsToTauTau/interface/GenLepton.h" // GenLepton

#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/MEMEvent.h"

class MEMEvent_dilepton : public MEMEvent
{
 public:
  MEMEvent_dilepton(const EventInfo & eventInfo, bool isSignal,
                    const mem::MeasuredParticle* measuredBJet1, const mem::MeasuredParticle* measuredBJet2, 
                    const mem::MeasuredParticle* measuredLepton1, const mem::MeasuredParticle* measuredLepton2, 
                    double measuredMEtPx, double measuredMEtPy, const TMatrixD& measuredMEtCov);
  ~MEMEvent_dilepton();

  void set_genLepton1(const GenLepton* genLepton1);
  void set_genLepton2(const GenLepton* genLepton2);

  const mem::MeasuredParticle* measuredLepton1() const;
  const GenLepton* genLepton1() const;
  const mem::MeasuredParticle* measuredLepton2() const;
  const GenLepton* genLepton2() const;
  int numMeasuredLeptons() const;
  int numGenLeptons() const;

 private:
  void countMeasuredLeptons();
  void countGenLeptons();

  const mem::MeasuredParticle* measuredLepton1_;
  const GenLepton* genLepton1_;
  const mem::MeasuredParticle* measuredLepton2_;
  const GenLepton* genLepton2_; 
  int numMeasuredLeptons_;
  int numGenLeptons_;
};

#endif // hhAnalysis_bbwwMEMPerformanceStudies_MEMEvent_dilepton_h

