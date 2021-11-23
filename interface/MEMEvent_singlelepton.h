#ifndef hhAnalysis_bbwwMEMPerformanceStudies_MEMEvent_dilepton_h
#define hhAnalysis_bbwwMEMPerformanceStudies_MEMEvent_dilepton_h

#include "tthAnalysis/HiggsToTauTau/interface/GenLepton.h"           // GenLepton

#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/MEMEvent.h" // MEMEvent, MEMEventInfo

class MEMEvent_singlelepton : public MEMEvent
{
 public:
  MEMEvent_singlelepton(const MEMEventInfo & eventInfo, bool isSignal,
                        const mem::MeasuredParticle* measuredBJet1, const mem::MeasuredParticle* measuredBJet2, 
                        const mem::MeasuredParticle* measuredWJet1, const mem::MeasuredParticle* measuredWJet2, 
                        const mem::MeasuredParticle* measuredLepton,
                        double measuredMEtPx, double measuredMEtPy, const TMatrixD& measuredMEtCov);
  ~MEMEvent_singlelepton();

  void set_genWJet1(const GenJet* genWJet1);
  void set_genWJet2(const GenJet* genWJet2);
  void set_numGenWJets(int numGenWJets);
  void set_isBoosted_Wjj(bool isBoosted_Wjj);
  
  void set_genLepton(const GenLepton* genLepton);
  void set_numGenLeptons(int numGenLeptons);

  const mem::MeasuredParticle* measuredWJet1() const;
  const GenJet* genWJet1() const;
  const mem::MeasuredParticle* measuredWJet2() const;
  const GenJet* genWJet2() const;
  bool isBoosted_Wjj() const;
  int numMeasuredWJets() const;
  int numGenWJets() const;

  const mem::MeasuredParticle* measuredLepton() const;
  const GenLepton* genLepton() const;
  int numMeasuredLeptons() const;
  int numGenLeptons() const;

 private:
  void countMeasuredWJets();
  void countGenWJets();

  void countMeasuredLeptons();
  void countGenLeptons();
 
  const mem::MeasuredParticle* measuredWJet1_;
  const GenJet* genWJet1_;
  const mem::MeasuredParticle* measuredWJet2_;
  const GenJet* genWJet2_;
  bool isBoosted_Wjj_;
  int numMeasuredWJets_;
  int numGenWJets_;

  const mem::MeasuredParticle* measuredLepton_;
  const GenLepton* genLepton_;
  int numMeasuredLeptons_;
  int numGenLeptons_;
};

#endif // hhAnalysis_bbwwMEMPerformanceStudies_MEMEvent_dilepton_h

