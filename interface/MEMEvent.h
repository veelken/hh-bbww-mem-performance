#ifndef hhAnalysis_bbwwMEMPerformanceStudies_MEMEvent_h
#define hhAnalysis_bbwwMEMPerformanceStudies_MEMEvent_h

#include "tthAnalysis/HiggsToTauTau/interface/EventInfo.h" // EventInfo
#include "tthAnalysis/HiggsToTauTau/interface/GenJet.h"    // GenJet

#include "hhAnalysis/bbwwMEM/interface/MeasuredParticle.h" // MeasuredParticle
#include "hhAnalysis/bbwwMEM/interface/MEMResult.h"        // MEMResultBase

#include <TMatrixD.h> // TMatrixD

class MEMEventInfo
{
 public:
  MEMEventInfo(UInt_t run, UInt_t lumi, ULong64_t event, Float_t genWeight)
    : run_(run)
    , lumi_(lumi)
    , event_(event)
    , genWeight_(genWeight)
  {}
  ~MEMEventInfo()
  {}

  UInt_t    run()       const { return run_;       }
  UInt_t    lumi()      const { return lumi_;      }
  ULong64_t event()     const { return event_;     }
  Float_t   genWeight() const { return genWeight_; }

 private:
  UInt_t    run_;       ///< run number
  UInt_t    lumi_;      ///< luminosity
  ULong64_t event_;     ///< event number
  Float_t   genWeight_; ///< generator-level weight (only if MC)
};

class MEMEvent
{
 public:
  MEMEvent(const MEMEventInfo & eventInfo, bool isSignal,
           const mem::MeasuredParticle* measuredBJet1, const mem::MeasuredParticle* measuredBJet2,
           double measuredMEtPx, double measuredMEtPy, const TMatrixD& measuredMEtCov);
  ~MEMEvent();

  void set_genBJet1(const GenJet* genBJet1);
  void set_genBJet2(const GenJet* genBJet2);
  void set_numGenBJets(int numGenBJets);

  void set_isBoosted_Hbb(bool isBoosted_Hbb);

  void set_genMEtPx(double genMEtPx);
  void set_genMEtPy(double genMEtPy);

  void set_numMeasuredBJets_loose(int numMeasuredBJets_loose);
  void set_numMeasuredBJets_medium(int numMeasuredBJets_medium);

  void set_memResult(const MEMResultBase& memResult);
  void set_memCpuTime(double memCpuTime);

  void set_barcode(int barcode);

  const MEMEventInfo & eventInfo() const;
  bool isSignal() const;

  const mem::MeasuredParticle* measuredBJet1() const;
  const GenJet* genBJet1() const;
  const mem::MeasuredParticle* measuredBJet2() const;
  const GenJet* genBJet2() const;
  bool isBoosted_Hbb() const;
  int numMeasuredBJets() const;
  int numGenBJets() const;

  int numMeasuredBJets_loose() const;
  int numMeasuredBJets_medium() const;

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

  MEMEventInfo eventInfo_;
  bool isSignal_;

  const mem::MeasuredParticle* measuredBJet1_;
  const GenJet* genBJet1_;
  const mem::MeasuredParticle* measuredBJet2_;
  const GenJet* genBJet2_;
  bool isBoosted_Hbb_;
  int numMeasuredBJets_;
  int numGenBJets_;

  int numMeasuredBJets_loose_;
  int numMeasuredBJets_medium_;

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

