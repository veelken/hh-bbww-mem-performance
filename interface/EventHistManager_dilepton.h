#ifndef hhAnalysis_bbwwMEMPerformanceStudies_EventHistManager_dilepton_h
#define hhAnalysis_bbwwMEMPerformanceStudies_EventHistManager_dilepton_h

#include "tthAnalysis/HiggsToTauTau/interface/HistManagerBase.h"       // HistManagerBase
#include "tthAnalysis/HiggsToTauTau/interface/histogramAuxFunctions.h" // fillWithOverFlow()

class EventHistManager_dilepton
  : public HistManagerBase
{
 public:
  EventHistManager_dilepton(const edm::ParameterSet& cfg);
  ~EventHistManager_dilepton();

  /// book and fill histograms
  void
  bookHistograms(TFileDirectory& dir) override;

  void
  fillHistograms(double mbb, double mll, double evtWeight);
  
  const TH1*
  getHistogram_EventCounter() const
  {
    return histogram_EventCounter_;
  }

 private:
  TH1* histogram_mbb_;
  TH1* histogram_mll_;
  TH2* histogram_mll_vs_mbb_;

  TH1* histogram_EventCounter_;
};

#endif // hhAnalysis_bbwwMEMPerformanceStudies_EventHistManager_dilepton_h
