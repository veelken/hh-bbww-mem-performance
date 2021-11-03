#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/EventHistManager_dilepton.h"

#include "tthAnalysis/HiggsToTauTau/interface/histogramAuxFunctions.h" // fillWithOverFlow, fillWithOverFlow2d

EventHistManager_dilepton::EventHistManager_dilepton(const edm::ParameterSet& cfg)
  : HistManagerBase(cfg)
{
  central_or_shiftOptions_["mbb"]          = { "central" };
  central_or_shiftOptions_["mll"]          = { "central" };
  central_or_shiftOptions_["mll_vs_mbb"]   = { "central" };
  central_or_shiftOptions_["EventCounter"] = { "*" };
}

EventHistManager_dilepton::~EventHistManager_dilepton()
{}

void
EventHistManager_dilepton::bookHistograms(TFileDirectory& dir)
{
  histogram_mbb_          = book1D(dir, "mbb",          200,  0.,  200.);
  histogram_mll_          = book1D(dir, "mll",          200,  0.,  200.);
  histogram_mll_vs_mbb_   = book2D(dir, "mll_vs_mbb",    40,  0.,  200.,  40,  0.,  200.);

  histogram_EventCounter_ = book1D(dir, "EventCounter",   1, -0.5, +0.5);
}

void
EventHistManager_dilepton::fillHistograms(double mbb, double mll, double evtWeight)
{
  const double evtWeightErr = 0.;

  fillWithOverFlow(histogram_mbb_,          mbb,      evtWeight, evtWeightErr);            
  fillWithOverFlow(histogram_mll_,               mll, evtWeight, evtWeightErr);    
  fillWithOverFlow2d(histogram_mll_vs_mbb_, mbb, mll, evtWeight, evtWeightErr);    

  fillWithOverFlow(histogram_EventCounter_,  0.,      evtWeight, evtWeightErr);
}
