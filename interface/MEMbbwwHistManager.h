#ifndef hhAnalysis_bbwwMEMPerformanceStudies_MEMbbwwHistManager_h
#define hhAnalysis_bbwwMEMPerformanceStudies_MEMbbwwHistManager_h

#include "tthAnalysis/HiggsToTauTau/interface/HistManagerBase.h" // HistManagerBase
#include "tthAnalysis/HiggsToTauTau/interface/histogramAuxFunctions.h" // fillWithOverFlow()

#include "hhAnalysis/bbwwMEM/interface/MEMResult.h" // MEMbbwwResultDilepton, MEMbbwwResultSingleLepton

template <class T>
class MEMbbwwHistManager
  : public HistManagerBase
{
 public:
  MEMbbwwHistManager(const edm::ParameterSet& cfg)
    : HistManagerBase(cfg)
  {
    central_or_shiftOptions_["log_memProb_signal"] = { "central" };
    central_or_shiftOptions_["log_memProbErr_signal"] = { "central" };
    central_or_shiftOptions_["log_memProb_background"] = { "central" };
    central_or_shiftOptions_["log_memProbErr_background"] = { "central" };
    central_or_shiftOptions_["memLR"] = { "*" };
    central_or_shiftOptions_["log_memLR_div_Err"] = { "central" };
    central_or_shiftOptions_["memScore"] = { "*" };
    central_or_shiftOptions_["memCpuTime"] = { "central" };
    central_or_shiftOptions_["EventCounter"] = { "*" };
  }
  ~MEMbbwwHistManager() {}

  /// book and fill histograms
  void
  bookHistograms(TFileDirectory& dir) override
  {
    histogram_log_memProb_signal_        = book1D(dir, "log_memProb_signal",        2000, -100., +100.);
    histogram_log_memProbErr_signal_     = book1D(dir, "log_memProbErr_signal",     2000, -100., +100.);
    histogram_log_memProb_background_    = book1D(dir, "log_memProb_background",    2000, -100., +100.);
    histogram_log_memProbErr_background_ = book1D(dir, "log_memProbErr_background", 2000, -100., +100.);
    histogram_memLR_                     = book1D(dir, "memLR",                     3600,    0.,    1.);
    histogram_log_memLR_div_Err_         = book1D(dir, "log_memLR_div_Err",         2000,  -10.,  +10.);
    histogram_memScore_                  = book1D(dir, "memScore",                  3600,  -18.,  +18.);
    histogram_memCpuTime_                = book1D(dir, "memCpuTime",                1000,    0., 1000.);

    histogram_EventCounter_              = book1D(dir, "EventCounter",                 1,   -0.5,  +0.5);
  }

  void
  fillHistograms(const T& memResult, double memCpuTime, double evtWeight)
  {
    const double evtWeightErr = 0.;

    fillWithOverFlow_logx(histogram_log_memProb_signal_,        memResult.getProb_signal(),        evtWeight, evtWeightErr);            
    fillWithOverFlow_logx(histogram_log_memProbErr_signal_,     memResult.getProbErr_signal(),     evtWeight, evtWeightErr);    
    fillWithOverFlow_logx(histogram_log_memProb_background_,    memResult.getProb_background(),    evtWeight, evtWeightErr);    
    fillWithOverFlow_logx(histogram_log_memProbErr_background_, memResult.getProbErr_background(), evtWeight, evtWeightErr);    
    fillWithOverFlow(histogram_memLR_,                          memResult.getLikelihoodRatio(),    evtWeight, evtWeightErr);           
    if ( memResult.getLikelihoodRatioErr() > 0. ) {
      double memLR_div_Err = memResult.getLikelihoodRatio()/memResult.getLikelihoodRatioErr();
      fillWithOverFlow_logx(histogram_log_memLR_div_Err_,       memLR_div_Err,                     evtWeight, evtWeightErr);  
    }
    fillWithOverFlow(histogram_memScore_,                       memResult.getScore(),              evtWeight, evtWeightErr);  
    fillWithOverFlow(histogram_memCpuTime_,                     memCpuTime,                        evtWeight, evtWeightErr);            

    fillWithOverFlow(histogram_EventCounter_,                   0.,                                evtWeight, evtWeightErr);
  }
  
  const TH1*
  getHistogram_EventCounter() const
  {
    return histogram_EventCounter_;
  }

 private:
  void fillWithOverFlow_logx(TH1* histogram, double x, double evtWeight, double evtWeightErr = 0.)
  {    
    const double nonzero = 1.e-30;
    fillWithOverFlow(histogram, TMath::Log(TMath::Max(nonzero, x)), evtWeight, evtWeightErr);
  }

  TH1* histogram_log_memProb_signal_;
  TH1* histogram_log_memProbErr_signal_;
  TH1* histogram_log_memProb_background_;
  TH1* histogram_log_memProbErr_background_;
  TH1* histogram_memLR_;
  TH1* histogram_log_memLR_div_Err_;
  TH1* histogram_memScore_;
  TH1* histogram_memCpuTime_;

  TH1* histogram_EventCounter_;
};

typedef MEMbbwwHistManager<MEMbbwwResultDilepton> MEMbbwwHistManagerDilepton;
typedef MEMbbwwHistManager<MEMbbwwResultSingleLepton> MEMbbwwHistManagerSingleLepton;

#endif // hhAnalysis_bbwwMEMPerformanceStudies_MEMbbwwHistManager_h
