#ifndef hhAnalysis_bbwwMEMPerformanceStudies_memNtupleAuxFunctions_dilepton_h
#define hhAnalysis_bbwwMEMPerformanceStudies_memNtupleAuxFunctions_dilepton_h

#include "tthAnalysis/HiggsToTauTau/interface/GenJet.h"     // GenJet
#include "tthAnalysis/HiggsToTauTau/interface/GenLepton.h"  // GenLepton

#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/MEMEvent_dilepton.h" // MEMEvent_dilepton

void
addGenMatches_dilepton(MEMEvent_dilepton & memEvent,
                       const std::vector<const GenJet*> & genBJets,
                       const std::vector<const GenLepton*> & genLeptons,
                       double genMEtPx, double genMEtPy);

#endif // hhAnalysis_bbwwMEMPerformanceStudies_MEMEvent_h


