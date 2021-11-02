#ifndef hhAnalysis_bbwwMEMPerformanceStudies_memNtupleAuxFunctions_singlelepton_h
#define hhAnalysis_bbwwMEMPerformanceStudies_memNtupleAuxFunctions_singlelepton_h

#include "tthAnalysis/HiggsToTauTau/interface/GenJet.h"     // GenJet
#include "tthAnalysis/HiggsToTauTau/interface/GenLepton.h"  // GenLepton

#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/MEMEvent_singlelepton.h" // MEMEvent_singlelepton

void
addGenMatches_singlelepton(MEMEvent_singlelepton & memEvent,
                           const std::vector<const GenJet*> & genBJets,
                           const std::vector<const GenJet*> & genWJets,
                           const std::vector<const GenLepton*> & genLeptons,
                           double genMEtPx, double genMEtPy);

#endif // hhAnalysis_bbwwMEMPerformanceStudies_memNtupleAuxFunctions_singlelepton_h


