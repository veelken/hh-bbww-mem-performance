#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/memNtupleAuxFunctions_dilepton.h"

#include "hhAnalysis/bbwwMEM/interface/measuredParticleAuxFunctions.h" // findGenMatch

void
addGenMatches_dilepton(MEMEvent_dilepton& memEvent,
                       const std::vector<const GenJet*>& genBJets,
                       const std::vector<const GenLepton*>& genLeptons,
                       double genMEtPx, double genMEtPy)
{
  memEvent.set_genBJet1(mem::findGenMatch(memEvent.measuredBJet1(), genBJets));
  memEvent.set_genBJet2(mem::findGenMatch(memEvent.measuredBJet2(), genBJets));
  memEvent.set_genLepton1(mem::findGenMatch(memEvent.measuredLepton1(), genLeptons));
  memEvent.set_genLepton2(mem::findGenMatch(memEvent.measuredLepton2(), genLeptons));
  memEvent.set_genMEtPx(genMEtPx);
  memEvent.set_genMEtPy(genMEtPy);
}
