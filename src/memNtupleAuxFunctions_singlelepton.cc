#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/memNtupleAuxFunctions_singlelepton.h"

#include "hhAnalysis/bbwwMEM/interface/measuredParticleAuxFunctions.h" // findGenMatch

void
addGenMatches_singlelepton(MEMEvent_singlelepton & memEvent,
                           const std::vector<const GenJet*> & genBJets,
                           const std::vector<const GenJet*> & genWJets,
                           const std::vector<const GenLepton*> & genLeptons,
                           double genMEtPx, double genMEtPy)
{
  memEvent.set_genBJet1(mem::findGenMatch(memEvent.measuredBJet1(), genBJets));
  memEvent.set_genBJet2(mem::findGenMatch(memEvent.measuredBJet2(), genBJets));
  memEvent.set_genWJet1(mem::findGenMatch(memEvent.measuredWJet1(), genWJets));
  memEvent.set_genWJet2(mem::findGenMatch(memEvent.measuredWJet2(), genWJets));
  memEvent.set_genLepton(mem::findGenMatch(memEvent.measuredLepton(), genLeptons));
  memEvent.set_genMEtPx(genMEtPx);
  memEvent.set_genMEtPy(genMEtPy);
}
