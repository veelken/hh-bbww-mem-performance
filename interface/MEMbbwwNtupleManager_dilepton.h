#ifndef hhAnalysis_bbwwMEMPerformanceStudies_MEMbbwwNtupleManager_dilepton_h
#define hhAnalysis_bbwwMEMPerformanceStudies_MEMbbwwNtupleManager_dilepton_h

#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/MEMbbwwNtupleManager.h" // MEMbbwwNtupleManager
#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/MEMEvent_dilepton.h"    // MEMEvent_dilepton

class MEMbbwwNtupleManager_dilepton : public MEMbbwwNtupleManager
{
 public:
  MEMbbwwNtupleManager_dilepton(const std::string & outputTreeName);
  ~MEMbbwwNtupleManager_dilepton();

  void initializeBranches();
  void read(const MEMEvent_dilepton & memEvent);
  void resetBranches();

 protected:
  measuredLeptonBranches lepton1_;
  measuredLeptonBranches lepton2_;
  Int_t nleptons_;
  genLeptonBranches gen_lepton1_;
  genLeptonBranches gen_lepton2_;
  Int_t gen_nleptons_;
  
  // CV: define auxiliary variables for BDT regression training
  Float_t ptww_;
  Float_t mww_;
  Float_t ptll_;
  Float_t drll_;
  Float_t dphill_;
  Float_t mll_;
  Float_t ptmiss_;
};

#endif // hhAnalysis_bbwwMEMPerformanceStudies_MEMbbwwNtupleManager_dilepton_h
