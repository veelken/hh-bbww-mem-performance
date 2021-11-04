#ifndef hhAnalysis_bbwwMEMPerformanceStudies_MEMbbwwNtupleManager_singlelepton_h
#define hhAnalysis_bbwwMEMPerformanceStudies_MEMbbwwNtupleManager_singlelepton_h

#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/MEMbbwwNtupleManager.h"  // MEMbbwwNtupleManager
#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/MEMEvent_singlelepton.h" // MEMEvent_singlelepton

class MEMbbwwNtupleManager_singlelepton : public MEMbbwwNtupleManager
{
 public:
  MEMbbwwNtupleManager_singlelepton(const std::string & outputDirectoryName, const std::string & outputTreeName);
  ~MEMbbwwNtupleManager_singlelepton();

  void initializeBranches();
  void read(const MEMEvent_singlelepton & memEvent);
  void resetBranches();

 protected:
  measuredJetBranches wjet1_;
  measuredJetBranches wjet2_;
  Int_t nwjets_;
  genJetBranches gen_wjet1_;
  genJetBranches gen_wjet2_;
  Int_t gen_nwjets_;

  measuredLeptonBranches lepton_;
  Int_t nleptons_;
  genLeptonBranches gen_lepton_;
  Int_t gen_nleptons_;

  // CV: define auxiliary variables for BDT regression training
  Float_t ptjj_;
  Float_t drjj_;
  Float_t mjj_;
  Float_t ptww_;
  Float_t mww_;
  Float_t mt_;
  Float_t ptmiss_;
};

#endif // hhAnalysis_bbwwMEMPerformanceStudies_MEMbbwwNtupleManager_singlelepton_h
