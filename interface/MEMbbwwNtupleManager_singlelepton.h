#ifndef hhAnalysis_bbwwMEMPerformanceStudies_MEMbbwwNtupleManager_singlelepton_h
#define hhAnalysis_bbwwMEMPerformanceStudies_MEMbbwwNtupleManager_singlelepton_h

#include "tthAnalysis/HiggsToTauTau/interface/GenJet.h"                          // GenJet
#include "tthAnalysis/HiggsToTauTau/interface/GenLepton.h"                       // GenLepton

#include "hhAnalysis/bbwwMEM/interface/MeasuredParticle.h"                       // MeasuredParticle
#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/MEMbbwwNtupleManager.h" // MEMbbwwNtupleManager

#include <TMatrixD.h>                                                            // TMatrixD

class MEMbbwwNtupleManager_singlelepton : public MEMbbwwNtupleManager
{
 public:
  MEMbbwwNtupleManager_singlelepton(const std::string & outputTreeName);
  ~MEMbbwwNtupleManager_singlelepton();

  void initializeBranches();
  using MEMbbwwNtupleManager::read;
  void read(const std::vector<mem::MeasuredParticle> & measuredParticles,
            double measuredMEtPx, double measuredMEtPy, const TMatrixD & measuredMEtCov);
  void read(const std::vector<const GenJet *> & genBJets,
            const std::vector<const GenJet *> & genWJets,
            const GenLepton * genLepton,
            double genMEtPx, double genMEtPy);
  void resetBranches();

 protected:
  jetBranches wjet1_;
  jetBranches wjet2_;
  Int_t nwjets_;
  jetBranches gen_wjet1_;
  jetBranches gen_wjet2_;
  Int_t gen_nwjets_;

  leptonBranches lepton_;
  Int_t nleptons_;
  leptonBranches gen_lepton_;
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
