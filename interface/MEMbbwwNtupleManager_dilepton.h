#ifndef hhAnalysis_bbwwMEMPerformanceStudies_MEMbbwwNtupleManager_dilepton_h
#define hhAnalysis_bbwwMEMPerformanceStudies_MEMbbwwNtupleManager_dilepton_h

#include "tthAnalysis/HiggsToTauTau/interface/GenJet.h"                          // GenJet
#include "tthAnalysis/HiggsToTauTau/interface/GenLepton.h"                       // GenLepton

#include "hhAnalysis/bbwwMEM/interface/MeasuredParticle.h"                       // MeasuredParticle
#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/MEMbbwwNtupleManager.h" // MEMbbwwNtupleManager

#include <TMatrixD.h>                                                            // TMatrixD

class MEMbbwwNtupleManager_dilepton : public MEMbbwwNtupleManager
{
 public:
  MEMbbwwNtupleManager_dilepton(const std::string & outputTreeName);
  ~MEMbbwwNtupleManager_dilepton();

  void initializeBranches();
  void resetBranches();
  void read(const std::vector<mem::MeasuredParticle> & measuredParticles,
            double measuredMEtPx, double measuredMEtPy, const TMatrixD & measuredMEtCov);
  void read(const std::vector<GenJet *> & genBJets,
            const std::vector<GenLepton *> & genLeptons,
            double genMEtPx, double genMEtPy);

 protected:
  leptonBranches lepton1_;
  leptonBranches lepton2_;
  Int_t nleptons_;
  leptonBranches gen_lepton1_;
  leptonBranches gen_lepton2_;
  Int_t gen_nleptons_;
  
  // CV: define auxiliary variables for BDT regression training
  Float_t ptww_;
  Float_t mww_;
  Float_t ptll_;
  Float_t mll_;
  Float_t ptmiss_;
};

#endif // hhAnalysis_bbwwMEMPerformanceStudies_MEMbbwwNtupleManager_dilepton_h
