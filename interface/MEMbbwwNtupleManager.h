#ifndef hhAnalysis_bbwwMEMPerformanceStudies_MEMbbwwNtupleManager_h
#define hhAnalysis_bbwwMEMPerformanceStudies_MEMbbwwNtupleManager_h

#include "CommonTools/Utils/interface/TFileDirectory.h"                // TFileDirectory
#include "DataFormats/Math/interface/deltaR.h"                         // deltaR

#include "tthAnalysis/HiggsToTauTau/interface/EventInfo.h"             // EventInfo
#include "tthAnalysis/HiggsToTauTau/interface/GenJet.h"                // GenJet
#include "tthAnalysis/HiggsToTauTau/interface/GenLepton.h"             // GenLepton
#include "tthAnalysis/HiggsToTauTau/interface/TypeTraits.h"            // Traits<>

#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/MEMEvent.h"   // MEMEvent

#include <TTree.h>                                                     // TTree
#include <TMatrixD.h>                                                  // TMatrixD

class MEMbbwwNtupleManager
{
public:
  MEMbbwwNtupleManager(const std::string & outputDirectoryName, const std::string & outputTreeName);
  virtual ~MEMbbwwNtupleManager();

  void makeTree(TFileDirectory & dir);

  virtual void initializeBranches();
  virtual void read(const MEMEvent & memEvent);
  void fill();
  virtual void resetBranches();

protected:

  std::string outputTreeName_;
  TTree* tree_;

  Int_t run_;
  Int_t ls_;
  Long64_t event_;

  Float_t genWeight_;

  Bool_t isSignal_;

  Double_t memProbS_;
  Double_t memProbSerr_;
  Double_t memProbB_;
  Double_t memProbBerr_;
  Double_t memLR_;
  Double_t memLRerr_;
  Float_t  memCpuTime_;

  struct genJetBranches
  {
    genJetBranches(const std::string & branchName)
      : branchName_(branchName)
    {
      resetBranches();
    }
    virtual ~genJetBranches()
    {}
    virtual void initializeBranches(TTree * tree)
    {
      tree->Branch(Form("%s_pt",   branchName_.c_str()), &pt_,   Form("%s_pt/%s",   branchName_.c_str(), Traits<Float_t>::TYPE_NAME));
      tree->Branch(Form("%s_eta",  branchName_.c_str()), &eta_,  Form("%s_eta/%s",  branchName_.c_str(), Traits<Float_t>::TYPE_NAME));
      tree->Branch(Form("%s_phi",  branchName_.c_str()), &phi_,  Form("%s_phi/%s",  branchName_.c_str(), Traits<Float_t>::TYPE_NAME));
      tree->Branch(Form("%s_mass", branchName_.c_str()), &mass_, Form("%s_mass/%s", branchName_.c_str(), Traits<Float_t>::TYPE_NAME));
    }
    virtual void resetBranches()
    {
      pt_   = 0.;
      eta_  = 0.;
      phi_  = 0.;
      mass_ = 0.;
    }
    virtual void read(const GenJet * jet)
    {
      if ( jet )
      {
        pt_   = jet->p4().pt();
        eta_  = jet->p4().eta();
        phi_  = jet->p4().phi();
        mass_ = jet->p4().mass();
      }
    }
    std::string branchName_;
    Float_t pt_;
    Float_t eta_;
    Float_t phi_;
    Float_t mass_;
  };
  struct measuredJetBranches : public genJetBranches
  {
    measuredJetBranches(const std::string & branchName)
      : genJetBranches(branchName)
    {
      resetBranches();
    }
    ~measuredJetBranches()
    {}
    void initializeBranches(TTree * tree)
    {
      genJetBranches::initializeBranches(tree);
      tree->Branch(Form("%s_isGenMatched", branchName_.c_str()), &isGenMatched_, Form("%s_isGenMatched/%s", branchName_.c_str(), Traits<Bool_t>::TYPE_NAME));
    }
    void resetBranches()
    {
      genJetBranches::resetBranches();
      isGenMatched_ = false;
    }
    void read(const mem::MeasuredParticle * jet, bool isGenMatched)
    {
      if ( jet )
      {
        pt_           = jet->pt();
        eta_          = jet->eta();
        phi_          = jet->phi();
        mass_         = jet->mass();
        isGenMatched_ = isGenMatched;
      }
    }
    bool isGenMatched_;
  };

  measuredJetBranches bjet1_;
  measuredJetBranches bjet2_;
  Int_t nbjets_;
  genJetBranches gen_bjet1_;
  genJetBranches gen_bjet2_;
  Int_t gen_nbjets_;

  struct genLeptonBranches
  {
    genLeptonBranches(const std::string & branchName)
      : branchName_(branchName)
    {
      resetBranches();
    }
    virtual ~genLeptonBranches()
    {}
    virtual void initializeBranches(TTree * tree)
    {
      tree->Branch(Form("%s_pt",    branchName_.c_str()), &pt_,    Form("%s_pt/%s",    branchName_.c_str(), Traits<Float_t>::TYPE_NAME));
      tree->Branch(Form("%s_eta",   branchName_.c_str()), &eta_,   Form("%s_eta/%s",   branchName_.c_str(), Traits<Float_t>::TYPE_NAME));
      tree->Branch(Form("%s_phi",   branchName_.c_str()), &phi_,   Form("%s_phi/%s",   branchName_.c_str(), Traits<Float_t>::TYPE_NAME));
      tree->Branch(Form("%s_pdgId", branchName_.c_str()), &pdgId_, Form("%s_pdgId/%s", branchName_.c_str(), Traits<Int_t>::TYPE_NAME));
    }
    virtual void resetBranches()
    {
      pt_    = 0.;
      eta_   = 0.;
      phi_   = 0.;
      pdgId_ = 0;
    }
    virtual void read(const GenLepton * lepton)
    {
      if ( lepton )
      {
        pt_    = lepton->p4().pt();
        eta_   = lepton->p4().eta();
        phi_   = lepton->p4().phi();
        pdgId_ = lepton->pdgId();
      }
    }
    std::string branchName_;
    Float_t pt_;
    Float_t eta_;
    Float_t phi_;
    Int_t pdgId_;
  };
  struct measuredLeptonBranches : public genLeptonBranches
  {
    measuredLeptonBranches(const std::string & branchName)
      : genLeptonBranches(branchName)
    {
      resetBranches();
    }
    ~measuredLeptonBranches()
    {}
    void initializeBranches(TTree * tree)
    {
      genLeptonBranches::initializeBranches(tree);
      tree->Branch(Form("%s_isGenMatched", branchName_.c_str()), &isGenMatched_, Form("%s_isGenMatched/%s", branchName_.c_str(), Traits<Bool_t>::TYPE_NAME));
    }
    void resetBranches()
    {
      genLeptonBranches::resetBranches();
      isGenMatched_ = false;
    }
    void read(const mem::MeasuredParticle * lepton, bool isGenMatched)
    {
      if ( lepton )
      {
        pt_  = lepton->pt();
        eta_ = lepton->eta();
        phi_ = lepton->phi();
        if      ( lepton->type() == mem::MeasuredParticle::kElectron && lepton->charge() < 0 ) pdgId_ = +11;
        else if ( lepton->type() == mem::MeasuredParticle::kElectron && lepton->charge() > 0 ) pdgId_ = -11;
        else if ( lepton->type() == mem::MeasuredParticle::kMuon     && lepton->charge() < 0 ) pdgId_ = +13;
        else if ( lepton->type() == mem::MeasuredParticle::kMuon     && lepton->charge() > 0 ) pdgId_ = -13;
        else assert(0);
        isGenMatched_ = isGenMatched;
      }
    }
    bool isGenMatched_;
  };

  struct metBranches
  {
    metBranches(const std::string & branchName)
      : branchName_(branchName)
    {}
    ~metBranches()
    {}
    void initializeBranches(TTree * tree)
    {
      tree->Branch(Form("%s_px",    branchName_.c_str()), &px_,    Form("%s_px/%s",    branchName_.c_str(), Traits<Float_t>::TYPE_NAME));
      tree->Branch(Form("%s_py",    branchName_.c_str()), &py_,    Form("%s_py/%s",    branchName_.c_str(), Traits<Float_t>::TYPE_NAME));
      tree->Branch(Form("%s_cov00", branchName_.c_str()), &cov00_, Form("%s_cov00/%s", branchName_.c_str(), Traits<Float_t>::TYPE_NAME));
      tree->Branch(Form("%s_cov01", branchName_.c_str()), &cov01_, Form("%s_cov01/%s", branchName_.c_str(), Traits<Float_t>::TYPE_NAME));
      tree->Branch(Form("%s_cov11", branchName_.c_str()), &cov11_, Form("%s_cov11/%s", branchName_.c_str(), Traits<Float_t>::TYPE_NAME));
    }
    void resetBranches()
    {
      px_    = 0.;
      py_    = 0.;
      cov00_ = 0.;
      cov01_ = 0.;
      cov11_ = 0.;
    }
    void read(double metPx, double metPy, const TMatrixD * metCov = nullptr)
    {
      px_    = metPx;
      py_    = metPy;
      if ( metCov )
      {
        cov00_ = metCov->operator()(0, 0);
        cov01_ = metCov->operator()(0, 1);
        cov11_ = metCov->operator()(1, 1);
      }
    }
    std::string branchName_;
    Float_t px_;
    Float_t py_;
    Float_t cov00_;
    Float_t cov01_;
    Float_t cov11_;
  };

  metBranches met_;
  metBranches gen_met_;

  int barcode_;

  // CV: define auxiliary variables for BDT regression training
  Float_t ptbb_;
  Float_t drbb_;
  Float_t mbb_;
};

#endif // hhAnalysis_bbwwMEMPerformanceStudies_MEMbbwwNtupleManager_h
