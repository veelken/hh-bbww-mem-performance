#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/MEMEvent_singlelepton.h"

MEMEvent_singlelepton::MEMEvent_singlelepton(const EventInfo & eventInfo, bool isSignal,
                                             const mem::MeasuredParticle* measuredBJet1, const mem::MeasuredParticle* measuredBJet2, 
                                             const mem::MeasuredParticle* measuredWJet1, const mem::MeasuredParticle* measuredWJet2, 
                                             const mem::MeasuredParticle* measuredLepton,
                                             double measuredMEtPx, double measuredMEtPy, const TMatrixD& measuredMEtCov)
  : MEMEvent(eventInfo, isSignal, measuredBJet1, measuredBJet2, measuredMEtPx, measuredMEtPy, measuredMEtCov)
  , measuredWJet1_(measuredWJet1)
  , genWJet1_(nullptr)
  , measuredWJet2_(measuredWJet2)
  , genWJet2_(nullptr)
  , isBoosted_Wjj_(false)
  , numMeasuredWJets_(0)
  , numGenWJets_(0)
  , measuredLepton_(measuredLepton)
  , genLepton_(nullptr)
{}

MEMEvent_singlelepton::~MEMEvent_singlelepton()
{}

void 
MEMEvent_singlelepton::set_genWJet1(const GenJet* genWJet1)
{
  genWJet1_ = genWJet1;
  countGenWJets();
}

void 
MEMEvent_singlelepton::set_genWJet2(const GenJet* genWJet2)
{
  genWJet2_ = genWJet2;
  countGenWJets();
}

void 
MEMEvent_singlelepton::set_isBoosted_Wjj(bool isBoosted_Wjj)
{
  isBoosted_Wjj_ = isBoosted_Wjj;
}

void 
MEMEvent_singlelepton::set_genLepton(const GenLepton* genLepton)
{
  genLepton_ = genLepton;
}

const mem::MeasuredParticle* 
MEMEvent_singlelepton::measuredWJet1() const
{
  return measuredWJet1_;
}
  
const GenJet* 
MEMEvent_singlelepton::genWJet1() const
{
  return genWJet1_;
}
  
const mem::MeasuredParticle* 
MEMEvent_singlelepton::measuredWJet2() const
{
  return measuredWJet2_;
}

const GenJet* 
MEMEvent_singlelepton::genWJet2() const
{
  return genWJet2_;
}

bool 
MEMEvent_singlelepton::isBoosted_Wjj() const
{
  return isBoosted_Wjj_;
}

int 
MEMEvent_singlelepton::numMeasuredWJets() const
{
  return numMeasuredWJets_;
}
  
int 
MEMEvent_singlelepton::numGenWJets() const
{
  return numGenWJets_;
}

const mem::MeasuredParticle* 
MEMEvent_singlelepton::measuredLepton() const
{
  return measuredLepton_;
}

const GenLepton* 
MEMEvent_singlelepton::genLepton() const
{
  return genLepton_;
}

int 
MEMEvent_singlelepton::numMeasuredLeptons() const
{
  return numMeasuredLeptons_;
}
  
int 
MEMEvent_singlelepton::numGenLeptons() const
{
  return numGenLeptons_;
}

void 
MEMEvent_singlelepton::countMeasuredWJets()
{
  numMeasuredWJets_ = 0;
  if ( measuredWJet1_ ) ++numMeasuredWJets_;
  if ( measuredWJet2_ ) ++numMeasuredWJets_;
}

void 
MEMEvent_singlelepton::countGenWJets()
{
  numGenWJets_ = 0;
  if ( genWJet1_ ) ++numGenWJets_;
  if ( genWJet2_ ) ++numGenWJets_;
}

void 
MEMEvent_singlelepton::countMeasuredLeptons()
{
  numMeasuredLeptons_ = 0;
  if ( measuredLepton_ ) ++numMeasuredLeptons_;
}

void 
MEMEvent_singlelepton::countGenLeptons()
{
  numGenLeptons_ = 0;
  if ( genLepton_ ) ++numGenLeptons_;
}

