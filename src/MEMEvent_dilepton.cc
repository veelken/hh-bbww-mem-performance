#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/MEMEvent_dilepton.h"

MEMEvent_dilepton::MEMEvent_dilepton(const EventInfo & eventInfo, bool isSignal,
                                     const mem::MeasuredParticle* measuredBJet1, const mem::MeasuredParticle* measuredBJet2, 
                                     const mem::MeasuredParticle* measuredLepton1, const mem::MeasuredParticle* measuredLepton2, 
                                     double measuredMEtPx, double measuredMEtPy, const TMatrixD& measuredMEtCov)
  : MEMEvent(eventInfo, isSignal, measuredBJet1, measuredBJet2, measuredMEtPx, measuredMEtPy, measuredMEtCov)
  , measuredLepton1_(measuredLepton1)
  , genLepton1_(nullptr)
  , measuredLepton2_(measuredLepton2)
  , genLepton2_(nullptr)
{}

MEMEvent_dilepton::~MEMEvent_dilepton()
{}

void 
MEMEvent_dilepton::set_genLepton1(const GenLepton* genLepton1)
{
  genLepton1_ = genLepton1;
}
  
void 
MEMEvent_dilepton::set_genLepton2(const GenLepton* genLepton2)
{
  genLepton2_ = genLepton2;
}

const mem::MeasuredParticle* 
MEMEvent_dilepton::measuredLepton1() const
{
  return measuredLepton1_;
}

const GenLepton* 
MEMEvent_dilepton::genLepton1() const
{
  return genLepton1_;
}

const mem::MeasuredParticle* 
MEMEvent_dilepton::measuredLepton2() const
{
  return measuredLepton2_;
}

const GenLepton* 
MEMEvent_dilepton::genLepton2() const
{
  return genLepton2_;
}

int 
MEMEvent_dilepton::numMeasuredLeptons() const
{
  return numMeasuredLeptons_;
}
  
int 
MEMEvent_dilepton::numGenLeptons() const
{
  return numGenLeptons_;
}

void 
MEMEvent_dilepton::countMeasuredLeptons()
{
  numMeasuredLeptons_ = 0;
  if ( measuredLepton1_ ) ++numMeasuredLeptons_;
  if ( measuredLepton2_ ) ++numMeasuredLeptons_;
}

void 
MEMEvent_dilepton::countGenLeptons()
{
  numGenLeptons_ = 0;
  if ( genLepton1_ ) ++numGenLeptons_;
  if ( genLepton2_ ) ++numGenLeptons_;
}

