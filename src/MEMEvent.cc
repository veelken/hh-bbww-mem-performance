#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/MEMEvent.h"

MEMEvent::MEMEvent(const EventInfo & eventInfo, bool isSignal,
                   const mem::MeasuredParticle* measuredBJet1, const mem::MeasuredParticle* measuredBJet2, 
                   double measuredMEtPx, double measuredMEtPy, const TMatrixD& measuredMEtCov)
    : eventInfo_(&eventInfo)
    , isSignal_(isSignal)
    , measuredBJet1_(measuredBJet1)
    , genBJet1_(nullptr)
    , measuredBJet2_(measuredBJet2)
    , genBJet2_(nullptr)
    , isBoosted_Hbb_(false)
    , numMeasuredBJets_(0)
    , numGenBJets_(0)
    , measuredMEtPx_(measuredMEtPx)
    , genMEtPx_(0.)
    , measuredMEtPy_(measuredMEtPy)
    , genMEtPy_(0.)
    , measuredMEtCov_(measuredMEtCov)
    , memCpuTime_(-1.)
    , barcode_(-1)
{
  countMeasuredBJets();
}

MEMEvent::~MEMEvent()
{}

void 
MEMEvent::set_genBJet1(const GenJet* genBJet1)
{
  genBJet1_ = genBJet1;
  countGenBJets();
}

void 
MEMEvent::set_genBJet2(const GenJet* genBJet2)
{
  genBJet2_ = genBJet2;
  countGenBJets();
}

void 
MEMEvent::set_isBoosted_Hbb(bool isBoosted_Hbb)
{
  isBoosted_Hbb_ = isBoosted_Hbb;
}

void 
MEMEvent::set_genMEtPx(double genMEtPx)
{
  genMEtPx_ = genMEtPx;
}

void 
MEMEvent::set_genMEtPy(double genMEtPy)
{
  genMEtPy_ = genMEtPy;
}

void 
MEMEvent::set_memResult(const MEMResultBase& memResult)
{
  memResult_ = memResult;
}

void 
MEMEvent::set_memCpuTime(double memCpuTime)
{
  memCpuTime_ = memCpuTime;
}

void 
MEMEvent::set_barcode(int barcode)
{
  barcode_ = barcode; 
}

const EventInfo &
MEMEvent::eventInfo() const
{
  return *eventInfo_;
}

bool 
MEMEvent::isSignal() const
{
  return isSignal_;
}

const mem::MeasuredParticle* 
MEMEvent::measuredBJet1() const
{
  return measuredBJet1_;
}
  
const GenJet* 
MEMEvent::genBJet1() const
{
  return genBJet1_;
}
  
const mem::MeasuredParticle* 
MEMEvent::measuredBJet2() const
{
  return measuredBJet2_;
}

const GenJet* 
MEMEvent::genBJet2() const
{
  return genBJet2_;
}

bool 
MEMEvent::isBoosted_Hbb() const
{
  return isBoosted_Hbb_;
}
 
int 
MEMEvent::numMeasuredBJets() const
{
  return numMeasuredBJets_;
}
  
int 
MEMEvent::numGenBJets() const
{
  return numGenBJets_;
}
  
double 
MEMEvent::measuredMEtPx() const
{
  return measuredMEtPx_;
}

double 
MEMEvent::genMEtPx() const
{
  return genMEtPx_;
}
  
double 
MEMEvent::measuredMEtPy() const
{
  return measuredMEtPy_;
}

double 
MEMEvent::genMEtPy() const
{
  return genMEtPy_;
}

const TMatrixD &
MEMEvent::measuredMEtCov() const
{
  return measuredMEtCov_;
}
  
const MEMResultBase &
MEMEvent::memResult() const
{
  return memResult_;
}
  
double 
MEMEvent::memCpuTime() const
{
  return memCpuTime_;
}
  
int 
MEMEvent::barcode() const
{
  return barcode_;
}

void 
MEMEvent::countMeasuredBJets()
{
  numMeasuredBJets_ = 0;
  if ( measuredBJet1_ ) ++numMeasuredBJets_;
  if ( measuredBJet2_ ) ++numMeasuredBJets_;
}

void 
MEMEvent::countGenBJets()
{
  numGenBJets_ = 0;
  if ( genBJet1_ ) ++numGenBJets_;
  if ( genBJet2_ ) ++numGenBJets_;
}
