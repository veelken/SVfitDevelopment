#ifndef TauAnalysis_SVfitDevelopment_SVfitEventHypothesisWatchdog_h
#define TauAnalysis_SVfitDevelopment_SVfitEventHypothesisWatchdog_h

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include <TMath.h>

#include <string>

class SVfitEventHypothesisWatchdog : public edm::EDAnalyzer 
{
 public:
  // constructor 
  explicit SVfitEventHypothesisWatchdog(const edm::ParameterSet&);
    
  // destructor
  ~SVfitEventHypothesisWatchdog() {}
    
 private:
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

  std::string moduleLabel_;

  edm::InputTag srcEventHypotheses1_;
  edm::InputTag srcEventHypotheses2_;
 
  edm::InputTag srcGenTauPairs_;

  int idxResonance_;

  double threshold1MassUp_;
  double threshold1MassDown_;
  double threshold1Sigma_;
  double threshold2MassUp_;
  double threshold2MassDown_;
  double threshold2Sigma_;

  long numEvents_processed_;
  long numEvents_watchdogTriggered1_;
  long numEvents_watchdogTriggered2_;
};

#endif   
