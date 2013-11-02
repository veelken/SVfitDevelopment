#ifndef TauAnalysis_SVfitDevelopment_SVfitBestEventHypothesisSelectorT_h
#define TauAnalysis_SVfitDevelopment_SVfitBestEventHypothesisSelectorT_h

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/normalizedPhi.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitResonanceHypothesisBase.h"
#include "AnalysisDataFormats/SVfit/interface/SVfitResonanceHypothesis.h"
#include "AnalysisDataFormats/SVfit/interface/SVfitResonanceHypothesisByIntegration.h"
#include "AnalysisDataFormats/SVfit/interface/SVfitTauDecayHypothesis.h"
#include "TauAnalysis/SVfit/interface/svFitAuxFunctions.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include <TMath.h>

#include <string>

template <typename T>
class SVfitBestEventHypothesisSelectorT : public edm::EDProducer 
{
  typedef std::vector<T> SVfitEventHypothesisCollection;

 public:
  // constructor 
  explicit SVfitBestEventHypothesisSelectorT(const edm::ParameterSet&);
    
  // destructor
  ~SVfitBestEventHypothesisSelectorT();
    
 private:
  void produce(edm::Event&, const edm::EventSetup&);

  std::string moduleLabel_;

  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag src_;

  const StringCutObjectSelector<T>* selector_;
  const StringObjectFunction<T>* rank_;

  int verbosity_;
};

#endif   
