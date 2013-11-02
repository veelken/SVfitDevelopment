#include "TauAnalysis/SVfitDevelopment/plugins/SVfitBestEventHypothesisSelectorT.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitResonanceHypothesisBase.h"
#include "AnalysisDataFormats/SVfit/interface/SVfitSingleParticleHypothesis.h"
#include "DataFormats/METReco/interface/PFMEtSignCovMatrix.h"

#include <TMatrixD.h>
#include <TVectorD.h>

template <typename T>
SVfitBestEventHypothesisSelectorT<T>::SVfitBestEventHypothesisSelectorT(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label")),
    selector_(0),
    rank_(0)
{
  src_ = cfg.getParameter<vInputTag>("src");
  
  if ( cfg.exists("selection") ) {
    std::string selection_string = cfg.getParameter<std::string>("selection");
    selector_ = new StringCutObjectSelector<T>(selection_string);
  }
  std::string rank_string = cfg.getParameter<std::string>("rank");
  rank_ = new StringObjectFunction<T>(rank_string);

  verbosity_ = cfg.exists("verbosity") ?
    cfg.getParameter<int>("verbosity") : 0;

  produces<SVfitEventHypothesisCollection>();
}

template <typename T>
SVfitBestEventHypothesisSelectorT<T>::~SVfitBestEventHypothesisSelectorT()
{
  delete selector_;
  delete rank_;
}

namespace
{
  template <typename T>
  void print(const std::string& label, const T& svFitEventHypothesis, double rank)
  {
    const SVfitResonanceHypothesisBase* svFitResonanceHypothesis = svFitEventHypothesis.resonance(0);
    bool svFitIsValidSolution = svFitResonanceHypothesis->isValidSolution();
    double svFitMass      = svFitResonanceHypothesis->mass();
    double svFitSigmaUp   = svFitResonanceHypothesis->massErrUp(); 
    double svFitSigmaDown = svFitResonanceHypothesis->massErrDown();
    std::cout << " " << label << " #" << svFitResonanceHypothesis->barcode() << " (rank = " << rank << "):" 
	      << " mass = " << svFitMass << " + " << svFitSigmaUp << " - " << svFitSigmaDown << " (isValidSolution = " << svFitIsValidSolution <<")" << std::endl;
  }
}

template <typename T>
void SVfitBestEventHypothesisSelectorT<T>::produce(edm::Event& evt, const edm::EventSetup& es)
{
  if ( verbosity_ ) {
    std::cout << "<SVfitBestEventHypothesisSelectorT::produce (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
  }

  std::vector<const T*> svFitEventHypotheses_output;
  std::vector<double> svFitEventHypotheses_rank;
  int numEventHypotheses = -1;

  if ( verbosity_ ) {
    std::cout << "input:" << std::endl;
  }
  for ( vInputTag::const_iterator srcEventHypotheses = src_.begin();
	srcEventHypotheses != src_.end(); ++srcEventHypotheses ) {
    edm::Handle<SVfitEventHypothesisCollection> svFitEventHypotheses;
    evt.getByLabel(*srcEventHypotheses, svFitEventHypotheses);
    
    if ( numEventHypotheses != -1 ) {
      if ( (int)svFitEventHypotheses->size() != numEventHypotheses )
	throw cms::Exception("SVfitBestEventHypothesisSelectorT::analyze") 
	  << "SVfitEventHypothesis collections must be of equal size !!\n";
    } else {      
      numEventHypotheses = svFitEventHypotheses->size();
      if ( numEventHypotheses >= 1 ) {
	svFitEventHypotheses_output.resize(numEventHypotheses, 0);
	svFitEventHypotheses_rank.resize(numEventHypotheses);
      }
    }

    int idx = 0;
    for ( typename SVfitEventHypothesisCollection::const_iterator svFitEventHypothesis = svFitEventHypotheses->begin();
	  svFitEventHypothesis != svFitEventHypotheses->end(); ++svFitEventHypothesis ) {

      if ( selector_ && !(*selector_)(*svFitEventHypothesis) ) continue;

      double rank = (*rank_)(*svFitEventHypothesis);
      if ( verbosity_ ) {
	print(srcEventHypotheses->label(), *svFitEventHypothesis, rank);
      }
      
      if ( svFitEventHypotheses_output[idx] == 0 || rank > svFitEventHypotheses_rank[idx] ) {
	svFitEventHypotheses_output[idx] = &(*svFitEventHypothesis);
	svFitEventHypotheses_rank[idx] = rank;
      }

      ++idx;
    }
  }

  std::auto_ptr<SVfitEventHypothesisCollection> output(new SVfitEventHypothesisCollection);

  if ( verbosity_ ) {
    std::cout << "output:" << std::endl;
  }
  int idx = 0;
  for ( typename std::vector<const T*>::const_iterator svFitEventHypothesis = svFitEventHypotheses_output.begin();
	svFitEventHypothesis != svFitEventHypotheses_output.end(); ++svFitEventHypothesis ) {
    if ( (*svFitEventHypothesis) != 0 ) {
      if ( verbosity_ ) {
	print(moduleLabel_, **svFitEventHypothesis, svFitEventHypotheses_rank[idx]);
      }
      output->push_back(**svFitEventHypothesis);
    }
    ++idx;
  }

  evt.put(output);
}

#include "AnalysisDataFormats/SVfit/interface/SVfitEventHypothesis.h"
#include "AnalysisDataFormats/SVfit/interface/SVfitEventHypothesisByIntegration.h"

typedef SVfitBestEventHypothesisSelectorT<SVfitEventHypothesis> SVfitBestEventHypothesisSelector;
typedef SVfitBestEventHypothesisSelectorT<SVfitEventHypothesisByIntegration> SVfitBestEventHypothesisByIntegrationSelector;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(SVfitBestEventHypothesisSelector);
DEFINE_FWK_MODULE(SVfitBestEventHypothesisByIntegrationSelector);





