#include "TauAnalysis/SVfitDevelopment/plugins/SVfitEventHypothesisWatchdog.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/View.h"

#include "AnalysisDataFormats/SVfit/interface/SVfitEventHypothesisBase.h"
#include "AnalysisDataFormats/SVfit/interface/SVfitResonanceHypothesisBase.h"
#include "AnalysisDataFormats/SVfit/interface/SVfitSingleParticleHypothesis.h"

SVfitEventHypothesisWatchdog::SVfitEventHypothesisWatchdog(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label")),
    numEvents_processed_(0),
    numEvents_watchdogTriggered1_(0),
    numEvents_watchdogTriggered2_(0)
{
  srcEventHypotheses1_ = cfg.getParameter<edm::InputTag>("srcEventHypotheses1");
  srcEventHypotheses2_ = cfg.getParameter<edm::InputTag>("srcEventHypotheses2");

  srcGenTauPairs_ = cfg.getParameter<edm::InputTag>("srcGenTauPairs");
  
  idxResonance_ = ( cfg.exists("idxResonance") ) ?
    cfg.getParameter<int>("idxResonance") : 0;

  threshold1MassUp_   = cfg.getParameter<double>("threshold1MassUp");
  threshold1MassDown_ = cfg.getParameter<double>("threshold1MassDown");
  threshold1Sigma_    = cfg.getParameter<double>("threshold1Sigma");
  threshold2MassUp_   = cfg.getParameter<double>("threshold2MassUp");
  threshold2MassDown_ = cfg.getParameter<double>("threshold2MassDown");
  threshold2Sigma_    = cfg.getParameter<double>("threshold2Sigma");
}

namespace
{
  double square(double x)
  {
    return x*x;
  }
}

void SVfitEventHypothesisWatchdog::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  typedef edm::View<SVfitEventHypothesisBase> SVfitEventHypothesisCollection;
  edm::Handle<SVfitEventHypothesisCollection> svFitEventHypotheses1;
  evt.getByLabel(srcEventHypotheses1_, svFitEventHypotheses1);
  
  edm::Handle<SVfitEventHypothesisCollection> svFitEventHypotheses2;
  evt.getByLabel(srcEventHypotheses2_, svFitEventHypotheses2);

  typedef edm::View<reco::Candidate> CandidateView;
  edm::Handle<CandidateView> genTauPairs;
  evt.getByLabel(srcGenTauPairs_, genTauPairs);

  if ( !(svFitEventHypotheses1->size() == svFitEventHypotheses2->size()) )
    throw cms::Exception("SVfitEventHypothesisAnalyzer") 
      << " Different number of entries in SVfitEventHypothesis collections " 
      << srcEventHypotheses1_.label() << " (" << svFitEventHypotheses1->size() << " entries) and"
      << srcEventHypotheses2_.label() << " (" << svFitEventHypotheses2->size() << " entries) !!\n";

  bool isTriggered1_event = false;
  bool isTriggered2_event = false;
  
  size_t numSVfitEventHypotheses = svFitEventHypotheses1->size();
  for ( size_t iSVfitEventHypothesis = 0; iSVfitEventHypothesis < numSVfitEventHypotheses; ++iSVfitEventHypothesis ) {
    const SVfitEventHypothesisBase* svFitEventHypothesis1 = &svFitEventHypotheses1->at(iSVfitEventHypothesis);
    if ( !((int)svFitEventHypothesis1->numResonances() > idxResonance_) )
      throw cms::Exception("SVfitEventHypothesisWatchdog") 
	<< " Failed to find resonance #" << idxResonance_ << " in first SVfitEventHypothesis !!\n";
    const SVfitResonanceHypothesisBase* svFitResonanceHypothesis1 = svFitEventHypothesis1->resonance(idxResonance_);
    bool svFitIsValidSolution1 = svFitResonanceHypothesis1->isValidSolution();
    double svFitMass1      = svFitResonanceHypothesis1->mass();
    double svFitSigmaUp1   = svFitResonanceHypothesis1->massErrUp(); 
    double svFitSigmaDown1 = svFitResonanceHypothesis1->massErrDown();
    double svFitSigma1     = TMath::Sqrt(square(svFitSigmaUp1) + square(svFitSigmaDown1));

    const SVfitEventHypothesisBase* svFitEventHypothesis2 = &svFitEventHypotheses2->at(iSVfitEventHypothesis);
    if ( !((int)svFitEventHypothesis2->numResonances() > idxResonance_) )
      throw cms::Exception("SVfitEventHypothesisWatchdog") 
	<< " Failed to find resonance #" << idxResonance_ << " in second SVfitEventHypothesis !!\n";
    const SVfitResonanceHypothesisBase* svFitResonanceHypothesis2 = svFitEventHypothesis2->resonance(idxResonance_);
    bool svFitIsValidSolution2 = svFitResonanceHypothesis2->isValidSolution();
    double svFitMass2      = svFitResonanceHypothesis2->mass();
    double svFitSigmaUp2   = svFitResonanceHypothesis2->massErrUp(); 
    double svFitSigmaDown2 = svFitResonanceHypothesis2->massErrDown();
    double svFitSigma2     = TMath::Sqrt(square(svFitSigmaUp2) + square(svFitSigmaDown2));

    assert(svFitResonanceHypothesis1->numDaughters() == 2);
    const SVfitSingleParticleHypothesis* svFitDaughter1 = dynamic_cast<const SVfitSingleParticleHypothesis*>(
      svFitResonanceHypothesis1->daughter(0));
    const reco::Candidate::LorentzVector& svFitDaughter1P4 = svFitDaughter1->p4();
    const SVfitSingleParticleHypothesis* svFitDaughter2 = dynamic_cast<const SVfitSingleParticleHypothesis*>(
      svFitResonanceHypothesis1->daughter(1));
    const reco::Candidate::LorentzVector& svFitDaughter2P4 = svFitDaughter2->p4();

    for ( CandidateView::const_iterator genTauPair = genTauPairs->begin();
	  genTauPair != genTauPairs->end(); ++genTauPair ) {

      double genMass = genTauPair->mass();
      
      bool isTriggered1 = false;
      bool isTriggered2 = false;
      std::string message;
      if        ( TMath::Abs(svFitMass1 - genMass) > (TMath::Abs(svFitMass2 - genMass) + threshold1MassUp_  ) && svFitMass1 > genMass ) {
	message = "mass1 too high";
	isTriggered1 = true;
      } else if ( TMath::Abs(svFitMass1 - genMass) > (TMath::Abs(svFitMass2 - genMass) + threshold1MassDown_) && svFitMass1 < genMass ) {
	message = "mass1 too low";
	isTriggered1 = true;
      } else if ( svFitSigma1 > (svFitSigma2 + threshold1Sigma_) ) {
	message = "sigma1 too high";
	isTriggered1 = true;
      } else if ( TMath::Abs(svFitMass2 - genMass) > (TMath::Abs(svFitMass1 - genMass) + threshold2MassUp_  ) && svFitMass2 > genMass ) {
	message = "mass2 too high";
	isTriggered2 = true;
      } else if ( TMath::Abs(svFitMass2 - genMass) > (TMath::Abs(svFitMass1 - genMass) + threshold2MassDown_) && svFitMass2 < genMass ) {
	message = "mass2 too low";
	isTriggered2 = true;
      } else if ( svFitSigma2 > (svFitSigma2 + threshold1Sigma_) ) {
	message = "sigma2 too high";
	isTriggered2 = true;
      } 
      if ( isTriggered1 || isTriggered2 ) {
	std::cout << "<SVfitEventHypothesisWatchdog::analyze (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
	std::cout << " " << message << "." << std::endl;
	std::cout << " src1 = " << srcEventHypotheses1_.label() << ": mass = " << svFitMass1 << " + " << svFitSigmaUp1 << " - " << svFitSigmaDown1 << std::endl;
	//svFitEventHypothesis1->print(std::cout);
	std::cout << " src2 = " << srcEventHypotheses2_.label() << ": mass = " << svFitMass2 << " + " << svFitSigmaUp2 << " - " << svFitSigmaDown2 << std::endl;
	//svFitEventHypothesis2->print(std::cout);
	if ( isTriggered1 ) isTriggered1_event = true;
	if ( isTriggered2 ) isTriggered2_event = true;
      }
    }
  }

  ++numEvents_processed_;
  if ( isTriggered1_event ) ++numEvents_watchdogTriggered1_;
  if ( isTriggered2_event ) ++numEvents_watchdogTriggered2_;
}

void SVfitEventHypothesisWatchdog::endJob()
{
  std::cout << "<SVfitEventHypothesisWatchdog::endJob (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
  std::cout << " num. Events processed = " << numEvents_processed_ << std::endl;
  std::cout << " num. Events watchdog triggered (" << srcEventHypotheses1_.label() << ") = " << numEvents_watchdogTriggered1_ << std::endl;
  std::cout << " num. Events watchdog triggered (" << srcEventHypotheses2_.label() << ") = " << numEvents_watchdogTriggered2_ << std::endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(SVfitEventHypothesisWatchdog);





