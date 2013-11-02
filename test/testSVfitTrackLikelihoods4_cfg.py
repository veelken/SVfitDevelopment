import FWCore.ParameterSet.Config as cms

process = cms.Process("testSVfitTrackLikelihoods4")

import os
import re

import TauAnalysis.Configuration.tools.castor as castor
import TauAnalysis.Configuration.tools.eos as eos
from TauAnalysis.Skimming.EventContent_cff import *

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.load('Configuration/Geometry/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START53_V15::All')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/data1/veelken/CMSSW_5_3_x/skims/selEvents_testSVfitTrackLikelihoods_Z_tautau_90_2013Aug08_AOD.root'
    ),
    dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
    inputCommands=cms.untracked.vstring(
        'keep *',
        'drop LHERunInfoProduct_*_*_*',                        
        'drop LHEEventProduct_*_*_*'
    ),
    ##eventsToProcess = cms.untracked.VEventRange(
    ##    '1:102093:40804761',
    ##)
)

#--------------------------------------------------------------------------------
# define configuration parameter default values

sample_type = 'Z'
#sample_type = 'Higgs'
#channel = 'etau'
#channel = 'mutau'
channel = 'tautau'
#channel = 'emu'
massPoint = 125
#massPoint = 300
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define "hooks" for replacing configuration parameters
# in case running jobs on the CERN batch system/grid
#
#__sample_type = '#sample_type#'
#__channel = '#channel#'
#__massPoint = #massPoint#
#
if sample_type == 'Z':
    massPoint = 90
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# set input files
inputFilePath = None
inputFile_regex = None
if sample_type == 'Z':
    if channel == 'mutau':
        inputFilePaths = [
            '/store/user/veelken/CMSSW_5_3_x/skims/SVfitStudies/AccLowPtThresholds/ZplusJets_mutau/'
        ]
    elif channel == 'tautau':
        inputFilePaths = [
            '/store/user/veelken/CMSSW_5_3_x/skims/SVfitStudies/AccLowPtThresholds/ZplusJets_tautau_highStats/'
        ]
    else:
        raise ValueError("Invalid channel = %s !!" % channel)
    inputFile_regex = \
      r"[a-zA-Z0-9_/:.]*genTauLeptonPairSkim_ZplusJets_%s_(?P<gridJob>\d*)(_(?P<gridTry>\d*))*_(?P<hash>[a-zA-Z0-9]*).root" % channel
elif sample_type == 'Higgs':
    if channel == 'mutau':
        inputFilePaths = [
            '/store/user/veelken/CMSSW_5_3_x/skims/SVfitStudies/AccLowPtThresholds/ggHiggs125_mutau/',
            '/store/user/veelken/CMSSW_5_3_x/skims/SVfitStudies/AccLowPtThresholds/vbfHiggs125_mutau/'
        ]
    elif channel == 'tautau':    
        inputFilePaths = [
            '/store/user/veelken/CMSSW_5_3_x/skims/SVfitStudies/AccLowPtThresholds/ggHiggs125_tautau_highStats/',
            '/store/user/veelken/CMSSW_5_3_x/skims/SVfitStudies/AccLowPtThresholds/vbfHiggs125_tautau_highStats/'
        ]
    else:
        raise ValueError("Invalid channel = %s !!" % channel)
    inputFile_regex = \
      r"[a-zA-Z0-9_/:.]*genTauLeptonPairSkim_(ggHiggs|ggPhi|vbfHiggs)%i_%s_(?P<gridJob>\d*)(_(?P<gridTry>\d*))*_(?P<hash>[a-zA-Z0-9]*).root" % (massPoint, channel)
else:
    raise ValueError("Invalid sample type = %s !!" % sample_type)

# check if name of inputFile matches regular expression
inputFileNames = []
for inputFilePath in inputFilePaths:
    files = None
    if inputFilePath.startswith('/castor/'):
        files = [ "".join([ "rfio:", file_info['path'] ]) for file_info in castor.nslsl(inputFilePath) ]
    elif inputFilePath.startswith('/store/'):
        files = [ file_info['path'] for file_info in eos.lsl(inputFilePath) ]
    else:
        files = [ "".join([ "file:", inputFilePath, file ]) for file in os.listdir(inputFilePath) ]
    for file in files:
        #print "file = %s" % file
        inputFile_matcher = re.compile(inputFile_regex)
        if inputFile_matcher.match(file):
           inputFileNames.append(file)
print "inputFileNames = %s" % inputFileNames 

process.source.fileNames = cms.untracked.vstring(inputFileNames)
#--------------------------------------------------------------------------------

process.testSVfitTrackLikelihoodProductionSequence = cms.Sequence()

#--------------------------------------------------------------------------------
# require event to contain generator level tau lepton pair,
# decaying to the specified decay mode
# matched to tau decay products on generator level

process.load("PhysicsTools/JetMCAlgos/TauGenJets_cfi")
process.testSVfitTrackLikelihoodProductionSequence += process.tauGenJets
process.load("PhysicsTools/JetMCAlgos/TauGenJetsDecayModeSelectorAllHadrons_cfi")
process.tauGenJetsSelectorAllHadrons.select = cms.vstring(
    'oneProng0Pi0', 
    'oneProng1Pi0', 
    'oneProng2Pi0', 
    'oneProngOther',
    'threeProng0Pi0', 
    'threeProng1Pi0', 
    'threeProngOther', 
    'rare'
)
process.tauGenJetsSelectorAllHadrons.filter = cms.bool(True)
process.testSVfitTrackLikelihoodProductionSequence += process.tauGenJetsSelectorAllHadrons

genElectronsFromTauDecays = None
genMuonsFromTauDecays = None
genTauJetsFromTauDecays = None
genTauPairs = None
genTaus = None
if sample_type == 'Z':
    process.load("TauAnalysis/GenSimTools/gen_decaysFromZs_cfi")
    ##process.genTauJetsFromZs.srcTauGenJets = cms.InputTag('tauGenJetsSelectorAllHadrons')
    process.testSVfitTrackLikelihoodProductionSequence += process.produceGenDecayProductsFromZs
    genElectronsFromTauDecays = 'genElectronsFromZtautauDecays'
    genMuonsFromTauDecays = 'genMuonsFromZtautauDecays'
    genTauJetsFromTauDecays = 'genHadronsFromZtautauDecays'
    genTauPairs = 'genZdecayToTaus'
    genTaus = 'genTausFromZs'
elif sample_type == 'Higgs':
    process.load("TauAnalysis/GenSimTools/gen_decaysFromAHs_cfi")
    ##process.genTauJetsFromAHs.srcTauGenJets = cms.InputTag('tauGenJetsSelectorAllHadrons')
    process.testSVfitTrackLikelihoodProductionSequence += process.produceGenDecayProductsFromAHs
    genElectronsFromTauDecays = 'genElectronsFromAHtautauDecays'
    genMuonsFromTauDecays = 'genMuonsFromAHtautauDecays'
    genTauJetsFromTauDecays = 'genHadronsFromAHtautauDecays'
    genTauPairs = 'genAHdecayToTaus'
    genTaus = 'genTausFromAHs'
else:
    raise ValueError("Invalid sample type = %s !!" % sample_type)

# CV: require that generator level tau lepton pair is within 0.70..1.30 of "nominal" mass point,
#     in order to cut low mass tail present in MSSM gg -> Phi Monte Carlo samples
process.genTauPairMassFilter = cms.EDFilter("CandViewSelector",
    src = cms.InputTag(genTauPairs),
    cut = cms.string('mass > %1.1f & mass < %1.1f' % (0.70*massPoint, 1.30*massPoint)),
    filter = cms.bool(True)
)
process.testSVfitTrackLikelihoodProductionSequence += process.genTauPairMassFilter

numElectrons = None
numMuons     = None
numTauJets   = None
if channel == 'mutau':
    numElectrons = 0
    numMuons     = 1
    numTauJets   = 1
elif channel == 'etau':
    numElectrons = 1
    numMuons     = 0
    numTauJets   = 1
elif channel == 'emu':
    numElectrons = 1
    numMuons     = 1
    numTauJets   = 0
elif channel == 'tautau':
    numElectrons = 0
    numMuons     = 0
    numTauJets   = 2  
else:
    raise ValueError("Invalid channel = %s !!" % channel)

if numElectrons > 0:
    process.genElectronFilter = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag(genElectronsFromTauDecays),
        minNumber = cms.uint32(numElectrons),
        maxNumber = cms.uint32(numElectrons)                      
    ) 
    process.testSVfitTrackLikelihoodProductionSequence += process.genElectronFilter

if numMuons > 0:    
    process.genMuonFilter = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag(genMuonsFromTauDecays),
        minNumber = cms.uint32(numMuons),
        maxNumber = cms.uint32(numMuons)                      
    )
    process.testSVfitTrackLikelihoodProductionSequence += process.genMuonFilter

if numTauJets > 0:
    process.genTauFilter = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag('tauGenJetsSelectorAllHadrons'),
        minNumber = cms.uint32(numTauJets),
        maxNumber = cms.uint32(numTauJets)                      
    )
    process.testSVfitTrackLikelihoodProductionSequence += process.genTauFilter
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# select collections of electrons, muons and tau-jets
# matching genuine tau -> e, tau -> mu and tau -> hadronic decays on generator level

process.genTauMatchedPFJets = cms.EDFilter("PFJetAntiOverlapSelector",
    src = cms.InputTag('ak5PFJets'),
    srcNotToBeFiltered = cms.VInputTag('tauGenJetsSelectorAllHadrons'),
    dRmin = cms.double(0.3),
    invert = cms.bool(True),
    filter = cms.bool(False)                                                          
)
process.testSVfitTrackLikelihoodProductionSequence += process.genTauMatchedPFJets

process.load("RecoTauTag/Configuration/RecoPFTauTag_cff")

process.ak5PFJetTracksAssociatorAtVertex.jets = cms.InputTag('genTauMatchedPFJets')
process.ak5PFJetsLegacyHPSPiZeros.jetSrc = cms.InputTag('genTauMatchedPFJets')
process.recoTauAK5PFJets08Region.src = cms.InputTag('genTauMatchedPFJets')
process.ak5PFJetsRecoTauChargedHadrons.jetSrc = cms.InputTag('genTauMatchedPFJets')
process.combinatoricRecoTaus.jetSrc = cms.InputTag('genTauMatchedPFJets')

##process.combinatoricRecoTaus.modifiers[2].verbosity = cms.int32(1)
##process.hpsSelectionDiscriminator.verbosity = cms.int32(1)
##process.hpsPFTauProducerSansRefs.cleaners[2].verbosity = cms.int32(1)
##process.hpsPFTauProducerSansRefs.verbosity = cms.int32(1)

process.testSVfitTrackLikelihoodProductionSequence += process.recoTauClassicHPSSequence

process.load("TauAnalysis.Skimming.goldenZmmSelectionVBTFnoMuonIsolation_cfi")
process.goodMuons.cut = cms.string(
    'pt > 1. & abs(eta) < 2.5 & isGlobalMuon && isPFMuon' \
    + ' && track.hitPattern.trackerLayersWithMeasurement > 5 & innerTrack.hitPattern.numberOfValidPixelHits > 0' \
    + ' && abs(dB) < 0.2 && globalTrack.normalizedChi2 < 10' \
    + ' && globalTrack.hitPattern.numberOfValidMuonHits > 0 && numberOfMatchedStations > 1'
)
process.goodMuons.filter = cms.bool(False)
process.goodIsoMuons.cut = cms.string(
    '(userIsolation("pat::User1Iso")' + \
    ' + max(0., userIsolation("pat::PfNeutralHadronIso") + userIsolation("pat::PfGammaIso")' + \
    '          - 0.5*userIsolation("pat::User2Iso"))) < 0.10*pt'
)
process.goodIsoMuons.filter = cms.bool(False)
process.muonSelectionSequence = cms.Sequence(
    process.pfNoPileUpSequence
   + process.pfParticleSelectionSequence
   + process.muonPFIsolationSequenceForPAT
   + process.patMuonsForGoldenZmmSelection
   + process.goodMuons
   + process.goodIsoMuons
)
process.testSVfitTrackLikelihoodProductionSequence += process.muonSelectionSequence

process.genMatchedPatMuons = cms.EDFilter("PATMuonAntiOverlapSelector",
    src = cms.InputTag('goodIsoMuons'),
    srcNotToBeFiltered = cms.VInputTag(genMuonsFromTauDecays),
    dRmin = cms.double(0.3),
    invert = cms.bool(True),
    filter = cms.bool(False)                                                          
)                        
process.testSVfitTrackLikelihoodProductionSequence += process.genMatchedPatMuons

process.load("PhysicsTools/PatAlgos/producersLayer1/electronProducer_cfi")
process.patElectronsForSVfitTrackLikelihood = process.patElectrons.clone(
    isoDeposits = cms.PSet(),
    addGenMatch = cms.bool(False),
    embedHighLevelSelection = cms.bool(True),
    usePV = cms.bool(False) # compute transverse impact parameter wrt. beamspot (not event vertex)
)
from TauAnalysis.RecoTools.recoVertexSelection_cff import *
process.selectedPrimaryVertexQualityForElectronId = selectedPrimaryVertexQuality.clone(
    src = cms.InputTag('offlinePrimaryVerticesWithBS'),
    cut = cms.string("isValid & ndof >= 7 & chi2 > 0 & tracksSize > 0")
)
process.selectedPrimaryVertexPositionForElectronId = selectedPrimaryVertexPosition.clone(
    src = cms.InputTag('selectedPrimaryVertexQualityForElectronId')
)
process.selectedPrimaryVertexHighestPtTrackSumForElectronId = selectedPrimaryVertexHighestPtTrackSum.clone(
    vertices = cms.InputTag('selectedPrimaryVertexPositionForElectronId')
)
process.load("TauAnalysis/RecoTools/patLeptonSelection_cff")
process.selectedPatElectronsForElecTauId.srcVertex = cms.InputTag('selectedPrimaryVertexHighestPtTrackSumForElectronId')
from TauAnalysis.CandidateTools.tools.objSelConfigurator import *
patElectronSelConfiguratorForSVfitTrackLikelihood = objSelConfigurator(
    [ process.selectedPatElectronsForElecTauId,
      process.selectedPatElectronsForElecTauAntiCrackCut,
      process.selectedPatElectronsForElecTauEta,
      process.selectedPatElectronsForElecTauIso,
      process.selectedPatElectronsForElecTauConversionVeto ],
    src = 'patElectronsForSVfitTrackLikelihood',
    doSelIndividual = False
)
process.selectPatElectronsForSVfitTrackLikelihood = patElectronSelConfiguratorForSVfitTrackLikelihood.configure(process = process)
process.electronSelectionSequence = cms.Sequence(
    process.patElectronsForSVfitTrackLikelihood
   + process.selectedPrimaryVertexQualityForElectronId + process.selectedPrimaryVertexPositionForElectronId + process.selectedPrimaryVertexHighestPtTrackSumForElectronId
   + process.selectPatElectronsForSVfitTrackLikelihood
)
process.testSVfitTrackLikelihoodProductionSequence += process.electronSelectionSequence

process.genMatchedPatElectrons = cms.EDFilter("PATElectronAntiOverlapSelector",
    src = cms.InputTag('selectedPatElectronsForElecTauConversionVetoCumulative'),                  
    srcNotToBeFiltered = cms.VInputTag(genElectronsFromTauDecays),
    dRmin = cms.double(0.3),
    invert = cms.bool(True),
    filter = cms.bool(False)                                                          
)
process.testSVfitTrackLikelihoodProductionSequence += process.genMatchedPatElectrons

from PhysicsTools.PatAlgos.producersLayer1.tauProducer_cfi import patTaus
process.patTausForSVfit = patTaus.clone(
    tauSource = cms.InputTag('hpsPFTauProducer'),
    isoDeposits = cms.PSet(),
    userIsolation = cms.PSet(),
    addTauID = cms.bool(True),
    tauIDSources = cms.PSet(
        decayModeFinding = cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding'),
        byLooseCombinedIsolationDeltaBetaCorr = cms.InputTag('hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr')
    ),
    addGenMatch = cms.bool(False),
    embedGenMatch = cms.bool(False),
    addGenJetMatch = cms.bool(False),
    embedGenJetMatch = cms.bool(False)
)
process.testSVfitTrackLikelihoodProductionSequence += process.patTausForSVfit

process.selectedPatTausForSVfit = cms.EDProducer("PATTauCleaner",
    src = cms.InputTag('patTausForSVfit'),
    ##preselection = cms.string("tauID('decayModeFinding') > 0.5 & tauID('byLooseCombinedIsolationDeltaBetaCorr') > 0.5"),                                             
    preselection = cms.string("tauID('decayModeFinding') > 0.5 & tauID('byLooseCombinedIsolationDeltaBetaCorr') > 0.5 & decayMode >= 10"),
    checkOverlaps = cms.PSet(),
    finalCut = cms.string('pt > 20. & abs(eta) < 2.3'),
)
process.testSVfitTrackLikelihoodProductionSequence += process.selectedPatTausForSVfit

process.genMatchedPatTaus = cms.EDFilter("PATTauAntiOverlapSelector",
    src = cms.InputTag('selectedPatTausForSVfit'),
    srcNotToBeFiltered = cms.VInputTag(genTauJetsFromTauDecays),
    dRmin = cms.double(0.3),
    invert = cms.bool(True),
    filter = cms.bool(False)                                                          
)
process.testSVfitTrackLikelihoodProductionSequence += process.genMatchedPatTaus
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# require event to contain reconstructed lepton pair,
# matched to tau decay products on generator level

if numElectrons > 0:
    process.recElectronFilter = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag('genMatchedPatElectrons'),
        minNumber = cms.uint32(numElectrons),
        maxNumber = cms.uint32(numElectrons)                      
    ) 
    process.testSVfitTrackLikelihoodProductionSequence += process.recElectronFilter

if numMuons > 0:    
    process.recMuonFilter = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag('genMatchedPatMuons'),
        minNumber = cms.uint32(numMuons),
        maxNumber = cms.uint32(numMuons)                      
    )
    process.testSVfitTrackLikelihoodProductionSequence += process.recMuonFilter

if numTauJets > 0:
    process.recTauFilter = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag('genMatchedPatTaus'),
        minNumber = cms.uint32(numTauJets),
        maxNumber = cms.uint32(numTauJets)                      
    )
    process.testSVfitTrackLikelihoodProductionSequence += process.recTauFilter
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# produce genMET
process.load("RecoMET.Configuration.GenMETParticles_cff")
process.testSVfitTrackLikelihoodProductionSequence += process.genParticlesForMETAllVisible

process.load("RecoMET.METProducers.genMetTrue_cfi")
process.genMetFromGenParticles = process.genMetTrue.clone(
    src = cms.InputTag('genParticlesForMETAllVisible'),
    alias = cms.string('genMetFromGenParticles')
)
process.testSVfitTrackLikelihoodProductionSequence += process.genMetFromGenParticles
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# reconstructed Type-1 corrected PFMET and its uncertainty 

process.load("JetMETCorrections/Type1MET/pfMETCorrectionType0_cfi")
process.testSVfitTrackLikelihoodProductionSequence += process.type0PFMEtCorrection

process.load("JetMETCorrections/Type1MET/pfMETsysShiftCorrections_cfi")
process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runAplusBvsNvtx_mc
process.testSVfitTrackLikelihoodProductionSequence += process.pfMEtSysShiftCorrSequence

process.load("JetMETCorrections/Type1MET/pfMETCorrections_cff")
process.pfType1CorrectedMet.applyType0Corrections = cms.bool(False)
process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag(
    cms.InputTag('pfJetMETcorr', 'type1'),
    cms.InputTag('pfMETcorrType0'),
    cms.InputTag('pfMEtSysShiftCorr')
)
process.testSVfitTrackLikelihoodProductionSequence += process.producePFMETCorrections

# CV: compute PFMET significance cov. matrix for uncorrected jets
#     in order to include pile-up jets
#    (which to a significant fraction get killed by L1Fastjet corrections)
process.ak5PFJetsNotOverlappingWithLeptons = cms.EDFilter("PFJetAntiOverlapSelector",
    src = cms.InputTag('ak5PFJets'),
    srcNotToBeFiltered = cms.VInputTag(
        'genMatchedPatElectrons',
        'genMatchedPatMuons',
        'genMatchedPatTaus'
    ),
    dRmin = cms.double(0.5),
    invert = cms.bool(False),
    filter = cms.bool(False)                                                          
)
process.testSVfitTrackLikelihoodProductionSequence += process.ak5PFJetsNotOverlappingWithLeptons

from JetMETCorrections.Type1MET.pfMETCorrections_cff import pfCandsNotInJet
process.pfCandsNotInJetForPFMEtSignCovMatrix = pfCandsNotInJet.clone()
process.testSVfitTrackLikelihoodProductionSequence += process.pfCandsNotInJetForPFMEtSignCovMatrix

from RecoMET.METProducers.METSigParams_cfi import *
process.pfMEtSignCovMatrix = cms.EDProducer("PFMEtSignCovMatrixProducer",
    METSignificance_params,                     
    src = cms.VInputTag(
        'genMatchedPatElectrons',                                                
        'genMatchedPatMuons',
        'genMatchedPatTaus',
        'ak5PFJetsNotOverlappingWithLeptons',                                        
        'pfCandsNotInJetForPFMEtSignCovMatrix'
    ),
    addJERcorr = cms.PSet(
        inputFileName = cms.FileInPath('PhysicsTools/PatUtils/data/pfJetResolutionMCtoDataCorrLUT.root'),
        lutName = cms.string('pfJetResolutionMCtoDataCorrLUT')
    )
)
process.testSVfitTrackLikelihoodProductionSequence += process.pfMEtSignCovMatrix

process.load("JetMETCorrections/METPUSubtraction/mvaPFMET_cff")
process.pfMEtMVA.srcLeptons = cms.VInputTag(
    'genMatchedPatElectrons',                                                
    'genMatchedPatMuons',
    'genMatchedPatTaus'
)
process.testSVfitTrackLikelihoodProductionSequence += process.pfMEtMVAsequence

process.pfMEtMVAwithPFMEtCov = cms.EDProducer("RecoPFMEtSignCovMatrixEmbedder",
    src = cms.InputTag('pfMEtMVA'),
    srcCov = cms.InputTag('pfMEtSignCovMatrix'),
    sf = cms.double(0.60*0.60)
)
process.testSVfitTrackLikelihoodProductionSequence += process.pfMEtMVAwithPFMEtCov

process.pfMEtMVAunityResponse = process.pfMEtMVA.clone(
    inputFileNames = cms.PSet(
        DPhi = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmetphi_53_June2013_type1.root'), # CV: same for unity and non-unity response training
        CovU2 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru2cov_53_Dec2012.root'),
        U = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmet_53_June2013_type1_UnityResponse.root'),
        CovU1 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru1cov_53_Dec2012.root')
    )
)
process.testSVfitTrackLikelihoodProductionSequence += process.pfMEtMVAunityResponse

process.pfMEtMVAunityResponseWithPFMEtCov = cms.EDProducer("RecoPFMEtSignCovMatrixEmbedder",
    src = cms.InputTag('pfMEtMVAunityResponse'),
    srcCov = cms.InputTag('pfMEtSignCovMatrix'),
    sf = cms.double(0.75*0.75)
)
process.testSVfitTrackLikelihoodProductionSequence += process.pfMEtMVAunityResponseWithPFMEtCov
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# select primary event vertex
process.load("TauAnalysis/RecoTools/recoVertexSelectionByLeptonTracks_cff")
process.selectedPrimaryVertexQuality.src = cms.InputTag('offlinePrimaryVerticesWithBS')
process.selectedPrimaryVertexByLeptonMatch.srcLeptons = cms.VInputTag(
    'genMatchedPatElectrons',    
    'genMatchedPatMuons',
    'genMatchedPatTaus'
)
process.testSVfitTrackLikelihoodProductionSequence += process.selectPrimaryVertexByLeptonTracks

# require event to have exactly one vertex associated to tracks of tau decay products
process.recEventVertexFilter = cms.EDFilter("VertexCountFilter",
    src = cms.InputTag('selectedPrimaryVertexByLeptonMatch'),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1)                                            
)
process.testSVfitTrackLikelihoodProductionSequence += process.recEventVertexFilter

process.genEventVertex = cms.EDProducer("GenVertexProducer",
    srcGenParticles = cms.InputTag(genTaus),
    pdgIds = cms.vint32(15)
)
process.testSVfitTrackLikelihoodProductionSequence += process.genEventVertex
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# dump tau-jet collection

process.dumpPATTaus = cms.EDAnalyzer("DumpPATTausForRaman",
    src = cms.InputTag('genMatchedPatTaus'),
    srcGenTauJets = cms.InputTag('tauGenJetsSelectorAllHadrons'),                                 
    srcVertex = cms.InputTag('selectedPrimaryVertexByLeptonMatch'),
    srcTracks = cms.InputTag('generalTracks'),
    minPt = cms.double(0.)                                     
)
##process.testSVfitTrackLikelihoodProductionSequence += process.dumpPATTaus
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# run SVfit with and without Track likelihoods

process.load("TauAnalysis.CandidateTools.svFitAlgorithmDiTau_cfi")

# CV: fix tau decay parameters to Monte Carlo truth values
process.svFitTauToElecBuilder.fixToGenVisEnFracX = cms.bool(False)
process.svFitTauToElecBuilder.fixToGenPhiLab = cms.bool(False)
process.svFitTauToElecBuilder.fixToGenNuInvMass = cms.bool(False)
process.svFitTauToElecBuilder.fixRecToGenVertex = cms.bool(False)
process.svFitTauToElecBuilder.fixToGenDecayDistance = cms.bool(False)
process.svFitTauToElecBuilder.fixToGenVisP4 = cms.bool(False)
process.svFitTauToElecBuilder.initializeToGen = cms.bool(False)
#process.svFitTauToElecBuilder.srcGenTaus = cms.InputTag('genParticles')
process.svFitTauToElecBuilder.dRmatch = cms.double(0.3)
process.svFitTauToElecBuilder.verbosity = cms.int32(0)

process.svFitTauToMuBuilder.fixToGenVisEnFracX = cms.bool(False)
process.svFitTauToMuBuilder.fixToGenPhiLab = cms.bool(False)
process.svFitTauToMuBuilder.fixToGenNuInvMass = cms.bool(False)
process.svFitTauToMuBuilder.fixRecToGenVertex = cms.bool(False)
process.svFitTauToMuBuilder.fixToGenDecayDistance = cms.bool(False)
process.svFitTauToMuBuilder.fixToGenVisP4 = cms.bool(False)
process.svFitTauToMuBuilder.initializeToGen = cms.bool(False)
#process.svFitTauToMuBuilder.srcGenTaus = cms.InputTag('genParticles')
process.svFitTauToMuBuilder.dRmatch = cms.double(0.3)
process.svFitTauToMuBuilder.verbosity = cms.int32(0)

process.svFitTauToHadBuilder.fixToGenVisEnFracX = cms.bool(False)
process.svFitTauToHadBuilder.fixToGenPhiLab = cms.bool(False)
process.svFitTauToHadBuilder.fixRecToGenVertex = cms.bool(False)
process.svFitTauToHadBuilder.fixToGenDecayDistance = cms.bool(False)
process.svFitTauToHadBuilder.fixToGenVisP4 = cms.bool(False)
process.svFitTauToHadBuilder.initializeToGen = cms.bool(False)
#process.svFitTauToHadBuilder.srcGenTaus = cms.InputTag('genParticles')
process.svFitTauToHadBuilder.dRmatch = cms.double(0.3)
process.svFitTauToHadBuilder.verbosity = cms.int32(0)

# CV: fix event vertex position to Monte Carlo truth value
process.svFitEventBuilder.fixToGenVertex = cms.bool(False)
#process.svFitEventBuilder.srcGenVertex = cms.InputTag('genEventVertex')
process.svFitEventBuilder.verbosity = cms.int32(0)

process.svFitElectronLikelihoodTrackInfo = cms.PSet(
    pluginName = cms.string("svFitTauToElecLikelihoodTrackInfo"),
    pluginType = cms.string("SVfitTauDecayLikelihoodTrackInfo"),
    useLifetimeConstraint = cms.bool(True),
    sfProdVertexCov = cms.double(1.0),
    sfDecayVertexCov = cms.double(1.0),
    verbosity = cms.int32(0)  
)

process.svFitMuonLikelihoodTrackInfo = cms.PSet(
    pluginName = cms.string("svFitTauToMuLikelihoodTrackInfo"),
    pluginType = cms.string("SVfitTauDecayLikelihoodTrackInfo"),
    useLifetimeConstraint = cms.bool(True),
    sfProdVertexCov = cms.double(1.0),
    sfDecayVertexCov = cms.double(1.0),
    verbosity = cms.int32(0)  
)

process.svFitTauLikelihoodTrackInfo = cms.PSet(
    pluginName = cms.string("svFitTauToHadLikelihoodTrackInfo"),
    pluginType = cms.string("SVfitTauDecayLikelihoodTrackInfo"),
    useLifetimeConstraint = cms.bool(False),
    sfProdVertexCov = cms.double(1.0),
    sfDecayVertexCov = cms.double(1.0),
    verbosity = cms.int32(0)  
)

process.svFitEventLikelihoodTrackInfo = cms.PSet(
    pluginName = cms.string("svFitEventLikelihoodTrackInfo"),
    pluginType = cms.string("SVfitEventLikelihoodTrackInfo"),
    verbosity = cms.int32(0)
)

srcGenLeg1                           = None
srcRecLeg1                           = None
svFitLikelihoodLeg1_kinematics      = None
svFitLikelihoodLeg1_dummyKinematics = None
svFitBuilderLeg1                    = None
svFitStandAloneTypeLeg1             = None
srcGenLeg2                           = None
srcRecLeg2                           = None
svFitLikelihoodLeg2_kinematics      = None
svFitLikelihoodLeg2_dummyKinematics = None
svFitBuilderLeg2                    = None
svFitStandAloneTypeLeg2             = None
if channel == 'mutau':
    srcGenLeg1                           = genMuonsFromTauDecays
    srcRecLeg1                           = 'genMatchedPatMuons'
    svFitLikelihoodLeg1_kinematics      = process.svFitMuonLikelihoodMatrixElement
    svFitLikelihoodLeg1_dummyKinematics = process.svFitMuonLikelihoodDummy
    svFitLikelihoodLeg1_trackinfo       = process.svFitMuonLikelihoodTrackInfo
    svFitBuilderLeg1                    = process.svFitTauToMuBuilder
    svFitStandAloneTypeLeg1             = "lep"
    srcGenLeg2                           = genTauJetsFromTauDecays
    srcRecLeg2                           = 'genMatchedPatTaus'
    svFitLikelihoodLeg2_kinematics      = process.svFitTauLikelihoodPhaseSpace
    svFitLikelihoodLeg2_dummyKinematics = process.svFitTauLikelihoodDummy
    svFitLikelihoodLeg2_trackinfo       = process.svFitTauLikelihoodTrackInfo
    svFitBuilderLeg2                    = process.svFitTauToHadBuilder
    svFitStandAloneTypeLeg2             = "had"
elif channel == 'etau':
    srcGenLeg1                           = genElectronsFromTauDecays
    srcRecLeg1                           = 'genMatchedPatElectrons'
    svFitLikelihoodLeg1_kinematics      = process.svFitElectronLikelihoodMatrixElement
    svFitLikelihoodLeg1_dummyKinematics = process.svFitElectronLikelihoodDummy
    svFitLikelihoodLeg1_trackinfo       = process.svFitElectronLikelihoodTrackInfo
    svFitBuilderLeg1                    = process.svFitTauToElecBuilder
    svFitStandAloneTypeLeg1             = "lep"
    srcGenLeg2                           = genTauJetsFromTauDecays
    srcRecLeg2                           = 'genMatchedPatTaus'
    svFitLikelihoodLeg2_kinematics      = process.svFitTauLikelihoodPhaseSpace
    svFitLikelihoodLeg2_dummyKinematics = process.svFitTauLikelihoodDummy
    svFitLikelihoodLeg2_trackinfo       = process.svFitTauLikelihoodTrackInfo
    svFitBuilderLeg2                    = process.svFitTauToHadBuilder
    svFitStandAloneTypeLeg2             = "had"
elif channel == 'emu':
    srcGenLeg1                           = genElectronsFromTauDecays
    srcRecLeg1                           = 'genMatchedPatElectrons'
    svFitLikelihoodLeg1_kinematics      = process.svFitElectronLikelihoodMatrixElement
    svFitLikelihoodLeg1_dummyKinematics = process.svFitElectronLikelihoodDummy
    svFitLikelihoodLeg1_trackinfo       = process.svFitElectronLikelihoodTrackInfo
    svFitBuilderLeg1                    = process.svFitTauToElecBuilder
    svFitStandAloneTypeLeg1             = "lep"
    srcGenLeg2                           = genMuonsFromTauDecays
    srcRecLeg2                           = 'genMatchedPatMuons'
    svFitLikelihoodLeg2_kinematics      = process.svFitMuonLikelihoodMatrixElement
    svFitLikelihoodLeg2_dummyKinematics = process.svFitMuonLikelihoodDummy
    svFitLikelihoodLeg2_trackinfo       = process.svFitMuonLikelihoodTrackInfo
    svFitBuilderLeg2                    = process.svFitTauToMuBuilder
    svFitStandAloneTypeLeg2             = "lep"
elif channel == 'tautau':
    srcGenLeg1                           = genTauJetsFromTauDecays
    srcRecLeg1                           = 'genMatchedPatTaus'
    svFitLikelihoodLeg1_kinematics      = process.svFitTauLikelihoodPhaseSpace
    svFitLikelihoodLeg1_dummyKinematics = process.svFitTauLikelihoodDummy
    svFitLikelihoodLeg1_trackinfo       = process.svFitTauLikelihoodTrackInfo
    svFitBuilderLeg1                    = process.svFitTauToHadBuilder
    svFitStandAloneTypeLeg1             = "had"
    srcGenLeg2                           = genTauJetsFromTauDecays
    srcRecLeg2                           = 'genMatchedPatTaus'
    svFitLikelihoodLeg2_kinematics      = process.svFitTauLikelihoodPhaseSpace
    svFitLikelihoodLeg2_dummyKinematics = process.svFitTauLikelihoodDummy
    svFitLikelihoodLeg2_trackinfo       = process.svFitTauLikelihoodTrackInfo
    svFitBuilderLeg2                    = process.svFitTauToHadBuilder
    svFitStandAloneTypeLeg2             = "had"
else:
    raise ValueError("Invalid channel = %s !!" % channel)

if svFitLikelihoodLeg1_kinematics:
    svFitLikelihoodLeg1_kinematics.verbosity = cms.int32(0)
if svFitLikelihoodLeg1_dummyKinematics:
    svFitLikelihoodLeg1_dummyKinematics.verbosity = cms.int32(0)
if svFitLikelihoodLeg2_kinematics:
    svFitLikelihoodLeg2_kinematics.verbosity = cms.int32(0)
if svFitLikelihoodLeg2_dummyKinematics:
    svFitLikelihoodLeg2_dummyKinematics.verbosity = cms.int32(0)
    
svFitProducerModuleNames = dict()
svFitAnalyzerModuleTypes = dict()

def customCapitalize(name):
    retVal = name[0:1].capitalize() + name[1:]
    return retVal

likelihoodOptions_base = [
    'decayKine',
    'decayKinePlusMEt',
    'trackWlifetimeWdecayVertexForLeg1Leg2PlusDecayKinePlusMEt'
]
likelihoodOptions = []
for likelihoodOption_base in likelihoodOptions_base:
    for logM_power in [ 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0 ]:
        logM_string = "logM%1.1f" % logM_power
        logM_string = logM_string.replace(".", "_")
        likelihoodOption = likelihoodOption_base + customCapitalize(logM_string)
        print "likelihoodOption = %s" % likelihoodOption
        likelihoodOptions.append(likelihoodOption)

svFitProducersVEGAS       = {} # key = integrationOption
svFitProducersMarkovChain = {} # key = integrationOption

for likelihoodOption in likelihoodOptions:
    for integrationOption in [ 'mean5sigmaWithinMax' ]:
        for applyJacobiFactors in [ True, False ]:
            for fixRecToGenVertex in [ True, False ]:

                applyJacobiFactors_string = None
                if applyJacobiFactors:
                    applyJacobiFactors_string = "WjacobiFactors"
                else:
                    applyJacobiFactors_string = "WOjacobiFactors"

                fixRecToGenVertex_string = None
                if fixRecToGenVertex:
                    fixRecToGenVertex_string = "genVertex"
                else:
                    fixRecToGenVertex_string = "recVertex"

                print "processing likelihood option = %s, integration option = %s, vertex = %s" % (likelihoodOption, integrationOption, fixRecToGenVertex_string)

                if not integrationOption in svFitProducersVEGAS.keys():
                    svFitProducersVEGAS[integrationOption] = {}
                if not fixRecToGenVertex_string in svFitProducersVEGAS[integrationOption].keys():
                    svFitProducersVEGAS[integrationOption][fixRecToGenVertex_string] = []
                if not integrationOption in svFitProducersMarkovChain.keys():
                    svFitProducersMarkovChain[integrationOption] = {}
                if not fixRecToGenVertex_string in svFitProducersMarkovChain[integrationOption].keys():
                    svFitProducersMarkovChain[integrationOption][fixRecToGenVertex_string] = []
        
                leg1Builder = copy.deepcopy(svFitBuilderLeg1)
                leg2Builder = copy.deepcopy(svFitBuilderLeg2)
                eventBuilder = copy.deepcopy(process.svFitEventBuilder)
                leg1LikelihoodFunctions = []
                leg2LikelihoodFunctions = []
                resonanceLikelihoodFunctions = []
                eventLikelihoodFunctions = []
                if likelihoodOption.find("track") != -1 or likelihoodOption.find("Track") != -1:
                    print " - enabling track likelihood"
                    leg1LikelihoodTrack = copy.deepcopy(svFitLikelihoodLeg1_trackinfo)
                    leg1LikelihoodTrack.ignore3Prongs = cms.bool(False)
                    leg1LikelihoodTrack.ignore1Prongs = cms.bool(True)
                    leg1LikelihoodFunctions.append(leg1LikelihoodTrack)
                    leg2LikelihoodTrack = copy.deepcopy(svFitLikelihoodLeg2_trackinfo)
                    leg2LikelihoodTrack.ignore3Prongs = cms.bool(False)
                    leg2LikelihoodTrack.ignore1Prongs = cms.bool(True)
                    leg2LikelihoodFunctions.append(leg2LikelihoodTrack)
                    ## CV: varying primary event vertex position makes algorithm numerically unstable !!
                    ##eventLikelihoodFunctions.append(process.svFitEventLikelihoodTrackInfo)
                    if likelihoodOption.find("Wlifetime") != -1:
                        print "  + enabling lifetime constraint"
                        leg1LikelihoodTrack.useLifetimeConstraint = cms.bool(True)
                        leg2LikelihoodTrack.useLifetimeConstraint = cms.bool(True)
                    elif likelihoodOption.find("WOlifetime") != -1:
                        print "  + disabling lifetime constraint"
                        leg1LikelihoodTrack.useLifetimeConstraint = cms.bool(False)
                        leg2LikelihoodTrack.useLifetimeConstraint = cms.bool(False)
                    else:
                        raise ValueError("Invalid likelihood option = %s !!" % likelihoodOption)
                    if likelihoodOption.find("WdecayVertexForLeg1Leg2") != -1:
                        print "  + enabling decay vertex reconstruction for leg1 and leg2"
                        leg1Builder.fitDecayVertex = cms.bool(True)
                        leg2Builder.fitDecayVertex = cms.bool(True)
                    elif likelihoodOption.find("WdecayVertexForLeg1") != -1:
                        print "  + enabling decay vertex reconstruction for leg1"
                        leg1Builder.fitDecayVertex = cms.bool(True)
                        leg2Builder.fitDecayVertex = cms.bool(False)
                    elif likelihoodOption.find("WdecayVertexForLeg2") != -1:
                        print "  + enabling decay vertex reconstruction for leg2"
                        leg1Builder.fitDecayVertex = cms.bool(False)
                        leg2Builder.fitDecayVertex = cms.bool(True)
                    else:
                        print "  + disabling decay vertex reconstruction"
                        leg1Builder.fitDecayVertex = cms.bool(False)
                        leg2Builder.fitDecayVertex = cms.bool(False)
                if likelihoodOption.find("decayKine") != -1 or likelihoodOption.find("DecayKine") != -1:
                    print " - enabling likelihood for tau decay kinematics"
                    leg1LikelihoodDecayKine = copy.deepcopy(svFitLikelihoodLeg1_kinematics)
                    leg1LikelihoodDecayKine.applySinThetaFactor = cms.bool(False)
                    leg1LikelihoodFunctions.append(leg1LikelihoodDecayKine)
                    leg2LikelihoodDecayKine = copy.deepcopy(svFitLikelihoodLeg2_kinematics)
                    leg2LikelihoodDecayKine.applySinThetaFactor = cms.bool(False)
                    leg2LikelihoodFunctions.append(leg2LikelihoodDecayKine)
                else:
                    leg1LikelihoodDecayKine = copy.deepcopy(svFitLikelihoodLeg1_dummyKinematics)
                    leg1LikelihoodDecayKine.verbosity = cms.int32(0)
                    leg1LikelihoodFunctions.append(leg1LikelihoodDecayKine)
                    leg2LikelihoodDecayKine = copy.deepcopy(svFitLikelihoodLeg2_dummyKinematics)
                    leg2LikelihoodDecayKine.verbosity = cms.int32(0)
                    leg2LikelihoodFunctions.append(leg2LikelihoodDecayKine)
                if likelihoodOption.find("met") != -1 or likelihoodOption.find("MEt") != -1:
                    print " - enabling MEt likelihood"
                    eventLikelihoodMEt = copy.deepcopy(process.svFitEventLikelihoodMEt2)
                    eventLikelihoodMEt.verbosity = cms.int32(0)
                    eventLikelihoodFunctions.append(eventLikelihoodMEt)
                if fixRecToGenVertex:
                    leg1Builder.fixRecToGenVertex = cms.bool(True)
                    leg1Builder.srcGenTaus = cms.InputTag('genParticles')
                    leg2Builder.fixRecToGenVertex = cms.bool(True)
                    leg2Builder.srcGenTaus = cms.InputTag('genParticles')
                    eventBuilder.srcGenVertex = cms.InputTag('genEventVertex')
                    eventBuilder.fixToGenVertex = cms.bool(True)
                else:
                    leg1Builder.fixRecToGenVertex = cms.bool(False)
                    leg2Builder.fixRecToGenVertex = cms.bool(False)
                    eventBuilder.fixToGenVertex = cms.bool(False)
                logM_string = None
                if likelihoodOption.find("logM") != -1:
                    logM_string = "logM"
                elif likelihoodOption.find("LogM") != -1:
                    logM_string = "LogM"
                if logM_string:
                    resonanceLikelihoodLogM = copy.deepcopy(process.svFitResonanceLikelihoodLogM)                    
                    power_string = likelihoodOption[likelihoodOption.find(logM_string) + 4:likelihoodOption.find(logM_string) + 7]
                    power_string = power_string.replace("_", ".")
                    power = float(power_string)
                    print "power_string = %s: power = %1.1f" % (power_string, power)
                    resonanceLikelihoodLogM.power = cms.double(power)
                    resonanceLikelihoodFunctions.append(resonanceLikelihoodLogM)

                # VEGAS integration of likelihood functions
                svFitProducerModuleVEGAS = process.svFitProducerByIntegration.clone()
                svFitProducerModuleVEGAS.config.event.resonances.A.daughters.leg1.src = cms.InputTag(srcRecLeg1)
                svFitProducerModuleVEGAS.config.event.resonances.A.daughters.leg1.builder = leg1Builder
                svFitProducerModuleVEGAS.config.event.resonances.A.daughters.leg1.likelihoodFunctions = cms.VPSet(leg1LikelihoodFunctions)
                svFitProducerModuleVEGAS.config.event.resonances.A.daughters.leg2.src = cms.InputTag(srcRecLeg2)
                svFitProducerModuleVEGAS.config.event.resonances.A.daughters.leg2.builder = leg2Builder
                svFitProducerModuleVEGAS.config.event.resonances.A.daughters.leg2.likelihoodFunctions = cms.VPSet(leg2LikelihoodFunctions)
                svFitProducerModuleVEGAS.config.event.resonances.A.likelihoodFunctions = cms.VPSet(resonanceLikelihoodFunctions)
                svFitProducerModuleVEGAS.config.event.builder = eventBuilder
                svFitProducerModuleVEGAS.config.event.likelihoodFunctions = cms.VPSet(eventLikelihoodFunctions)
                svFitProducerModuleVEGAS.config.event.srcMEt = cms.InputTag('pfType1CorrectedMet')
                svFitProducerModuleVEGAS.config.event.srcPrimaryVertex = cms.InputTag('selectedPrimaryVertexByLeptonMatch')
                svFitProducerModuleVEGAS.algorithm.max_or_median = cms.string(integrationOption)
                pluginNameVEGAS = "svFitProducerByIntegrationW%s%s" % (likelihoodOption, integrationOption)
                svFitProducerModuleVEGAS.algorithm.pluginName = cms.string(pluginNameVEGAS)
                svFitProducerModuleVEGAS.algorithm.verbosity = cms.int32(1)
                svFitProducerModuleNameVEGAS = "svFitProducerByIntegrationW%s%s%s" % (likelihoodOption, customCapitalize(integrationOption), customCapitalize(fixRecToGenVertex_string))
                svFitAnalyzerModuleTypeVEGAS = "SVfitEventHypothesisByIntegrationAnalyzer"
                print "adding module '%s' of type VEGAS" % svFitProducerModuleNameVEGAS
                svFitProducersVEGAS[integrationOption][fixRecToGenVertex_string].append(svFitProducerModuleNameVEGAS)
    
                svFitProducerModuleNames["%s%s%sVEGAS" % (likelihoodOption, customCapitalize(integrationOption), customCapitalize(fixRecToGenVertex_string))] = svFitProducerModuleNameVEGAS
                svFitAnalyzerModuleTypes["%s%s%sVEGAS" % (likelihoodOption, customCapitalize(integrationOption), customCapitalize(fixRecToGenVertex_string))] = svFitAnalyzerModuleTypeVEGAS

                setattr(process, svFitProducerModuleNameVEGAS, svFitProducerModuleVEGAS)
                process.testSVfitTrackLikelihoodProductionSequence += svFitProducerModuleVEGAS

                # Markov Chain integration of likelihood functions
                numCallsMarkovChain = 250
                ##numCallsMarkovChain = 1
                svFitProducerModuleMarkovChain = process.svFitProducerByIntegration2.clone()
                svFitProducerModuleMarkovChain.config.event = svFitProducerModuleVEGAS.config.event
                svFitProducerModuleMarkovChain.algorithm = process.svFitProducerByIntegration2.algorithm.clone()
                svFitProducerModuleMarkovChain.algorithm.markovChainOptions.initMode = cms.string("none")
                svFitProducerModuleMarkovChain.algorithm.markovChainOptions.numIterBurnin = cms.uint32(numCallsMarkovChain*100)
                svFitProducerModuleMarkovChain.algorithm.markovChainOptions.numIterSampling = cms.uint32(numCallsMarkovChain*1000)
                svFitProducerModuleMarkovChain.algorithm.markovChainOptions.numIterSimAnnealingPhase1 = cms.uint32(numCallsMarkovChain*20)
                svFitProducerModuleMarkovChain.algorithm.markovChainOptions.numIterSimAnnealingPhase2 = cms.uint32(numCallsMarkovChain*60)
                svFitProducerModuleMarkovChain.algorithm.markovChainOptions.alpha = cms.double(1.0 - 0.1/numCallsMarkovChain)
                svFitProducerModuleMarkovChain.algorithm.markovChainOptions.numChains = cms.uint32(1)
                svFitProducerModuleMarkovChain.algorithm.markovChainOptions.numBatches = cms.uint32(1)
                svFitProducerModuleMarkovChain.algorithm.markovChainOptions.epsilon0 = cms.double(1.e-2)    
                svFitProducerModuleMarkovChain.algorithm.monitorFilePath = cms.string('/data1/veelken/tmp/svFitStudies/')
                svFitProducerModuleMarkovChain.algorithm.max_or_median = cms.string(integrationOption)
                pluginNameMarkovChain = "svFitProducerByIntegration2W%s%s%s" % (likelihoodOption, customCapitalize(integrationOption), customCapitalize(fixRecToGenVertex_string))
                svFitProducerModuleMarkovChain.algorithm.pluginName = cms.string(pluginNameMarkovChain)
                svFitProducerModuleMarkovChain.algorithm.monitorMarkovChain = cms.bool(False)
                svFitProducerModuleMarkovChain.algorithm.verbosity = cms.int32(1)
                svFitProducerModuleNameMarkovChain = "svFitProducerByIntegration2W%s%s%s" % (likelihoodOption, customCapitalize(integrationOption), customCapitalize(fixRecToGenVertex_string))
                svFitAnalyzerModuleTypeMarkovChain = "SVfitEventHypothesisAnalyzer"
                print "adding module '%s' of type Markov-Chain" % svFitProducerModuleNameMarkovChain
                svFitProducersMarkovChain[integrationOption][fixRecToGenVertex_string].append(svFitProducerModuleNameMarkovChain)

                svFitProducerModuleNames["%s%s%sMarkovChain" % (likelihoodOption, customCapitalize(integrationOption), customCapitalize(fixRecToGenVertex_string))] = svFitProducerModuleNameMarkovChain
                svFitAnalyzerModuleTypes["%s%s%sMarkovChain" % (likelihoodOption, customCapitalize(integrationOption), customCapitalize(fixRecToGenVertex_string))] = svFitAnalyzerModuleTypeMarkovChain

                setattr(process, svFitProducerModuleNameMarkovChain, svFitProducerModuleMarkovChain)
                process.testSVfitTrackLikelihoodProductionSequence += svFitProducerModuleMarkovChain

##print "svFitProducerModuleNames:"
##print svFitProducerModuleNames

# CV: select among different likelihoodOptions
#      the SVfit mass with the best resolution (minimal sigma)
for integrationOption in svFitProducersVEGAS.keys():
    for fixRecToGenVertex_string in svFitProducersVEGAS[integrationOption].keys():
        svFitSelectorModuleVEGAS = cms.EDProducer("SVfitBestEventHypothesisByIntegrationSelector",
            src = cms.VInputTag(svFitProducersVEGAS[integrationOption][fixRecToGenVertex_string]),
            selection = cms.string("resonance(0).isValidSolution"),
            # CV: take width of likelihood function around minimum as measure of resolution,
            #     i.e. sum sigmaUp and sigmaDown linearly
            rank = cms.string("-(resonance(0).massErrUp() + resonance(0).massErrDown())"),
            ##rank = cms.string("-(resonance(0).massErrUp()*resonance(0).massErrUp() + resonance(0).massErrDown()*resonance(0).massErrDown())")
           verbosity = cms.int32(1)                                               
        )
        svFitSelectorModuleNameVEGAS = "svFitSelectorByIntegration%s%s%s" % ("MinSigma", customCapitalize(integrationOption), customCapitalize(fixRecToGenVertex_string)) 
        svFitAnalyzerModuleTypeVEGAS = "SVfitEventHypothesisByIntegrationAnalyzer"
        print "adding module '%s' of type VEGAS" % svFitSelectorModuleNameVEGAS
        svFitProducersVEGAS[integrationOption][fixRecToGenVertex_string].append(svFitSelectorModuleNameVEGAS)
    
        svFitProducerModuleNames["%s%s%sVEGAS" % ("MinSigma", customCapitalize(integrationOption), customCapitalize(fixRecToGenVertex_string))] = svFitSelectorModuleNameVEGAS
        svFitAnalyzerModuleTypes["%s%s%sVEGAS" % ("MinSigma", customCapitalize(integrationOption), customCapitalize(fixRecToGenVertex_string))] = svFitAnalyzerModuleTypeVEGAS

        setattr(process, svFitSelectorModuleNameVEGAS, svFitSelectorModuleVEGAS)
        process.testSVfitTrackLikelihoodProductionSequence += svFitSelectorModuleVEGAS
    
        svFitSelectorModuleMarkovChain = cms.EDProducer("SVfitBestEventHypothesisSelector",
            src = cms.VInputTag(svFitProducersMarkovChain[integrationOption][fixRecToGenVertex_string]),
            selection = svFitSelectorModuleVEGAS.selection,
            rank = svFitSelectorModuleVEGAS.rank,
            verbosity = cms.int32(1)  
        )
        svFitSelectorModuleNameMarkovChain = "svFitSelectorByIntegration2%s%s%s" % ("MinSigma", customCapitalize(integrationOption), customCapitalize(fixRecToGenVertex_string))
        svFitAnalyzerModuleTypeMarkovChain = "SVfitEventHypothesisAnalyzer"
        print "adding module '%s' of type Markov-Chain" % svFitSelectorModuleNameMarkovChain
        svFitProducersMarkovChain[integrationOption][fixRecToGenVertex_string].append(svFitSelectorModuleNameMarkovChain)
    
        svFitProducerModuleNames["%s%s%sMarkovChain" % ("MinSigma", customCapitalize(integrationOption), customCapitalize(fixRecToGenVertex_string))] = svFitSelectorModuleNameMarkovChain
        svFitAnalyzerModuleTypes["%s%s%sMarkovChain" % ("MinSigma", customCapitalize(integrationOption), customCapitalize(fixRecToGenVertex_string))] = svFitAnalyzerModuleTypeMarkovChain

        setattr(process, svFitSelectorModuleNameMarkovChain, svFitSelectorModuleMarkovChain)
        process.testSVfitTrackLikelihoodProductionSequence += svFitSelectorModuleMarkovChain
#--------------------------------------------------------------------------------

process.testSVfitTrackLikelihoodAnalysisSequence = cms.Sequence()

#--------------------------------------------------------------------------------
# fill histograms of mass reconstructed by running SVfit with different integration options

for option in svFitAnalyzerModuleTypes.keys():

    svFitProducerModuleName = svFitProducerModuleNames[option]
    svFitAnalyzerModuleType = svFitAnalyzerModuleTypes[option]

    if not hasattr(process, svFitProducerModuleName):
        continue
    
    svFitAnalyzerModule = cms.EDAnalyzer(svFitAnalyzerModuleType,
        srcEventHypotheses = cms.InputTag(svFitProducerModuleName),
        srcGenTauPairs = cms.InputTag(genTauPairs),
        srcGenLeg1 = cms.InputTag(srcGenLeg1),
        srcGenLeg2 = cms.InputTag(srcGenLeg2),
        srcGenMEt = cms.InputTag('genMetFromGenParticles'),
        srcPFMEtCovMatrix = cms.InputTag('pfMEtSignCovMatrix'),
        srcWeights = cms.VInputTag(),
        idxResonance = cms.int32(0),
        numBinsSVfitMass = cms.int32(1000),
        svFitMassMax = cms.double(1000.),
        numBinsSVfitSigma = cms.int32(250),
        svFitSigmaMax = cms.double(250.),
        dqmDirectory = cms.string(svFitProducerModuleName)
    )
    svFitAnalyzerModuleName = svFitProducerModuleName
    svFitAnalyzerModuleName = svFitAnalyzerModuleName.replace("Producer", "Analyzer")
    svFitAnalyzerModuleName = svFitAnalyzerModuleName.replace("Selector", "Analyzer")
    setattr(process, svFitAnalyzerModuleName, svFitAnalyzerModule)
    process.testSVfitTrackLikelihoodAnalysisSequence += svFitAnalyzerModule
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# make plots of SVfit input variables, separetely for 
# events in which SVfit mass reconstructed using track information is better/worse
# than SVfit mass reconstructed without using track information
# in likelihood model

for idxOption1, option1 in enumerate(likelihoodOptions):
    for idxOption2, option2 in enumerate(likelihoodOptions):
                    
        if not (idxOption2 > idxOption1):
            continue

        for integrationMethod in [ 'VEGAS', 'MarkovChain' ]:
            for integrationOption in [ 'mean5sigmaWithinMax' ]:
                for fixRecToGenVertex in [ True, False ]:
                    
                    fixRecToGenVertex_string = None
                    if fixRecToGenVertex:
                        fixRecToGenVertex_string = "genVertex"
                    else:
                        fixRecToGenVertex_string = "recVertex"
                    

            ##print "option1 = %s, option2 = %s, intMethod = %s" % (option1, option2, intMethod)

            option1Name = "%s%s%s%s" % (option1, customCapitalize(integrationOption), customCapitalize(fixRecToGenVertex_string), integrationMethod)
            if not (option1Name in svFitProducerModuleNames):
                print "Failed to find option = %s in svFitProducerModuleNames --> skipping !!" % option1Name
                continue
            svFitProducerModuleName1 = svFitProducerModuleNames[option1Name]
            ##print "svFitProducerModuleName1['%s'] = %s" % (option1Name, svFitProducerModuleName1)

            option2Name = "%s%s%s%s" % (option2, customCapitalize(integrationOption), customCapitalize(fixRecToGenVertex_string), integrationMethod)
            if not (option2Name in svFitProducerModuleNames):
                print "Failed to find option = %s in svFitProducerModuleNames --> skipping !!" % option2Name
                continue
            svFitProducerModuleName2 = svFitProducerModuleNames[option2Name]
            ##print "svFitProducerModuleName2['%s'] = %s" % (option2Name, svFitProducerModuleName2)

            if not (hasattr(process, svFitProducerModuleName1) and hasattr(process, svFitProducerModuleName2)):
                continue
            
            analyzerCorrelationOption1vs2 = cms.EDAnalyzer("SVfitEventHypothesisCorrelationAnalyzer",
                srcEventHypotheses1 = cms.InputTag(svFitProducerModuleName1),
                srcEventHypotheses2 = cms.InputTag(svFitProducerModuleName2),
                srcWeights = cms.VInputTag(),
                dqmDirectory = cms.string("correlation%sVs%s" % (svFitProducerModuleName1, svFitProducerModuleName2)) 
            )
            analyzerNameCorrelationOption1vs2 = "correlation%sVs%s" % (svFitProducerModuleName1, svFitProducerModuleName2)
            print "adding module '%s' of type SVfitEventHypothesisCorrelationAnalyzer" % analyzerNameCorrelationOption1vs2
            setattr(process, analyzerNameCorrelationOption1vs2, analyzerCorrelationOption1vs2)
            process.testSVfitTrackLikelihoodAnalysisSequence += analyzerCorrelationOption1vs2

            watchdogOption1vs2 = cms.EDAnalyzer("SVfitEventHypothesisWatchdog",
                srcEventHypotheses1 = cms.InputTag(svFitProducerModuleName1),
                srcEventHypotheses2 = cms.InputTag(svFitProducerModuleName2),
                srcGenTauPairs = cms.InputTag(genTauPairs),
                # CV: only report differences in case adding information
                #     actually makes the mass resolution worse
                #    (Note: this assumes that the likelihoodOptions are sorted in the order [ 'decayKine', 'decayKinePlusMEt', 'trackWlifetimeWdecayVertexForLeg1Leg2PlusDecayKinePlusMEt' ])
                threshold1MassUp = cms.double(1.e+3),
                threshold1MassDown = cms.double(1.e+3),
                threshold1Sigma = cms.double(1.e+3),
                threshold2MassUp = cms.double(10.),
                threshold2MassDown = cms.double(10.),
                threshold2Sigma = cms.double(1.e+3)
            )
            watchdogNameOption1vs2 = "watchdog%sVs%s" % (svFitProducerModuleName1, svFitProducerModuleName2)
            print "adding watchdog '%s'" % watchdogNameOption1vs2
            setattr(process, watchdogNameOption1vs2, watchdogOption1vs2)
            process.testSVfitTrackLikelihoodAnalysisSequence += watchdogOption1vs2

            filterTypeOption1betterThan2 = None
            if integrationMethod == "VEGAS":
                filterTypeOption1betterThan2 = "BestMatchFilterCandidateToSVfitEventHypothesisByIntegration"
            elif integrationMethod == "MarkovChain":
                filterTypeOption1betterThan2 = "BestMatchFilterCandidateToSVfitEventHypothesis"
            else:
                raise ValueError("No BestMatchFilter type defined for option = %i" % option1)
            filterOption1betterThan2 = cms.EDFilter(filterTypeOption1betterThan2,
                srcRef = cms.InputTag(genTauPairs),
                expressionRef = cms.string('mass()'),
                srcTest1 = cms.InputTag(svFitProducerModuleName1),
                srcTest2 = cms.InputTag(svFitProducerModuleName2),
                expressionTest = cms.string('resonance(0).mass()')
            )
            filterNameOption1betterThan2 = "filter%sbetterThan%s" % (svFitProducerModuleName1, svFitProducerModuleName2)
            setattr(process, filterNameOption1betterThan2, filterOption1betterThan2)
        
            analyzerEventOption1betterThan2 = cms.EDAnalyzer("SVfitEventInputAnalyzer",
                srcGenParticles = cms.InputTag('genParticles'),
                srcElectrons = cms.InputTag('genMatchedPatElectrons'),                                            
                srcMuons = cms.InputTag('genMatchedPatMuons'),
                srcTaus = cms.InputTag('genMatchedPatTaus'),
                srcMEt = cms.InputTag('pfType1CorrectedMet'),
                srcMEtCov = cms.InputTag('pfMEtSignCovMatrix'),
                srcGenVertex = cms.InputTag('genEventVertex'),
                srcRecVertex = cms.InputTag('selectedPrimaryVertexByLeptonMatch'),
                srcBeamSpot = cms.InputTag('offlineBeamSpot'),
                algorithm = cms.string("AdaptiveVertexFitter"),
                applyBeamSpotConstraint = cms.bool(False),   
                srcWeights = cms.VInputTag(),
                dqmDirectory = cms.string("%sbetterThan%s/event" % (svFitProducerModuleName1, svFitProducerModuleName2))                                             
            )
            analyzerEventNameOption1betterThan2 = "analyzerEvent%sbetterThan%s" % (svFitProducerModuleName1, svFitProducerModuleName2)
            setattr(process, analyzerEventNameOption1betterThan2, analyzerEventOption1betterThan2)

            analyzerElectronOption1betterThan2 = cms.EDAnalyzer("SVfitTauToElecInputAnalyzer",
                srcGenParticles = cms.InputTag('genParticles'),
                srcElectrons = cms.InputTag('genMatchedPatElectrons'),                                            
                srcMuons = cms.InputTag('genMatchedPatMuons'),
                srcTaus = cms.InputTag('genMatchedPatTaus'),
                srcRecVertex = cms.InputTag('selectedPrimaryVertexByLeptonMatch'),
                srcBeamSpot = cms.InputTag('offlineBeamSpot'),
                algorithm = cms.string("AdaptiveVertexFitter"),
                applyBeamSpotConstraint = cms.bool(False),   
                srcWeights = cms.VInputTag(),
                dqmDirectory = cms.string("%sbetterThan%s/leg1" % (svFitProducerModuleName1, svFitProducerModuleName2))                                             
            )
            analyzerElectronNameOption1betterThan2 = "analyzerElectron%sbetterThan%s" % (svFitProducerModuleName1, svFitProducerModuleName2)
            setattr(process, analyzerElectronNameOption1betterThan2, analyzerElectronOption1betterThan2)

            analyzerMuonOption1betterThan2 = cms.EDAnalyzer("SVfitTauToMuInputAnalyzer",
                srcGenParticles = cms.InputTag('genParticles'),
                srcElectrons = cms.InputTag('genMatchedPatElectrons'),                                            
                srcMuons = cms.InputTag('genMatchedPatMuons'),
                srcTaus = cms.InputTag('genMatchedPatTaus'),
                srcRecVertex = cms.InputTag('selectedPrimaryVertexByLeptonMatch'),
                srcBeamSpot = cms.InputTag('offlineBeamSpot'),
                algorithm = cms.string("AdaptiveVertexFitter"),
                applyBeamSpotConstraint = cms.bool(False),   
                srcWeights = cms.VInputTag(),
                dqmDirectory = cms.string("%sbetterThan%s/leg1" % (svFitProducerModuleName1, svFitProducerModuleName2))                                             
            )
            analyzerMuonNameOption1betterThan2 = "analyzerMuon%sbetterThan%s" % (svFitProducerModuleName1, svFitProducerModuleName2)
            setattr(process, analyzerMuonNameOption1betterThan2, analyzerMuonOption1betterThan2)

            analyzerTauOption1betterThan2 = cms.EDAnalyzer("SVfitTauToHadInputAnalyzer",
                srcGenParticles = cms.InputTag('genParticles'),
                srcElectrons = cms.InputTag('genMatchedPatElectrons'),                                            
                srcMuons = cms.InputTag('genMatchedPatMuons'),
                srcTaus = cms.InputTag('genMatchedPatTaus'),
                srcRecVertex = cms.InputTag('selectedPrimaryVertexByLeptonMatch'),
                srcBeamSpot = cms.InputTag('offlineBeamSpot'),
                algorithm = cms.string("AdaptiveVertexFitter"),
                applyBeamSpotConstraint = cms.bool(False),                                              
                srcWeights = cms.VInputTag(),
                dqmDirectory = cms.string("%sbetterThan%s/leg2" % (svFitProducerModuleName1, svFitProducerModuleName2))                                             
            )
            analyzerTauNameOption1betterThan2 = "analyzerTau%sbetterThan%s" % (svFitProducerModuleName1, svFitProducerModuleName2)
            setattr(process, analyzerTauNameOption1betterThan2, analyzerTauOption1betterThan2)

            pathOption1betterThan2 = cms.Path(
                process.testSVfitTrackLikelihoodProductionSequence
               + filterOption1betterThan2
               + analyzerEventOption1betterThan2 + analyzerMuonOption1betterThan2 + analyzerTauOption1betterThan2
            )
            pathNameOption1betterThan2 = "p%sbetterThan%s" % (svFitProducerModuleName1, svFitProducerModuleName2)
            setattr(process, pathNameOption1betterThan2, pathOption1betterThan2)

            filterOption2betterThan1 = filterOption1betterThan2.clone(
                srcTest1 = cms.InputTag(svFitProducerModuleName2),
                srcTest2 = cms.InputTag(svFitProducerModuleName1)
            )
            filterNameOption2betterThan1 = "filter%sbetterThan%s" % (svFitProducerModuleName2, svFitProducerModuleName1)
            setattr(process, filterNameOption2betterThan1, filterOption2betterThan1)
        
            analyzerEventOption2betterThan1 = analyzerEventOption1betterThan2.clone(
                dqmDirectory = cms.string("%sbetterThan%s/event" % (svFitProducerModuleName2, svFitProducerModuleName1))                                             
            )
            analyzerEventNameOption2betterThan1 = "analyzerEvent%sbetterThan%s" % (svFitProducerModuleName2, svFitProducerModuleName1)
            setattr(process, analyzerEventNameOption2betterThan1, analyzerEventOption2betterThan1)

            analyzerElectronOption2betterThan1 = analyzerElectronOption1betterThan2.clone(
                dqmDirectory = cms.string("%sbetterThan%s/leg1" % (svFitProducerModuleName2, svFitProducerModuleName1))                                             
            )
            analyzerElectronNameOption2betterThan1 = "analyzerElectron%sbetterThan%s" % (svFitProducerModuleName2, svFitProducerModuleName1)
            setattr(process, analyzerElectronNameOption2betterThan1, analyzerElectronOption2betterThan1)

            analyzerMuonOption2betterThan1 = analyzerMuonOption1betterThan2.clone(
                dqmDirectory = cms.string("%sbetterThan%s/leg1" % (svFitProducerModuleName2, svFitProducerModuleName1))                                             
            )
            analyzerMuonNameOption2betterThan1 = "analyzerMuon%sbetterThan%s" % (svFitProducerModuleName2, svFitProducerModuleName1)
            setattr(process, analyzerMuonNameOption2betterThan1, analyzerMuonOption2betterThan1)

            analyzerTauOption2betterThan1 = analyzerTauOption1betterThan2.clone(
                dqmDirectory = cms.string("%sbetterThan%s/leg2" % (svFitProducerModuleName2, svFitProducerModuleName1))                                             
            )
            analyzerTauNameOption2betterThan1 = "analyzerTau%sbetterThan%s" % (svFitProducerModuleName2, svFitProducerModuleName1)
            setattr(process, analyzerTauNameOption2betterThan1, analyzerTauOption2betterThan1)

            pathOption2betterThan1 = cms.Path(
                process.testSVfitTrackLikelihoodProductionSequence
               + filterOption2betterThan1
               + analyzerEventOption2betterThan1 + analyzerMuonOption2betterThan1 + analyzerTauOption2betterThan1
            )
            pathNameOption2betterThan1 = "p%sbetterThan%s" % (svFitProducerModuleName2, svFitProducerModuleName1)
            setattr(process, pathNameOption2betterThan1, pathOption2betterThan1)
#--------------------------------------------------------------------------------

process.testSVfitTrackLikelihoodSequence = cms.Sequence(process.testSVfitTrackLikelihoodProductionSequence + process.testSVfitTrackLikelihoodAnalysisSequence)

process.p = cms.Path(process.testSVfitTrackLikelihoodSequence)

#--------------------------------------------------------------------------------
# save plots

process.DQMStore = cms.Service("DQMStore")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string("/data1/veelken/tmp/svFitStudies/testSVfitTrackLikelihoods4_plots_%s_%s_%i_2013Aug08.root" % (sample_type, channel, massPoint))
)

process.q = cms.EndPath(process.savePlots)
#--------------------------------------------------------------------------------

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

processDumpFile = open('testSVfitTrackLikelihoods4.dump' , 'w')
print >> processDumpFile, process.dumpPython()

