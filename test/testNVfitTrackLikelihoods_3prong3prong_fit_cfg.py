import FWCore.ParameterSet.Config as cms

process = cms.Process("testSVfitTrackLikelihoods")

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
    input = cms.untracked.int32(-1)
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
    ##    '1:3369:1346645',
    ##    '1:89621:35819785'
    ##    '1:23167:9259467'
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
#massPoint = 125
massPoint = 300
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
    if channel == 'tautau':
        inputFilePaths = [
            '/data2/veelken/CMSSW_5_3_x/skims/user/veelken/CMSSW_5_3_x/skims/SVfitStudies/AccLowPtThresholds/ZplusJets_tautau_highStats/'
        ]
    elif channel == 'mutau':
        inputFilePaths = [
            '/data2/veelken/CMSSW_5_3_x/skims/user/veelken/CMSSW_5_3_x/skims/SVfitStudies/AccLowPtThresholds/ZplusJets_mutau/'
        ]
    elif channel == 'emu':
        inputFilePaths = [
            '/data2/veelken/CMSSW_5_3_x/skims/user/veelken/CMSSW_5_3_x/skims/SVfitStudies/AccLowPtThresholds/ZplusJets_emu/'
        ]
    else:
        raise ValueError("Invalid channel = %s !!" % channel)    
    inputFile_regex = \
      r"[a-zA-Z0-9_/:.]*genTauLeptonPairSkim_ZplusJets_%s_(?P<gridJob>\d*)(_(?P<gridTry>\d*))*_(?P<hash>[a-zA-Z0-9]*).root" % channel
elif sample_type == 'Higgs':
    if channel == 'tautau':
        if massPoint == 125:
            inputFilePaths = [
                '/data2/veelken/CMSSW_5_3_x/skims/user/veelken/CMSSW_5_3_x/skims/SVfitStudies/AccLowPtThresholds/ggHiggs125_tautau_highStats/',
                '/data2/veelken/CMSSW_5_3_x/skims/user/veelken/CMSSW_5_3_x/skims/SVfitStudies/AccLowPtThresholds/vbfHiggs125_tautau_highStats/'
            ]
        elif massPoint == 300:
            inputFilePaths = [
                '/data2/veelken/CMSSW_5_3_x/skims/user/veelken/CMSSW_5_3_x/skims/SVfitStudies/AccLowPtThresholds/ggPhi300_tautau_highStats/'
            ]
        else:
            raise ValueError("Invalid mass-point = %i !!" % massPoint)
    elif channel == 'mutau':
        if massPoint == 125:
            inputFilePaths = [
                '/data2/veelken/CMSSW_5_3_x/skims/user/veelken/CMSSW_5_3_x/skims/SVfitStudies/AccLowPtThresholds/ggHiggs125_mutau/',
                '/data2/veelken/CMSSW_5_3_x/skims/user/veelken/CMSSW_5_3_x/skims/SVfitStudies/AccLowPtThresholds/vbfHiggs125_mutau/'
            ]
        elif massPoint == 300:
            inputFilePaths = [
                '/data2/veelken/CMSSW_5_3_x/skims/user/veelken/CMSSW_5_3_x/skims/SVfitStudies/AccLowPtThresholds/ggPhi300_mutau/'
            ]
        else:
            raise ValueError("Invalid mass-point = %i !!" % massPoint)
    elif channel == 'emu':   
        if massPoint == 125:
            inputFilePaths = [
                '/data2/veelken/CMSSW_5_3_x/skims/user/veelken/CMSSW_5_3_x/skims/SVfitStudies/AccLowPtThresholds/ggHiggs125_emu/',
                '/data2/veelken/CMSSW_5_3_x/skims/user/veelken/CMSSW_5_3_x/skims/SVfitStudies/AccLowPtThresholds/vbfHiggs125_emu/'
            ]
        elif massPoint == 300:
            inputFilePaths = [
                '/data2/veelken/CMSSW_5_3_x/skims/user/veelken/CMSSW_5_3_x/skims/SVfitStudies/AccLowPtThresholds/ggPhi300_emu/'
            ]
        else:
            raise ValueError("Invalid mass-point = %i !!" % massPoint)
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
    #'oneProng0Pi0', 
    #'oneProng1Pi0', 
    #'oneProng2Pi0', 
    #'oneProngOther',
    'threeProng0Pi0', 
    #'threeProng1Pi0', 
    #'threeProngOther', 
    #'rare'
)
if channel == 'mutau' or channel == 'etau' or channel == 'tautau':
    process.tauGenJetsSelectorAllHadrons.filter = cms.bool(True)
elif channel == 'emu':
    process.tauGenJetsSelectorAllHadrons.filter = cms.bool(False)
else:
    raise ValueError("Invalid channel = %s !!" % channel)    
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

process.load("RecoTauTag/ImpactParameter/PFTau3ProngReco_cfi")
process.PFTau3ProngReco.PFTauTag = cms.InputTag('hpsPFTauProducer')
process.PFTau3ProngReco.PFTauTIPTag = cms.InputTag('hpsPFTauTransverseImpactParameters')
process.PFTau3ProngReco.Algorithm = cms.int32(0)
process.PFTau3ProngReco.discriminators = cms.VPSet()
process.PFTau3ProngReco.cut = cms.string("")
process.testSVfitTrackLikelihoodProductionSequence += process.PFTau3ProngReco

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
process.svFitTauToElecBuilder.initializeToGen = cms.bool(True)
process.svFitTauToElecBuilder.srcGenTaus = cms.InputTag('genParticles')
process.svFitTauToElecBuilder.dRmatch = cms.double(0.3)
process.svFitTauToElecBuilder.verbosity = cms.int32(0)

process.svFitTauToMuBuilder.fixToGenVisEnFracX = cms.bool(False)
process.svFitTauToMuBuilder.fixToGenPhiLab = cms.bool(False)
process.svFitTauToMuBuilder.fixToGenNuInvMass = cms.bool(False)
process.svFitTauToMuBuilder.fixRecToGenVertex = cms.bool(False)
process.svFitTauToMuBuilder.fixToGenDecayDistance = cms.bool(False)
process.svFitTauToMuBuilder.fixToGenVisP4 = cms.bool(False)
process.svFitTauToMuBuilder.initializeToGen = cms.bool(True)
process.svFitTauToMuBuilder.srcGenTaus = cms.InputTag('genParticles')
process.svFitTauToMuBuilder.dRmatch = cms.double(0.3)
process.svFitTauToMuBuilder.verbosity = cms.int32(0)

process.svFitTauToHadBuilder.fixToGenVisEnFracX = cms.bool(False)
process.svFitTauToHadBuilder.fixToGenPhiLab = cms.bool(False)
process.svFitTauToHadBuilder.fixRecToGenVertex = cms.bool(False)
process.svFitTauToHadBuilder.fixToGenDecayDistance = cms.bool(False)
process.svFitTauToHadBuilder.fixToGenVisP4 = cms.bool(False)
process.svFitTauToHadBuilder.initializeToGen = cms.bool(True)
process.svFitTauToHadBuilder.srcGenTaus = cms.InputTag('genParticles')
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

likelihoodOptions = [
    '',
    'trackWOlifetimeWdecayVertexForLeg1Leg2',
    'trackWOlifetimeWdecayVertexForLeg1Leg2PlusDecayKine',
    'trackWOlifetimeWdecayVertexForLeg1Leg2PlusMEt',
    'trackWOlifetimeWdecayVertexForLeg1Leg2PlusDecayKinePlusMEt',
    'trackWOlifetimeWdecayVertexForLeg1PlusDecayKinePlusMEt',
    'trackWOlifetimeWdecayVertexForLeg2PlusDecayKinePlusMEt',
    'trackWOlifetimePlusDecayKinePlusMEt',
    'trackWlifetimeWdecayVertexForLeg1Leg2',
    'trackWlifetimeWdecayVertexForLeg1Leg2PlusDecayKine',
    'trackWlifetimeWdecayVertexForLeg1Leg2PlusMEt',
    'trackWlifetimeWdecayVertexForLeg1Leg2PlusDecayKinePlusMEt',
    'trackWlifetimeWdecayVertexForLeg1PlusDecayKinePlusMEt',
    'trackWlifetimeWdecayVertexForLeg2PlusDecayKinePlusMEt',
    'trackWlifetimePlusDecayKinePlusMEt',
    'decayKine',
    'decayKinePlusMEt',
    'met'
]

def customCapitalize(name):
    retVal = name[0:1].capitalize() + name[1:]
    return retVal

for likelihoodOption in likelihoodOptions:
    for fixRecToGenVertex in [ False ]:
        for initializeToGen in [ True, False ]:

            fixRecToGenVertex_string = None
            if fixRecToGenVertex:
                fixRecToGenVertex_string = "genVertex"
            else:
                fixRecToGenVertex_string = "recVertex"

            initializeToGen_string = None
            if initializeToGen:
                initializeToGen_string = "WgenInitialization"
            else:
                initializeToGen_string = "WOgenInitialization"

            print "processing likelihood option = %s, vertex = %s, initialization = %s" % (likelihoodOption, fixRecToGenVertex_string, initializeToGen_string)

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
                eventBuilder.fixToGenVertex = cms.bool(True)
                eventBuilder.srcGenVertex = cms.InputTag('genEventVertex')
            else:
                leg1Builder.fixRecToGenVertex = cms.bool(False)
                leg2Builder.fixRecToGenVertex = cms.bool(False)
                eventBuilder.fixToGenVertex = cms.bool(False)
            if initializeToGen:
                leg1Builder.initializeToGen = cms.bool(True)
                leg2Builder.initializeToGen = cms.bool(True)
            else:
                leg1Builder.initializeToGen = cms.bool(False)
                leg2Builder.initializeToGen = cms.bool(False)
            resonanceLikelihoodLogM = copy.deepcopy(process.svFitResonanceLikelihoodLogM)
            resonanceLikelihoodFunctions.append(resonanceLikelihoodLogM)

            # MINUIT maximization of likelihood functions
            svFitProducerModuleMINUIT = process.svFitProducerByLikelihoodMaximization.clone()
            svFitProducerModuleMINUIT.config.event.resonances.A.daughters.leg1.src = cms.InputTag(srcRecLeg1)
            svFitProducerModuleMINUIT.config.event.resonances.A.daughters.leg1.builder = leg1Builder
            svFitProducerModuleMINUIT.config.event.resonances.A.daughters.leg1.likelihoodFunctions = cms.VPSet(leg1LikelihoodFunctions)
            svFitProducerModuleMINUIT.config.event.resonances.A.daughters.leg2.src = cms.InputTag(srcRecLeg2)
            svFitProducerModuleMINUIT.config.event.resonances.A.daughters.leg2.builder = leg2Builder
            svFitProducerModuleMINUIT.config.event.resonances.A.daughters.leg2.likelihoodFunctions = cms.VPSet(leg2LikelihoodFunctions)
            svFitProducerModuleMINUIT.config.event.resonances.A.likelihoodFunctions = cms.VPSet(resonanceLikelihoodFunctions)
            svFitProducerModuleMINUIT.config.event.builder = eventBuilder
            svFitProducerModuleMINUIT.config.event.likelihoodFunctions = cms.VPSet(eventLikelihoodFunctions)
            svFitProducerModuleMINUIT.config.event.srcMEt = cms.InputTag('pfType1CorrectedMet')
            svFitProducerModuleMINUIT.config.event.srcPrimaryVertex = cms.InputTag('selectedPrimaryVertexByLeptonMatch')
            pluginNameMINUIT = "svFitProducerByLikelihoodMaximizationW%s%s%s" % \
              (likelihoodOption, initializeToGen_string, customCapitalize(fixRecToGenVertex_string))
            svFitProducerModuleMINUIT.algorithm.pluginName = cms.string(pluginNameMINUIT)
            svFitProducerModuleMINUIT.algorithm.verbosity = cms.int32(1)
            svFitProducerModuleNameMINUIT = "svFitProducerByLikelihoodMaximizationW%s%s%s" % \
              (likelihoodOption, initializeToGen_string, customCapitalize(fixRecToGenVertex_string))
            svFitAnalyzerModuleTypeMINUIT = "SVfitEventHypothesisAnalyzer"
            print "adding module '%s' of type MINUIT" % svFitProducerModuleNameMINUIT

            keyMINUIT = "%s%s%sMINUIT" % (likelihoodOption, initializeToGen_string, customCapitalize(fixRecToGenVertex_string))
            svFitProducerModuleNames[keyMINUIT] = svFitProducerModuleNameMINUIT
            svFitAnalyzerModuleTypes[keyMINUIT] = svFitAnalyzerModuleTypeMINUIT
            
            setattr(process, svFitProducerModuleNameMINUIT, svFitProducerModuleMINUIT)
            process.testSVfitTrackLikelihoodProductionSequence += svFitProducerModuleMINUIT

##print "svFitProducerModuleNames:"
##print svFitProducerModuleNames            
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
    svFitAnalyzerModuleName = svFitProducerModuleName.replace("Producer", "Analyzer")
    setattr(process, svFitAnalyzerModuleName, svFitAnalyzerModule)
    process.testSVfitTrackLikelihoodAnalysisSequence += svFitAnalyzerModule
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# fill histograms of mass reconstructed by Ian's 3-prong kinematic reconstruction tool

process.threeProngKineReconstructionAnalyzerModule = cms.EDAnalyzer("PATTau3ProngRecoAnalyzer",    
    srcPATTaus = cms.InputTag('genMatchedPatTaus'),
    srcPFTaus = cms.InputTag('hpsPFTauProducer'),
    srcPFTau3ProngReco = cms.InputTag('PFTau3ProngReco'),
    srcGenTauPairs = cms.InputTag(genTauPairs),
    srcPFMEt = cms.InputTag('pfType1CorrectedMet'),  
    srcPFMEtCovMatrix = cms.InputTag('pfMEtSignCovMatrix'),
    srcGenMEt = cms.InputTag('genMetFromGenParticles'),
    srcWeights = cms.VInputTag(),                                                    
    dqmDirectory = cms.string('threeProngKineReconstructionAnalyzerModule'),
    verbosity = cms.int32(1)
)                                     
process.testSVfitTrackLikelihoodAnalysisSequence += process.threeProngKineReconstructionAnalyzerModule

process.dumpGenTaus = cms.EDAnalyzer("DumpGenTaus",
    src = cms.InputTag('genParticles')
)
##process.testSVfitTrackLikelihoodAnalysisSequence += process.dumpGenTaus

process.dumpPFTaus = cms.EDAnalyzer("DumpPFTaus",
    src = cms.InputTag('hpsPFTauProducer'),
    srcTauRef = cms.InputTag('hpsPFTauProducer'),
    srcTauRefDiscriminators = cms.VInputTag(),
    srcVertex = cms.InputTag('selectedPrimaryVertexByLeptonMatch')                                
)
##process.testSVfitTrackLikelihoodAnalysisSequence += process.dumpPFTaus

process.dumpPFTau3ProngReco = cms.EDAnalyzer("DumpPFTau3ProngReco",
    srcTau3ProngReco = cms.InputTag('PFTau3ProngReco'),
    srcTau = cms.InputTag('hpsPFTauProducer')
)
process.testSVfitTrackLikelihoodAnalysisSequence += process.dumpPFTau3ProngReco
#--------------------------------------------------------------------------------

process.testSVfitTrackLikelihoodSequence = cms.Sequence(process.testSVfitTrackLikelihoodProductionSequence + process.testSVfitTrackLikelihoodAnalysisSequence)

process.p = cms.Path(process.testSVfitTrackLikelihoodSequence)

#--------------------------------------------------------------------------------
# save plots

process.DQMStore = cms.Service("DQMStore")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string("/data1/veelken/tmp/svFitStudies/testSVfitTrackLikelihoods_3prong3prong_fit_plots_%s_%s_%i_2013Aug14.root" % (sample_type, channel, massPoint))
)

process.q = cms.EndPath(process.savePlots)
#--------------------------------------------------------------------------------

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

processDumpFile = open('testSVfitTrackLikelihoods_3prong3prong.dump' , 'w')
print >> processDumpFile, process.dumpPython()

