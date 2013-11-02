import FWCore.ParameterSet.Config as cms

process = cms.Process("testSVfitTrackLikelihoods2")

import os
import re

import TauAnalysis.Configuration.tools.castor as castor
from TauAnalysis.Skimming.EventContent_cff import *

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load('Configuration/Geometry/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START53_V15::All')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/data1/veelken/CMSSW_5_2_x/skims/goldenZmumuEvents_ZplusJets_madgraph2_2012Apr12_AOD_9_1_cSC.root'
    ),
    dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
    inputCommands=cms.untracked.vstring(
        'keep *',
        'drop LHERunInfoProduct_*_*_*',                        
        'drop LHEEventProduct_*_*_*'
    ),
    ##eventsToProcess = cms.untracked.VEventRange(
    ##    '1:150755:60253795'
    ##)
)

#--------------------------------------------------------------------------------
# define configuration parameter default values

#sample_type = 'Z'
sample_type = 'Higgs'
#channel = 'etau'
channel = 'mutau'
#channel = 'emu'
massPoint = '125'
#massPoint = '300'
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define "hooks" for replacing configuration parameters
# in case running jobs on the CERN batch system/grid
#
#__sample_type = '#sample_type#'
#__channel = '#channel#'
#__massPoint = '#massPoint#'
#
if sample_type == 'Z':
    massPoint = '90'
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# set input files
inputFilePath = None
inputFile_regex = None
if sample_type == 'Z':
    inputFilePath = '/data1/veelken/CMSSW_5_2_x/skims/genHtautauLeptonPairAcc/user/veelken/CMSSW_5_2_x/skims/'
    inputFile_regex = \
      r"[a-zA-Z0-9_/:.]*genTauLeptonsPairAccSkim_ZplusJets_%s_(?P<gridJob>\d*)(_(?P<gridTry>\d*))*_(?P<hash>[a-zA-Z0-9]*).root" % channel
elif sample_type == 'Higgs':
    inputFilePath = '/data1/veelken/CMSSW_5_2_x/skims/genHtautauLeptonPairAcc/user/v/veelken/CMSSW_5_2_x/skims/'
    inputFile_regex = \
      r"[a-zA-Z0-9_/:.]*genTauLeptonsPairAccSkim_(ggHiggs|ggPhi|vbfHiggs)%s_%s_(?P<gridJob>\d*)(_(?P<gridTry>\d*))*_(?P<hash>[a-zA-Z0-9]*).root" % (massPoint, channel)
else:
    raise ValueError("Invalid sample type = %s !!" % sample_type)

# check if name of inputFile matches regular expression
inputFileNames = []
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
#print "inputFileNames = %s" % inputFileNames 

process.source.fileNames = cms.untracked.vstring(inputFileNames)
#--------------------------------------------------------------------------------

process.testSVfitTrackLikelihoodProductionSequence = cms.Sequence()

#--------------------------------------------------------------------------------
# select collections of electrons, muons and tau-jets
# matching genuine tau -> e, tau -> mu and tau -> hadronic decays on generator level

process.load("RecoTauTag/Configuration/RecoPFTauTag_cff")
process.testSVfitTrackLikelihoodProductionSequence += process.recoTauClassicHPSSequence

process.load("PhysicsTools/JetMCAlgos/TauGenJets_cfi")
process.testSVfitTrackLikelihoodProductionSequence += process.tauGenJets

genElectronsFromTauDecays = None
genMuonsFromTauDecays = None
genTauJetsFromTauDecays = None
genTauPairs = None
genTaus = None
if sample_type == 'Z':
    process.load("TauAnalysis/GenSimTools/gen_decaysFromZs_cfi")
    process.testSVfitTrackLikelihoodProductionSequence += process.produceGenDecayProductsFromZs
    genElectronsFromTauDecays = 'genElectronsFromZtautauDecays'
    genMuonsFromTauDecays = 'genMuonsFromZtautauDecays'
    genTauJetsFromTauDecays = 'genHadronsFromZtautauDecays'
    genTauPairs = 'genZdecayToTaus'
    genTaus = 'genTausFromZs'
elif sample_type == 'Higgs':
    process.load("TauAnalysis/GenSimTools/gen_decaysFromAHs_cfi")
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
    cut = cms.string('mass > %1.1f & mass < %1.1f' % (0.70*massPoint, 1.30*massPoint))
    filter = cms.bool(True)
)
process.testSVfitTrackLikelihoodProductionSequence += process.process.genTauPairMassFilter

process.load("TauAnalysis/Skimming/goldenZmmSelectionVBTFnoMuonIsolation_cfi")
process.goodMuons. cut = cms.string(
    'pt > 1. & abs(eta) < 2.5 & isGlobalMuon' \
    + ' & innerTrack.hitPattern.numberOfValidTrackerHits > 9 & innerTrack.hitPattern.numberOfValidPixelHits > 0' \
    + ' & abs(dB) < 0.2 & globalTrack.normalizedChi2 < 10' \
    + ' & globalTrack.hitPattern.numberOfValidMuonHits > 0 & numberOfMatches > 1' 
)
process.goodIsoMuons.cut = cms.string(
    '(userIsolation("pat::User1Iso")' + \
    ' + max(0., userIsolation("pat::PfNeutralHadronIso") + userIsolation("pat::PfGammaIso")' + \
    '          - 0.5*userIsolation("pat::User2Iso"))) < 0.10*pt'
)
process.muonSelectionSequence = cms.Sequence(
    process.pfNoPileUpSequence
   + process.pfParticleSelectionSequence
   + process.muonPFIsolationDepositsSequence
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

process.selectedTaus = cms.EDFilter("PFTauSelector",
    src = cms.InputTag('hpsPFTauProducer'),
    discriminators = cms.VPSet(
        cms.PSet(
            discriminator = cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding'),
            selectionCut = cms.double(0.5)
        ),
        cms.PSet(
            discriminator = cms.InputTag('hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr'),
            selectionCut = cms.double(0.5)
        )                        
    ),
    cut = cms.string("pt > 20. & abs(eta) < 2.3")                        
)
process.testSVfitTrackLikelihoodProductionSequence += process.selectedTaus

from PhysicsTools.PatAlgos.producersLayer1.tauProducer_cfi import patTaus
process.patTausForSVfit = patTaus.clone(
    tauSource = cms.InputTag("selectedTaus"),
    isoDeposits = cms.PSet(),
    userIsolation = cms.PSet(),
    addTauID = cms.bool(False),
    tauIDSources = cms.PSet(),
    addGenMatch = cms.bool(False),
    embedGenMatch = cms.bool(False),
    addGenJetMatch = cms.bool(False),
    embedGenJetMatch = cms.bool(False)
)
process.testSVfitTrackLikelihoodProductionSequence += process.patTausForSVfit

process.genMatchedPatTaus = cms.EDFilter("PATTauAntiOverlapSelector",
    src = cms.InputTag('patTausForSVfit'),
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
elif channel == 'ditau':
    numElectrons = 0
    numMuons     = 0
    numTauJets   = 2  
else:
    raise ValueError("Invalid channel = %s !!" % channel)

if numElectrons > 0:
    process.electronFilter = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag('genMatchedPatElectrons'),
        minNumber = cms.uint32(numElectrons),
        maxNumber = cms.uint32(numElectrons)                      
    ) 
    process.testSVfitTrackLikelihoodProductionSequence += process.electronFilter

if numMuons > 0:    
    process.muonFilter = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag('genMatchedPatMuons'),
        minNumber = cms.uint32(numMuons),
        maxNumber = cms.uint32(numMuons)                      
    )
    process.testSVfitTrackLikelihoodProductionSequence += process.muonFilter

if numTauJets > 0:
    process.tauFilter = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag('genMatchedPatTaus'),
        minNumber = cms.uint32(numTauJets),
        maxNumber = cms.uint32(numTauJets)                      
    )
    process.testSVfitTrackLikelihoodProductionSequence += process.tauFilter
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
# reconstructed Type 1 corrected PFMET and its uncertainty 

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
        DPhi = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmetphi_53.root'), # CV: same for unity and non-unity response training
        CovU2 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru2cov_53_UnityResponse.root'),
        U = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmet_53_UnityResponse.root'),
        CovU1 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbru1cov_53_UnityResponse.root')
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
# run SVfit with and without Track likelihoods

process.load("TauAnalysis.CandidateTools.svFitAlgorithmDiTau_cfi")

# CV: fix tau decay parameters to Monte Carlo truth values
process.svFitTauToElecBuilder.fixToGenVisEnFracX = cms.bool(False)
process.svFitTauToElecBuilder.fixToGenPhiLab = cms.bool(False)
process.svFitTauToElecBuilder.fixToGenNuInvMass = cms.bool(False)
process.svFitTauToElecBuilder.fixToGenDeltaR = cms.bool(False)
process.svFitTauToElecBuilder.fixToGenVisP4 = cms.bool(False)
process.svFitTauToElecBuilder.initializeToGen = cms.bool(False)
#process.svFitTauToElecBuilder.srcGenTaus = cms.InputTag('genParticles')
process.svFitTauToElecBuilder.dRmatch = cms.double(0.3)
process.svFitTauToElecBuilder.verbosity = cms.int32(0)

process.svFitTauToMuBuilder.fixToGenVisEnFracX = cms.bool(False)
process.svFitTauToMuBuilder.fixToGenPhiLab = cms.bool(False)
process.svFitTauToMuBuilder.fixToGenNuInvMass = cms.bool(False)
process.svFitTauToMuBuilder.fixToGenDeltaR = cms.bool(False)
process.svFitTauToMuBuilder.fixToGenVisP4 = cms.bool(False)
process.svFitTauToMuBuilder.initializeToGen = cms.bool(False)
#process.svFitTauToMuBuilder.srcGenTaus = cms.InputTag('genParticles')
process.svFitTauToMuBuilder.dRmatch = cms.double(0.3)
process.svFitTauToMuBuilder.verbosity = cms.int32(0)

process.svFitTauToHadBuilder.fixToGenVisEnFracX = cms.bool(False)
process.svFitTauToHadBuilder.fixToGenPhiLab = cms.bool(False)
process.svFitTauToHadBuilder.fixToGenDeltaR = cms.bool(False)
process.svFitTauToHadBuilder.fixToGenVisP4 = cms.bool(False)
process.svFitTauToHadBuilder.initializeToGen = cms.bool(False)
#process.svFitTauToHadBuilder.srcGenTaus = cms.InputTag('genParticles')
process.svFitTauToHadBuilder.dRmatch = cms.double(0.3)
process.svFitTauToHadBuilder.verbosity = cms.int32(0)

# CV: fix event vertex position to Monte Carlo truth value
process.svFitEventBuilder.fixToGenVertex = cms.bool(False)
#process.svFitEventBuilder.srcGenVertex = cms.InputTag('genEventVertex')

srcGenLeg1                      = None
srcRecLeg1                      = None
svFitLikelihoodLeg1_kinematics = None
svFitBuilderLeg1               = None
svFitStandAloneTypeLeg1        = None
srcGenLeg2                      = None
srcRecLeg2                      = None
svFitLikelihoodLeg2_kinematics = None
svFitBuilderLeg2               = None
svFitStandAloneTypeLeg2        = None
if channel == 'mutau':
    srcGenLeg1                      = genMuonsFromTauDecays
    srcRecLeg1                      = 'genMatchedPatMuons'
    svFitLikelihoodLeg1_kinematics = process.svFitMuonLikelihoodMatrixElement
    svFitBuilderLeg1               = process.svFitTauToMuBuilder
    svFitStandAloneTypeLeg1        = "lep"
    srcGenLeg2                      = genTauJetsFromTauDecays
    srcRecLeg2                      = 'genMatchedPatTaus'
    svFitLikelihoodLeg2_kinematics = process.svFitTauLikelihoodPhaseSpace
    svFitBuilderLeg2               = process.svFitTauToHadBuilder
    svFitStandAloneTypeLeg2        = "had"
elif channel == 'etau':
    srcGenLeg1                      = genElectronsFromTauDecays
    srcRecLeg1                      = 'genMatchedPatElectrons'
    svFitLikelihoodLeg1_kinematics = process.svFitElectronLikelihoodMatrixElement
    svFitBuilderLeg1               = process.svFitTauToElecBuilder
    svFitStandAloneTypeLeg1        = "lep"
    srcGenLeg2                      = genTauJetsFromTauDecays
    srcRecLeg2                      = 'genMatchedPatTaus'
    svFitLikelihoodLeg2_kinematics = process.svFitTauLikelihoodPhaseSpace
    svFitBuilderLeg2               = process.svFitTauToHadBuilder
    svFitStandAloneTypeLeg2        = "had"
elif channel == 'emu':
    srcGenLeg1                      = genElectronsFromTauDecays
    srcRecLeg1                      = 'genMatchedPatElectrons'
    svFitLikelihoodLeg1_kinematics = process.svFitElectronLikelihoodMatrixElement
    svFitBuilderLeg1               = process.svFitTauToElecBuilder
    svFitStandAloneTypeLeg1        = "lep"
    srcGenLeg2                      = genMuonsFromTauDecays
    srcRecLeg2                      = 'genMatchedPatMuons'
    svFitLikelihoodLeg2_kinematics = process.svFitMuonLikelihoodMatrixElement
    svFitBuilderLeg2               = process.svFitTauToMuBuilder
    svFitStandAloneTypeLeg2        = "lep"
elif channel == 'ditau':
    srcGenLeg1                      = genTauJetsFromTauDecays
    srcRecLeg1                      = 'genMatchedPatTaus'
    svFitLikelihoodLeg1_kinematics = process.svFitTauLikelihoodPhaseSpace
    svFitBuilderLeg1               = process.svFitTauToHadBuilder
    svFitStandAloneTypeLeg1        = "had"
    srcGenLeg2                      = genTauJetsFromTauDecays
    srcRecLeg2                      = 'genMatchedPatTaus'
    svFitLikelihoodLeg2_kinematics = process.svFitTauLikelihoodPhaseSpace
    svFitBuilderLeg2               = process.svFitTauToHadBuilder
    svFitStandAloneTypeLeg2        = "had"
else:
    raise ValueError("Invalid channel = %s !!" % channel)

svFitProducerModuleNames = dict()
svFitAnalyzerModuleTypes = dict()

# option 1: VEGAS integration of likelihood functions, no tracking information used
svFitProducerModule = process.svFitProducerByIntegration.clone()
svFitProducerModule.config.event.resonances.A.daughters.leg1.src = cms.InputTag(srcRecLeg1)
svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions = cms.VPSet(svFitLikelihoodLeg1_kinematics)
svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions[0].verbosity = cms.int32(0)
svFitProducerModule.config.event.resonances.A.daughters.leg1.builder = svFitBuilderLeg1
svFitProducerModule.config.event.resonances.A.daughters.leg2.src = cms.InputTag(srcRecLeg2)
svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions = cms.VPSet(svFitLikelihoodLeg2_kinematics)
svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[0].verbosity = cms.int32(0)
svFitProducerModule.config.event.resonances.A.daughters.leg2.builder = svFitBuilderLeg2
svFitProducerModule.config.event.resonances.A.likelihoodFunctions = cms.VPSet()
svFitProducerModule.config.event.srcMEt = cms.InputTag('pfType1CorrectedMet')
svFitProducerModule.config.event.srcPrimaryVertex = cms.InputTag('selectedPrimaryVertexByLeptonMatch')
svFitProducerModule.algorithm.pluginName = cms.string("svFitProducerByIntegrationWOtracks")
svFitProducerModule.algorithm.verbosity = cms.int32(0)
svFitProducerModuleName = "svFitProducerByIntegrationWOtracks"
svFitAnalyzerModuleType = "SVfitEventHypothesisByIntegrationAnalyzer"

svFitProducerModuleNames["VEGAS"] = svFitProducerModuleName
svFitAnalyzerModuleTypes["VEGAS"] = svFitAnalyzerModuleType

setattr(process, svFitProducerModuleName, svFitProducerModule)
process.testSVfitTrackLikelihoodProductionSequence += svFitProducerModule

metAlgorithms = {
    'PFMEt' : {
        'srcMEt'        : 'pfType1CorrectedMet',
        'srcMEtCov'     : 'pfMEtSignCovMatrix',
        'metLikelihood' : "svFitEventLikelihoodMEt2"
    },
##    'MVAMEtUnityResponse' : {
##        'srcMEt'        : 'pfMEtMVAunityResponse',
##        'srcMEtCov'     : 'pfMEtMVAunityResponse',
##        'metLikelihood' : "svFitEventLikelihoodMEt2"
##     },
##     'MVAMEtUnityResponsePFMEtCov' : {
##         'srcMEt'       : 'pfMEtMVAunityResponseWithPFMEtCov',
##         'srcMEtCov'    : 'pfMEtMVAunityResponseWithPFMEtCov',
##        'metLikelihood' : "svFitEventLikelihoodMEt2"
##     },
    'MVAMEtNonUnityResponse2' : {
        'srcMEt'        : 'pfMEtMVA',
        'srcMEtCov'     : 'pfMEtMVA',
        'metLikelihood' : "svFitEventLikelihoodMEt2"
        
    },
    'MVAMEtNonUnityResponse2pfMEtCov' : {
        'srcMEt'        : 'pfMEtMVAwithPFMEtCov',
        'srcMEtCov'     : 'pfMEtMVAwithPFMEtCov',
        'metLikelihood' : "svFitEventLikelihoodMEt2"
    },
    'MVAMEtNonUnityResponse3' : {
        'srcMEt'        : 'pfMEtMVA',
        'srcMEtCov'     : 'pfMEtMVA',
        'metLikelihood' : "svFitEventLikelihoodMEt3"
    }
}

##for numCalls in [ 10, 20, 50, 100, 250, 500, 1000 ]:
##for numCalls in [ 5000 ]:
for numCalls in [ 100 ]:
    for metAlgorithmName in metAlgorithms.keys():
        
        # Markov Chain integration of likelihood functions, no tracking information used
        svFitProducerModule = process.svFitProducerByIntegration2.clone()
        svFitProducerModule.config.event.resonances.A.daughters.leg1.src = cms.InputTag(srcRecLeg1)
        svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions = cms.VPSet(svFitLikelihoodLeg1_kinematics)
        svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
        svFitProducerModule.config.event.resonances.A.daughters.leg1.builder = svFitBuilderLeg1
        svFitProducerModule.config.event.resonances.A.daughters.leg2.src = cms.InputTag(srcRecLeg2)
        svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions = cms.VPSet(svFitLikelihoodLeg2_kinematics)
        svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
        svFitProducerModule.config.event.resonances.A.daughters.leg2.builder = svFitBuilderLeg2
        svFitProducerModule.config.event.resonances.A.likelihoodFunctions = cms.VPSet()
        svFitProducerModule.config.event.srcMEt = cms.InputTag(metAlgorithms[metAlgorithmName]['srcMEt'])
        svFitProducerModule.config.event.likelihoodFunctions[0] = getattr(process, metAlgorithms[metAlgorithmName]['metLikelihood']).clone()
        svFitProducerModule.config.event.likelihoodFunctions[0].srcMEtCovMatrix = cms.InputTag(metAlgorithms[metAlgorithmName]['srcMEtCov'])
        svFitProducerModule.config.event.likelihoodFunctions[0].verbosity = cms.int32(0)
        svFitProducerModule.config.event.srcPrimaryVertex = cms.InputTag('selectedPrimaryVertexByLeptonMatch')
        svFitProducerModule.config.event.builder = process.svFitEventBuilder
        svFitProducerModule.algorithm.markovChainOptions.initMode = cms.string("none")
        svFitProducerModule.algorithm.markovChainOptions.numIterBurnin = cms.uint32(numCalls*100)
        svFitProducerModule.algorithm.markovChainOptions.numIterSampling = cms.uint32(numCalls*1000)
        svFitProducerModule.algorithm.markovChainOptions.numIterSimAnnealingPhase1 = cms.uint32(numCalls*20)
        svFitProducerModule.algorithm.markovChainOptions.numIterSimAnnealingPhase2 = cms.uint32(numCalls*60)
        svFitProducerModule.algorithm.markovChainOptions.alpha = cms.double(1.0 - 0.1/numCalls)
        svFitProducerModule.algorithm.markovChainOptions.numChains = cms.uint32(1)
        svFitProducerModule.algorithm.markovChainOptions.numBatches = cms.uint32(1)
        svFitProducerModule.algorithm.markovChainOptions.epsilon0 = cms.double(1.e-2)
        ##svFitProducerModule.algorithm.monitorMarkovChain = cms.bool(True)
        svFitProducerModule.algorithm.monitorMarkovChain = cms.bool(False)
        svFitProducerModule.algorithm.monitorFilePath = cms.string('/data1/veelken/tmp/svFitStudies/')
        svFitProducerModule.algorithm.max_or_median = cms.string("max")
        svFitProducerModule.algorithm.pluginName = cms.string("svFitProducerByIntegration2WOtracks%sMax%ikCalls" % (metAlgorithmName, numCalls))
        svFitProducerModule.algorithm.verbosity = cms.int32(1)
        svFitProducerModuleName = "svFitProducerByIntegration2WOtracks%sMax%ikCalls" % (metAlgorithmName, numCalls)
        svFitAnalyzerModuleType = "SVfitEventHypothesisAnalyzer"

        svFitProducerModuleNames["MC%s%ik" % (metAlgorithmName, numCalls)] = svFitProducerModuleName
        svFitAnalyzerModuleTypes["MC%s%ik" % (metAlgorithmName, numCalls)] = svFitAnalyzerModuleType

        setattr(process, svFitProducerModuleName, svFitProducerModule)
        process.testSVfitTrackLikelihoodProductionSequence += svFitProducerModule
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# CV: run "stand-alone" version of SVfit for comparison

process.svFitStandaloneAnalyzerFit = cms.EDAnalyzer("SVfitStandaloneTestAnalyzer",
    doGenPlots = cms.bool(True),                                                     
    srcGenTauPairs = cms.InputTag(genTauPairs),
    srcGenLeg1 = cms.InputTag(srcGenLeg1),
    srcGenLeg2 = cms.InputTag(srcGenLeg2),
    srcGenMEt = cms.InputTag('genMetFromGenParticles'),
    doRecPlots = cms.bool(True),                                                       
    srcRecLeg1 = cms.InputTag(srcRecLeg1),
    srcRecLeg2 = cms.InputTag(srcRecLeg2),                                              
    srcRecMEt = cms.InputTag('pfType1CorrectedMet'),                                              
    srcRecMEtCov = cms.InputTag('pfMEtSignCovMatrix'),
    srcWeights = cms.VInputTag(),
    typeLeg1 = cms.string(svFitStandAloneTypeLeg1),
    typeLeg2 = cms.string(svFitStandAloneTypeLeg2),
    mode = cms.string("fit"),                                                     
    redoMEtCov = cms.bool(False),
    paramsMEtCov = process.svFitEventLikelihoodMEt2,
    fillHistograms = cms.bool(True),                                                     
    dqmDirectory = cms.string("svFitStandaloneAnalyzerFit"),
    verbosity = cms.int32(0)                                                 
)                                                  
process.testSVfitTrackLikelihoodProductionSequence += process.svFitStandaloneAnalyzerFit
process.svFitStandaloneAnalyzerInt = process.svFitStandaloneAnalyzerFit.clone(
    mode = cms.string("int"),
    dqmDirectory = cms.string("svFitStandaloneAnalyzerInt"),
    verbosity = cms.int32(0)
)
process.testSVfitTrackLikelihoodProductionSequence += process.svFitStandaloneAnalyzerInt
process.svFitStandaloneAnalyzerInt2 = process.svFitStandaloneAnalyzerFit.clone(
    mode = cms.string("int2"),
    dqmDirectory = cms.string("svFitStandaloneAnalyzerInt2"),
    verbosity = cms.int32(1)
)
process.testSVfitTrackLikelihoodProductionSequence += process.svFitStandaloneAnalyzerInt2
#--------------------------------------------------------------------------------

process.testSVfitTrackLikelihoodAnalysisSequence = cms.Sequence()

#--------------------------------------------------------------------------------
# fill histograms of mass reconstructed by running SVfit with different integration options

for option in svFitAnalyzerModuleTypes.keys():

    svFitProducerModuleName = svFitProducerModuleNames[option]
    svFitAnalyzerModuleType = svFitAnalyzerModuleTypes[option]
    
    svFitAnalyzerModule = cms.EDAnalyzer(svFitAnalyzerModuleType,
        srcEventHypotheses = cms.InputTag(svFitProducerModuleName),
        srcGenTauPairs = cms.InputTag(genTauPairs),
        srcGenLeg1 = cms.InputTag(srcGenLeg1),
        srcGenLeg2 = cms.InputTag(srcGenLeg2),
        srcGenMEt = cms.InputTag('genMetFromGenParticles'),
        srcPFMEtCovMatrix = cms.InputTag('pfMEtSignCovMatrix'),
        srcWeights = cms.VInputTag(),
        idxResonance = cms.int32(0),
        numBinsSVfitMass = cms.int32(500),
        svFitMassMax = cms.double(500.),
        numBinsSVfitSigma = cms.int32(250),
        svFitSigmaMax = cms.double(250.),
        dqmDirectory = cms.string(svFitProducerModuleName)
    )
    svFitAnalyzerModuleName = svFitProducerModuleName.replace("Producer", "Analyzer")
    setattr(process, svFitAnalyzerModuleName, svFitAnalyzerModule)
    process.testSVfitTrackLikelihoodAnalysisSequence += svFitAnalyzerModule
#--------------------------------------------------------------------------------

process.testSVfitTrackLikelihoodSequence = cms.Sequence(process.testSVfitTrackLikelihoodProductionSequence + process.testSVfitTrackLikelihoodAnalysisSequence)

#--------------------------------------------------------------------------------
# produce NeuralMtautau Ntuple

process.load("TauAnalysis/CandidateTools/neuralMtautau_cff")
process.neuralMtautauNtupleProducer.srcGenTauPair = cms.InputTag(genTauPairs)
process.neuralMtautauNtupleProducer.srcGenMEt = cms.InputTag('genMetFromGenParticles')
process.neuralMtautauNtupleProducer.srcGenTaus = cms.InputTag(genTaus)
process.neuralMtautauNtupleProducer.srcRecPFJets = cms.InputTag('ak5PFJetsNotOverlappingWithLeptons')
process.neuralMtautauNtupleProducer.srcRecPFCandidatesNotWithinJets = cms.InputTag('pfCandsNotInJetForPFMEtSignCovMatrix')
process.neuralMtautauNtupleProducer.srcRecLeg1 = cms.InputTag(srcRecLeg1)
process.neuralMtautauNtupleProducer.srcRecLeg2 = cms.InputTag(srcRecLeg2)
process.neuralMtautauNtupleProducer.srcRecMEt = cms.InputTag('pfType1CorrectedMet')
process.neuralMtautauNtupleProducer.srcPFMEtCovMatrix = cms.InputTag('pfMEtSignCovMatrix')
process.testSVfitTrackLikelihoodSequence += process.neuralMtautauNtupleProducer

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("/data1/veelken/tmp/svFitStudies/testSVfitTrackLikelihoods3_ntuple_%s_%s_%s_2012Dec16.root" % (sample_type, channel, massPoint))
)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# fill histograms of observables sensitive to tau polarization
# (to measure Higgs quantum numbers)

process.diTauPolarizationAnalyzer = cms.EDAnalyzer("DiCandidatePairPolarizationAnalyzer",
    srcGenParticles = cms.InputTag('genParticles'),
    srcRecLeg1 = cms.InputTag(srcRecLeg1),                                              
    srcRecLeg2 = cms.InputTag(srcRecLeg2),
    srcWeights = cms.VInputTag(),                                               
    dqmDirectory = cms.string("diTauPolarizationAnalyzer")
)
process.testSVfitTrackLikelihoodSequence += process.diTauPolarizationAnalyzer
#--------------------------------------------------------------------------------

process.p = cms.Path(process.testSVfitTrackLikelihoodSequence)

#--------------------------------------------------------------------------------
# save plots

process.DQMStore = cms.Service("DQMStore")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string("/data1/veelken/tmp/svFitStudies/testSVfitTrackLikelihoods3_plots_%s_%s_%s_2012Dec16.root" % (sample_type, channel, massPoint))
)

process.q = cms.EndPath(process.savePlots)
#--------------------------------------------------------------------------------

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

processDumpFile = open('testSVfitTrackLikelihoods3.dump' , 'w')
print >> processDumpFile, process.dumpPython()

