import FWCore.ParameterSet.Config as cms

process = cms.Process("testSVfitTrackLikelihoods2")

import os
import re

import TauAnalysis.Configuration.tools.castor as castor
import TauAnalysis.Configuration.tools.eos as eos
from TauAnalysis.Skimming.EventContent_cff import *

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load('Configuration/Geometry/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START53_V15::All')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2500)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        ##'file:/data1/veelken/CMSSW_5_2_x/skims/goldenZmumuEvents_ZplusJets_madgraph2_2012Apr12_AOD_9_1_cSC.root'
        'file:/tmp/veelken/genTauLeptonsPairAccSkim_ggHiggs125_mutau_2_1_4nD.root'                
    ),
    dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
    inputCommands=cms.untracked.vstring(
        'keep *',
        'drop LHERunInfoProduct_*_*_*',                        
        'drop LHEEventProduct_*_*_*'
    ),
    ##eventsToProcess = cms.untracked.VEventRange(
    ##    '1:1673:501726',
    ##    '1:1673:501684',
    ##    '1:1673:501676'                                
    ##)
)

#--------------------------------------------------------------------------------
# define configuration parameter default values

#sample_type = 'Z'
sample_type = 'Higgs'
#channel = 'etau'
channel = 'mutau'
#channel = 'emu'
#massPoint = '125'
massPoint = '300'
#qTmin = 50.
#qTmin = 20.
qTmin = -1.
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
    inputFilePath = '/store/user/veelken/CMSSW_5_3_x/skims/SVfitStudies/AccLowPtThresholds/ZplusJets_mutau'
    inputFile_regex = \
      r"[a-zA-Z0-9_/:.]*genTauLeptonPairSkim_ZplusJets_%s_(?P<gridJob>\d*)(_(?P<gridTry>\d*))*_(?P<hash>[a-zA-Z0-9]*).root" % channel
elif sample_type == 'Higgs':
    inputFilePath = '/store/user/veelken/CMSSW_5_3_x/skims/SVfitStudies/AccLowPtThresholds/ggPhi300_mutau/'
    inputFile_regex = \
      r"[a-zA-Z0-9_/:.]*genTauLeptonPairSkim_(ggHiggs|ggPhi|vbfHiggs)%s_%s_(?P<gridJob>\d*)(_(?P<gridTry>\d*))*_(?P<hash>[a-zA-Z0-9]*).root" % (massPoint, channel)
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
process.testSVfitTrackLikelihoodProductionSequence += process.PFTau

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

process.boostedDiTaus = cms.EDFilter("CandViewSelector",
    src = cms.InputTag(genTauPairs),
    cut = cms.string("pt > %f" % qTmin)
)
process.testSVfitTrackLikelihoodProductionSequence += process.boostedDiTaus
process.boostedDiTauFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag('boostedDiTaus'),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1000)                      
)
process.testSVfitTrackLikelihoodProductionSequence += process.boostedDiTauFilter
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

process.svFitTauToMuBuilder.fixToGenVisEnFracX = cms.bool(False)
process.svFitTauToMuBuilder.fixToGenPhiLab = cms.bool(False)
process.svFitTauToMuBuilder.fixToGenNuInvMass = cms.bool(False)
process.svFitTauToMuBuilder.fixToGenDeltaR = cms.bool(False)
process.svFitTauToMuBuilder.fixToGenVisP4 = cms.bool(False)
process.svFitTauToMuBuilder.initializeToGen = cms.bool(False)
#process.svFitTauToMuBuilder.srcGenTaus = cms.InputTag('genParticles')
process.svFitTauToMuBuilder.dRmatch = cms.double(0.3)

process.svFitTauToHadBuilder.fixToGenVisEnFracX = cms.bool(False)
process.svFitTauToHadBuilder.fixToGenPhiLab = cms.bool(False)
process.svFitTauToHadBuilder.fixToGenDeltaR = cms.bool(False)
process.svFitTauToHadBuilder.fixToGenVisP4 = cms.bool(False)
process.svFitTauToHadBuilder.initializeToGen = cms.bool(False)
#process.svFitTauToHadBuilder.srcGenTaus = cms.InputTag('genParticles')
process.svFitTauToHadBuilder.dRmatch = cms.double(0.3)

# CV: fix event vertex position to Monte Carlo truth value
process.svFitEventBuilder.fixToGenVertex = cms.bool(False)
#process.svFitEventBuilder.srcGenVertex = cms.InputTag('genEventVertex')

process.svFitElectronLikelihoodTrackInfo = cms.PSet(
    pluginName = cms.string("svFitTauToElecLikelihoodTrackInfo"),
    pluginType = cms.string("SVfitTauDecayLikelihoodTrackInfo"),
    useLifetimeConstraint = cms.bool(True),
    sfProdVertexCov = cms.double(2.0),
    sfDecayVertexCov = cms.double(2.0),
    verbosity = cms.int32(0)  
)

process.svFitMuonLikelihoodTrackInfo = cms.PSet(
    pluginName = cms.string("svFitTauToMuLikelihoodTrackInfo"),
    pluginType = cms.string("SVfitTauDecayLikelihoodTrackInfo"),
    useLifetimeConstraint = cms.bool(True),
    sfProdVertexCov = cms.double(2.0),
    sfDecayVertexCov = cms.double(2.0),
    verbosity = cms.int32(0)  
)

process.svFitTauLikelihoodTrackInfo = cms.PSet(
    pluginName = cms.string("svFitTauToHadLikelihoodTrackInfo"),
    pluginType = cms.string("SVfitTauDecayLikelihoodTrackInfo"),
    useLifetimeConstraint = cms.bool(True),
    sfProdVertexCov = cms.double(2.0),
    sfDecayVertexCov = cms.double(2.0),
    verbosity = cms.int32(0)  
)

process.svFitEventLikelihoodTrackInfo = cms.PSet(
    pluginName = cms.string("svFitEventLikelihoodTrackInfo"),
    pluginType = cms.string("SVfitEventLikelihoodTrackInfo"),
    verbosity = cms.int32(0)
)

srcGenLeg1                      = None
srcRecLeg1                      = None
svFitLikelihoodLeg1_kinematics = None
svFitLikelihoodLeg1_trackinfo  = None
svFitBuilderLeg1               = None
svFitStandAloneTypeLeg1        = None
srcGenLeg2                      = None
srcRecLeg2                      = None
svFitLikelihoodLeg2_kinematics = None
svFitLikelihoodLeg2_trackinfo  = None
svFitBuilderLeg2               = None
svFitStandAloneTypeLeg2        = None
if channel == 'mutau':
    srcGenLeg1                      = genMuonsFromTauDecays
    srcRecLeg1                      = 'genMatchedPatMuons'
    svFitLikelihoodLeg1_kinematics = process.svFitMuonLikelihoodMatrixElement
    svFitLikelihoodLeg1_trackinfo  = process.svFitMuonLikelihoodTrackInfo
    svFitBuilderLeg1               = process.svFitTauToMuBuilder
    svFitStandAloneTypeLeg1        = "lep"
    srcGenLeg2                      = genTauJetsFromTauDecays
    srcRecLeg2                      = 'genMatchedPatTaus'
    svFitLikelihoodLeg2_kinematics = process.svFitTauLikelihoodPhaseSpace
    svFitLikelihoodLeg2_trackinfo  = process.svFitTauLikelihoodTrackInfo
    svFitBuilderLeg2               = process.svFitTauToHadBuilder
    svFitStandAloneTypeLeg2        = "had"
elif channel == 'etau':
    srcGenLeg1                      = genElectronsFromTauDecays
    srcRecLeg1                      = 'genMatchedPatElectrons'
    svFitLikelihoodLeg1_kinematics = process.svFitElectronLikelihoodMatrixElement
    svFitLikelihoodLeg1_trackinfo  = process.svFitElectronLikelihoodTrackInfo
    svFitBuilderLeg1               = process.svFitTauToElecBuilder
    svFitStandAloneTypeLeg1        = "lep"
    srcGenLeg2                      = genTauJetsFromTauDecays
    srcRecLeg2                      = 'genMatchedPatTaus'
    svFitLikelihoodLeg2_kinematics = process.svFitTauLikelihoodPhaseSpace
    svFitLikelihoodLeg2_trackinfo  = process.svFitTauLikelihoodTrackInfo
    svFitBuilderLeg2               = process.svFitTauToHadBuilder
    svFitStandAloneTypeLeg2        = "had"
elif channel == 'emu':
    srcGenLeg1                      = genElectronsFromTauDecays
    srcRecLeg1                      = 'genMatchedPatElectrons'
    svFitLikelihoodLeg1_kinematics = process.svFitElectronLikelihoodMatrixElement
    svFitLikelihoodLeg1_trackinfo  = process.svFitElectronLikelihoodTrackInfo
    svFitBuilderLeg1               = process.svFitTauToElecBuilder
    svFitStandAloneTypeLeg1        = "lep"
    srcGenLeg2                      = genMuonsFromTauDecays
    srcRecLeg2                      = 'genMatchedPatMuons'
    svFitLikelihoodLeg2_kinematics = process.svFitMuonLikelihoodMatrixElement
    svFitLikelihoodLeg2_trackinfo  = process.svFitMuonLikelihoodTrackInfo
    svFitBuilderLeg2               = process.svFitTauToMuBuilder
    svFitStandAloneTypeLeg2        = "lep"
elif channel == 'ditau':
    srcGenLeg1                      = genTauJetsFromTauDecays
    srcRecLeg1                      = 'genMatchedPatTaus'
    svFitLikelihoodLeg1_kinematics = process.svFitTauLikelihoodPhaseSpace
    svFitLikelihoodLeg1_trackinfo  = process.svFitTauLikelihoodTrackInfo
    svFitBuilderLeg1               = process.svFitTauToHadBuilder
    svFitStandAloneTypeLeg1        = "had"
    srcGenLeg2                      = genTauJetsFromTauDecays
    srcRecLeg2                      = 'genMatchedPatTaus'
    svFitLikelihoodLeg2_kinematics = process.svFitTauLikelihoodPhaseSpace
    svFitLikelihoodLeg2_trackinfo  = process.svFitTauLikelihoodTrackInfo
    svFitBuilderLeg2               = process.svFitTauToHadBuilder
    svFitStandAloneTypeLeg2        = "had"
else:
    raise ValueError("Invalid channel = %s !!" % channel)

svFitProducerModuleNames = dict()
svFitAnalyzerModuleTypes = dict()

for option in [ '5a', '5b' ]:
        
    svFitProducerModule = None
    svFitProducerModuleName = None
    svFitAnalyzerModuleType = None
    
    if option == '1':
        # option 1: VEGAS integration of likelihood functions, no tracking information used
        svFitProducerModule = process.svFitProducerByIntegration.clone()
        svFitProducerModule.config.event.resonances.A.daughters.leg1.src = cms.InputTag(srcRecLeg1)
        svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions = cms.VPSet(svFitLikelihoodLeg1_kinematics)
        svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
        svFitProducerModule.config.event.resonances.A.daughters.leg1.builder = svFitBuilderLeg1
        svFitProducerModule.config.event.resonances.A.daughters.leg2.src = cms.InputTag(srcRecLeg2)
        svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions = cms.VPSet(svFitLikelihoodLeg2_kinematics)
        svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
        svFitProducerModule.config.event.resonances.A.daughters.leg2.builder = svFitBuilderLeg2
        svFitProducerModule.config.event.resonances.A.likelihoodFunctions = cms.VPSet()
        svFitProducerModule.config.event.srcMEt = cms.InputTag('pfType1CorrectedMet')
        svFitProducerModule.config.event.srcPrimaryVertex = cms.InputTag('selectedPrimaryVertexByLeptonMatch')
        svFitProducerModule.algorithm.pluginName = cms.string("svFitProducerByIntegrationWOtracks")
        svFitProducerModule.algorithm.verbosity = cms.int32(1)
        svFitProducerModuleName = "svFitProducerByIntegrationWOtracks"
        svFitAnalyzerModuleType = "SVfitEventHypothesisByIntegrationAnalyzer"
    elif option == '2':
        if channel in [ "mutau", "etau", "ditau" ]:
            # option 2: VEGAS integration of likelihood functions,
            #           tracking information used for 3-prongs only,
            #           usage of tau life-time information enabled
            svFitProducerModule = process.svFitProducerByIntegration.clone()
            svFitProducerModule.config.event.resonances.A.daughters.leg1.src = cms.InputTag(srcRecLeg1)
            if channel in [ "ditau" ]:
                svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions = cms.VPSet(svFitLikelihoodLeg1_kinematics, svFitLikelihoodLeg1_trackinfo.clone())
                svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
                svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions[1].ignore3Prongs = cms.bool(False)
                svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions[1].ignore1Prongs = cms.bool(True)
            else:
                svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions = cms.VPSet(svFitLikelihoodLeg1_kinematics)
                svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
            svFitProducerModule.config.event.resonances.A.daughters.leg1.builder = svFitBuilderLeg1
            svFitProducerModule.config.event.resonances.A.daughters.leg2.src = cms.InputTag(srcRecLeg2)
            svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions = cms.VPSet(svFitLikelihoodLeg2_kinematics, svFitLikelihoodLeg2_trackinfo.clone())
            svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
            svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[1].ignore3Prongs = cms.bool(False)
            svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[1].ignore1Prongs = cms.bool(True)
            svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[1].useLifetimeConstraint = cms.bool(True)
            svFitProducerModule.config.event.resonances.A.daughters.leg2.builder = svFitBuilderLeg2
            svFitProducerModule.config.event.resonances.A.likelihoodFunctions = cms.VPSet()
            ## CV: varying primary event vertex position makes algorithm numerically unstable !!
            ##svFitProducerModule.config.event.likelihoodFunctions = cms.VPSet(process.svFitEventLikelihoodMEt2, process.svFitEventLikelihoodTrackInfo)
            svFitProducerModule.config.event.likelihoodFunctions = cms.VPSet(process.svFitEventLikelihoodMEt2)
            svFitProducerModule.config.event.srcMEt = cms.InputTag('pfType1CorrectedMet')
            svFitProducerModule.config.event.srcPrimaryVertex = cms.InputTag('selectedPrimaryVertexByLeptonMatch')
            svFitProducerModule.algorithm.pluginName = cms.string("svFitProducerByIntegrationWtracksFor3prongsOnlyWlifetimeConstraint")
            svFitProducerModule.algorithm.verbosity = cms.int32(1)
            svFitProducerModuleName = "svFitProducerByIntegrationWtracksFor3prongsOnlyWlifetimeConstraint"
            svFitAnalyzerModuleType = "SVfitEventHypothesisByIntegrationAnalyzer"
    elif option == '3':
        if channel in [ "mutau", "etau", "ditau" ]:
            # option 3: VEGAS integration of likelihood functions,
            #           tracking information used for 3-prongs only,
            #           usage of tau life-time information disabled
            svFitProducerModule = process.svFitProducerByIntegration.clone()
            svFitProducerModule.config.event.resonances.A.daughters.leg1.src = cms.InputTag(srcRecLeg1)
            if channel in [ "diTau" ]:
                svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions = cms.VPSet(svFitLikelihoodLeg1_kinematics, svFitLikelihoodLeg1_trackinfo.clone())
                svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
                svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions[1].ignore3Prongs = cms.bool(False)
                svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions[1].ignore1Prongs = cms.bool(True)
            else:
                svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions = cms.VPSet(svFitLikelihoodLeg1_kinematics)
                svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
            svFitProducerModule.config.event.resonances.A.daughters.leg1.builder = svFitBuilderLeg1
            svFitProducerModule.config.event.resonances.A.daughters.leg2.src = cms.InputTag(srcRecLeg2)
            svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions = cms.VPSet(svFitLikelihoodLeg2_kinematics, svFitLikelihoodLeg2_trackinfo.clone())
            svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
            svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[1].ignore3Prongs = cms.bool(False)
            svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[1].ignore1Prongs = cms.bool(True)
            svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[1].useLifetimeConstraint = cms.bool(False)
            svFitProducerModule.config.event.resonances.A.daughters.leg2.builder = svFitBuilderLeg2
            svFitProducerModule.config.event.resonances.A.likelihoodFunctions = cms.VPSet()
            ## CV: varying primary event vertex position makes algorithm numerically unstable !!
            ##svFitProducerModule.config.event.likelihoodFunctions = cms.VPSet(process.svFitEventLikelihoodMEt2, process.svFitEventLikelihoodTrackInfo)
            svFitProducerModule.config.event.likelihoodFunctions = cms.VPSet(process.svFitEventLikelihoodMEt2)
            svFitProducerModule.config.event.srcMEt = cms.InputTag('pfType1CorrectedMet')
            svFitProducerModule.config.event.srcPrimaryVertex = cms.InputTag('selectedPrimaryVertexByLeptonMatch')
            svFitProducerModule.algorithm.pluginName = cms.string("svFitProducerByIntegrationWtracksFor3prongsOnlyWOlifetimeConstraint")
            svFitProducerModule.algorithm.verbosity = cms.int32(1)
            svFitProducerModuleName = "svFitProducerByIntegrationWtracksFor3prongsOnlyWOlifetimeConstraint"
            svFitAnalyzerModuleType = "SVfitEventHypothesisByIntegrationAnalyzer"        
    elif option == '4':
        # option 4: VEGAS integration of likelihood functions,
        #           tracking information used for all tau decays
        svFitProducerModule = process.svFitProducerByIntegration.clone()
        svFitProducerModule.config.event.resonances.A.daughters.leg1.src = cms.InputTag(srcRecLeg1)
        svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions = cms.VPSet(svFitLikelihoodLeg1_kinematics, svFitLikelihoodLeg1_trackinfo)
        svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
        svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions[1].ignore3Prongs = cms.bool(False)
        svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions[1].ignore1Prongs = cms.bool(False)
        svFitProducerModule.config.event.resonances.A.daughters.leg1.builder = svFitBuilderLeg1
        svFitProducerModule.config.event.resonances.A.daughters.leg2.src = cms.InputTag(srcRecLeg2)
        svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions = cms.VPSet(svFitLikelihoodLeg2_kinematics, svFitLikelihoodLeg2_trackinfo)
        svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
        svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[1].ignore3Prongs = cms.bool(False)
        svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[1].ignore1Prongs = cms.bool(False)
        svFitProducerModule.config.event.resonances.A.daughters.leg2.builder = svFitBuilderLeg2
        svFitProducerModule.config.event.resonances.A.likelihoodFunctions = cms.VPSet()
        ## CV: varying primary event vertex position makes algorithm numerically unstable !!
        ##svFitProducerModule.config.event.likelihoodFunctions = cms.VPSet(process.svFitEventLikelihoodMEt2, process.svFitEventLikelihoodTrackInfo)
        svFitProducerModule.config.event.likelihoodFunctions = cms.VPSet(process.svFitEventLikelihoodMEt2)
        svFitProducerModule.config.event.srcMEt = cms.InputTag('pfType1CorrectedMet')
        svFitProducerModule.config.event.srcPrimaryVertex = cms.InputTag('selectedPrimaryVertexByLeptonMatch')
        svFitProducerModule.algorithm.pluginName = cms.string("svFitProducerByIntegrationWtracksFor1prongsAnd3prongs")
        svFitProducerModule.algorithm.verbosity = cms.int32(1)
        svFitProducerModuleName = "svFitProducerByIntegrationWtracksFor1prongsAnd3prongs"
        svFitAnalyzerModuleType = "SVfitEventHypothesisByIntegrationAnalyzer"        
    elif option == '5a' or option == '5b':
        # option 5: Markov Chain integration of likelihood functions, no tracking information used
        svFitLikelihoodLeg1_kinematics_cloned = svFitLikelihoodLeg1_kinematics.clone()
        svFitLikelihoodLeg2_kinematics_cloned = svFitLikelihoodLeg2_kinematics.clone()
        if channel == 'mutau' or channel == 'etau':
            if option == '5a':
                svFitLikelihoodLeg2_kinematics_cloned.applyVisMassFactor = cms.bool(False)
            elif option == '5b':
                svFitLikelihoodLeg2_kinematics_cloned.applyVisMassFactor = cms.bool(True)
        elif channel == 'ditau':
            if option == '5a':
                svFitLikelihoodLeg1_kinematics_cloned.applyVisMassFactor = cms.bool(False)
                svFitLikelihoodLeg2_kinematics_cloned.applyVisMassFactor = cms.bool(False)
            elif option == '5b':
                svFitLikelihoodLeg1_kinematics_cloned.applyVisMassFactor = cms.bool(True)
                svFitLikelihoodLeg2_kinematics_cloned.applyVisMassFactor = cms.bool(True)
        svFitProducerModule = process.svFitProducerByIntegration2.clone()
        svFitProducerModule.config.event.resonances.A.daughters.leg1.src = cms.InputTag(srcRecLeg1)
        svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions = cms.VPSet(svFitLikelihoodLeg1_kinematics_cloned)
        svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
        svFitProducerModule.config.event.resonances.A.daughters.leg1.builder = svFitBuilderLeg1
        svFitProducerModule.config.event.resonances.A.daughters.leg2.src = cms.InputTag(srcRecLeg2)
        svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions = cms.VPSet(svFitLikelihoodLeg2_kinematics_cloned)
        svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
        svFitProducerModule.config.event.resonances.A.daughters.leg2.builder = svFitBuilderLeg2
        svFitProducerModule.config.event.resonances.A.likelihoodFunctions = cms.VPSet()
        svFitProducerModule.config.event.likelihoodFunctions = cms.VPSet(process.svFitEventLikelihoodMEt2)
        ##svFitProducerModule.config.event.likelihoodFunctions[0].monitorMEtUncertainty = cms.bool(True)
        ##svFitProducerModule.config.event.likelihoodFunctions[0].verbosity = cms.int32(1)
        svFitProducerModule.config.event.likelihoodFunctions[0].monitorMEtUncertainty = cms.bool(False)
        svFitProducerModule.config.event.likelihoodFunctions[0].verbosity = cms.int32(0)
        svFitProducerModule.config.event.srcMEt = cms.InputTag('pfType1CorrectedMet')
        svFitProducerModule.config.event.srcPrimaryVertex = cms.InputTag('selectedPrimaryVertexByLeptonMatch')
        svFitProducerModule.config.event.builder = process.svFitEventBuilder
        svFitProducerModule.algorithm.markovChainOptions.initMode = cms.string("none")
        svFitProducerModule.algorithm.markovChainOptions.numIterBurnin = cms.uint32(50000)
        svFitProducerModule.algorithm.markovChainOptions.numIterSampling = cms.uint32(500000)
        svFitProducerModule.algorithm.markovChainOptions.numIterSimAnnealingPhase1 = cms.uint32(10000)
        svFitProducerModule.algorithm.markovChainOptions.numIterSimAnnealingPhase2 = cms.uint32(30000)
        svFitProducerModule.algorithm.markovChainOptions.alpha = cms.double(0.9999)
        svFitProducerModule.algorithm.markovChainOptions.numChains = cms.uint32(1)
        svFitProducerModule.algorithm.markovChainOptions.numBatches = cms.uint32(1)
        svFitProducerModule.algorithm.markovChainOptions.epsilon0 = cms.double(1.e-2)
        svFitProducerModule.algorithm.monitorMarkovChain = cms.bool(False)
        svFitProducerModule.algorithm.max_or_median = cms.string("max")
        if option == '5a':
            svFitProducerModule.algorithm.pluginName = cms.string("svFitProducerByIntegration2WOtracksMaxWOvisMassFactor")
        elif option == '5b':
            svFitProducerModule.algorithm.pluginName = cms.string("svFitProducerByIntegration2WOtracksMaxWvisMassFactor")
        svFitProducerModule.algorithm.verbosity = cms.int32(1)
        if option == '5a':
            svFitProducerModuleName = "svFitProducerByIntegration2WOtracksMaxWOvisMassFactor"
        elif option == '5b':
            svFitProducerModuleName = "svFitProducerByIntegration2WOtracksMaxWvisMassFactor"
        svFitAnalyzerModuleType = "SVfitEventHypothesisAnalyzer"
    elif option == '6':
        if channel in [ "mutau", "etau", "ditau" ]:
            # option 6: Markov Chain integration of likelihood functions, tracking information used for 3-prongs only,
            #           usage of tau life-time information enabled
            svFitProducerModule = process.svFitProducerByIntegration2.clone()
            svFitProducerModule.config.event.resonances.A.daughters.leg1.src = cms.InputTag(srcRecLeg1)
            if channel in [ "ditau" ]:
                svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions = cms.VPSet(svFitLikelihoodLeg1_kinematics, svFitLikelihoodLeg1_trackinfo.clone())
                svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
                svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions[1].ignore3Prongs = cms.bool(False)
                svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions[1].ignore1Prongs = cms.bool(True)
            else:
                svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions = cms.VPSet(svFitLikelihoodLeg1_kinematics)
                svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
            svFitProducerModule.config.event.resonances.A.daughters.leg1.builder = svFitBuilderLeg1
            svFitProducerModule.config.event.resonances.A.daughters.leg2.src = cms.InputTag(srcRecLeg2)
            svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions = cms.VPSet(svFitLikelihoodLeg2_kinematics, svFitLikelihoodLeg2_trackinfo.clone())
            svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
            svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[1].ignore3Prongs = cms.bool(False)
            svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[1].ignore1Prongs = cms.bool(True)
            svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[1].useLifetimeConstraint = cms.bool(True)
            svFitProducerModule.config.event.resonances.A.daughters.leg2.builder = svFitBuilderLeg2
            svFitProducerModule.config.event.resonances.A.likelihoodFunctions = cms.VPSet()
            ## CV: varying primary event vertex position makes algorithm numerically unstable !!
            ##svFitProducerModule.config.event.likelihoodFunctions = cms.VPSet(process.svFitEventLikelihoodMEt2, process.svFitEventLikelihoodTrackInfo)
            svFitProducerModule.config.event.likelihoodFunctions = cms.VPSet(process.svFitEventLikelihoodMEt2)
            svFitProducerModule.config.event.srcMEt = cms.InputTag('pfType1CorrectedMet')
            svFitProducerModule.config.event.srcPrimaryVertex = cms.InputTag('selectedPrimaryVertexByLeptonMatch')
            svFitProducerModule.config.event.builder = process.svFitEventBuilder
            svFitProducerModule.algorithm.markovChainOptions.initMode = cms.string("none")
            svFitProducerModule.algorithm.markovChainOptions.numIterBurnin = cms.uint32(50000)
            svFitProducerModule.algorithm.markovChainOptions.numIterSampling = cms.uint32(500000)
            svFitProducerModule.algorithm.markovChainOptions.numIterSimAnnealingPhase1 = cms.uint32(10000)
            svFitProducerModule.algorithm.markovChainOptions.numIterSimAnnealingPhase2 = cms.uint32(30000)
            svFitProducerModule.algorithm.markovChainOptions.alpha = cms.double(0.9999)
            svFitProducerModule.algorithm.markovChainOptions.numChains = cms.uint32(1)
            svFitProducerModule.algorithm.markovChainOptions.numBatches = cms.uint32(1)
            svFitProducerModule.algorithm.markovChainOptions.epsilon0 = cms.double(1.e-2)
            svFitProducerModule.algorithm.monitorMarkovChain = cms.bool(False)
            svFitProducerModule.algorithm.pluginName = cms.string("svFitProducerByIntegration2WtracksFor3prongsOnlyWlifetimeConstraint")
            svFitProducerModule.algorithm.verbosity = cms.int32(1)
            svFitProducerModuleName = "svFitProducerByIntegration2WtracksFor3prongsOnlyWlifetimeConstraint"
            svFitAnalyzerModuleType = "SVfitEventHypothesisAnalyzer"
    elif option == '7':
        if channel in [ "mutau", "etau", "ditau" ]:
            # option 7: Markov Chain integration of likelihood functions, tracking information used for 3-prongs only,
            #           usage of tau life-time information disabled
            svFitProducerModule = process.svFitProducerByIntegration2.clone()
            svFitProducerModule.config.event.resonances.A.daughters.leg1.src = cms.InputTag(srcRecLeg1)
            svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions = cms.VPSet(svFitLikelihoodLeg1_kinematics)
            svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
            svFitProducerModule.config.event.resonances.A.daughters.leg1.builder = svFitBuilderLeg1
            svFitProducerModule.config.event.resonances.A.daughters.leg2.src = cms.InputTag(srcRecLeg2)
            svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions = cms.VPSet(svFitLikelihoodLeg2_kinematics, svFitLikelihoodLeg2_trackinfo.clone())
            svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
            svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[1].ignore3Prongs = cms.bool(False)
            svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[1].ignore1Prongs = cms.bool(True)
            svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[1].useLifetimeConstraint = cms.bool(False)
            svFitProducerModule.config.event.resonances.A.daughters.leg2.builder = svFitBuilderLeg2
            svFitProducerModule.config.event.resonances.A.likelihoodFunctions = cms.VPSet()
            ## CV: varying primary event vertex position makes algorithm numerically unstable !!
            ##svFitProducerModule.config.event.likelihoodFunctions = cms.VPSet(process.svFitEventLikelihoodMEt2, process.svFitEventLikelihoodTrackInfo)
            svFitProducerModule.config.event.likelihoodFunctions = cms.VPSet(process.svFitEventLikelihoodMEt2)
            svFitProducerModule.config.event.srcMEt = cms.InputTag('pfType1CorrectedMet')
            svFitProducerModule.config.event.srcPrimaryVertex = cms.InputTag('selectedPrimaryVertexByLeptonMatch')
            svFitProducerModule.config.event.builder = process.svFitEventBuilder
            svFitProducerModule.algorithm.markovChainOptions.initMode = cms.string("none")
            svFitProducerModule.algorithm.markovChainOptions.numIterBurnin = cms.uint32(50000)
            svFitProducerModule.algorithm.markovChainOptions.numIterSampling = cms.uint32(500000)
            svFitProducerModule.algorithm.markovChainOptions.numIterSimAnnealingPhase1 = cms.uint32(10000)
            svFitProducerModule.algorithm.markovChainOptions.numIterSimAnnealingPhase2 = cms.uint32(30000)
            svFitProducerModule.algorithm.markovChainOptions.alpha = cms.double(0.9999)
            svFitProducerModule.algorithm.markovChainOptions.numChains = cms.uint32(1)
            svFitProducerModule.algorithm.markovChainOptions.numBatches = cms.uint32(1)
            svFitProducerModule.algorithm.markovChainOptions.epsilon0 = cms.double(1.e-2)
            svFitProducerModule.algorithm.monitorMarkovChain = cms.bool(False)
            svFitProducerModule.algorithm.pluginName = cms.string("svFitProducerByIntegration2WtracksFor3prongsOnlyWOlifetimeConstraint")
            svFitProducerModule.algorithm.verbosity = cms.int32(1)
            svFitProducerModuleName = "svFitProducerByIntegration2WtracksFor3prongsOnlyWOlifetimeConstraint"
            svFitAnalyzerModuleType = "SVfitEventHypothesisAnalyzer"        
    elif option == '8':
        # option 8: Markov Chain integration of likelihood functions, tracking information used for all tau decays
        svFitProducerModule = process.svFitProducerByIntegration2.clone()
        svFitProducerModule.config.event.resonances.A.daughters.leg1.src = cms.InputTag(srcRecLeg1)
        svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions = cms.VPSet(svFitLikelihoodLeg1_kinematics, svFitLikelihoodLeg1_trackinfo)
        svFitProducerModule.config.event.resonances.A.daughters.leg1.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
        svFitProducerModule.config.event.resonances.A.daughters.leg1.builder = svFitBuilderLeg1
        svFitProducerModule.config.event.resonances.A.daughters.leg2.src = cms.InputTag(srcRecLeg2)
        svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions = cms.VPSet(svFitLikelihoodLeg2_kinematics, svFitLikelihoodLeg2_trackinfo)
        svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
        svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[1].ignore3Prongs = cms.bool(False)
        svFitProducerModule.config.event.resonances.A.daughters.leg2.likelihoodFunctions[1].ignore1Prongs = cms.bool(False)
        svFitProducerModule.config.event.resonances.A.daughters.leg2.builder = svFitBuilderLeg2
        svFitProducerModule.config.event.resonances.A.likelihoodFunctions = cms.VPSet()
        ## CV: varying primary event vertex position makes algorithm numerically unstable !!
        ##svFitProducerModule.config.event.likelihoodFunctions = cms.VPSet(process.svFitEventLikelihoodMEt2, process.svFitEventLikelihoodTrackInfo)
        svFitProducerModule.config.event.likelihoodFunctions = cms.VPSet(process.svFitEventLikelihoodMEt2)
        svFitProducerModule.config.event.srcMEt = cms.InputTag('pfType1CorrectedMet')
        svFitProducerModule.config.event.srcPrimaryVertex = cms.InputTag('selectedPrimaryVertexByLeptonMatch')
        svFitProducerModule.config.event.builder = process.svFitEventBuilder
        svFitProducerModule.algorithm.markovChainOptions.initMode = cms.string("none")
        svFitProducerModule.algorithm.markovChainOptions.numIterBurnin = cms.uint32(50000)
        svFitProducerModule.algorithm.markovChainOptions.numIterSampling = cms.uint32(500000)
        svFitProducerModule.algorithm.markovChainOptions.numIterSimAnnealingPhase1 = cms.uint32(10000)
        svFitProducerModule.algorithm.markovChainOptions.numIterSimAnnealingPhase2 = cms.uint32(30000)
        svFitProducerModule.algorithm.markovChainOptions.alpha = cms.double(0.9999)
        svFitProducerModule.algorithm.markovChainOptions.numChains = cms.uint32(1)
        svFitProducerModule.algorithm.markovChainOptions.numBatches = cms.uint32(1)
        svFitProducerModule.algorithm.markovChainOptions.epsilon0 = cms.double(1.e-2)
        svFitProducerModule.algorithm.monitorMarkovChain = cms.bool(False)
        svFitProducerModule.algorithm.max_or_median = cms.string("max")
        svFitProducerModule.algorithm.pluginName = cms.string("svFitProducerByIntegration2WtracksFor1prongsAnd3prongsMax")
        svFitProducerModule.algorithm.verbosity = cms.int32(1)
        svFitProducerModuleName = "svFitProducerByIntegration2WtracksFor1prongsAnd3prongsMax"
        svFitAnalyzerModuleType = "SVfitEventHypothesisAnalyzer"
    elif option == '9':
        # option 9: Markov Chain integration of likelihood functions, no tracking information used;
        #           new SVfitEventLikelihoodMEt3 used
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
        svFitProducerModule.config.event.likelihoodFunctions = cms.VPSet(process.svFitEventLikelihoodMEt3)
        ##svFitProducerModule.config.event.likelihoodFunctions[0].monitorMEtUncertainty = cms.bool(True)
        ##svFitProducerModule.config.event.likelihoodFunctions[0].verbosity = cms.int32(1)
        svFitProducerModule.config.event.likelihoodFunctions[0].monitorMEtUncertainty = cms.bool(False)
        svFitProducerModule.config.event.likelihoodFunctions[0].verbosity = cms.int32(0)
        svFitProducerModule.config.event.srcMEt = cms.InputTag('pfType1CorrectedMet')
        svFitProducerModule.config.event.srcPrimaryVertex = cms.InputTag('selectedPrimaryVertexByLeptonMatch')
        svFitProducerModule.config.event.builder = process.svFitEventBuilder
        svFitProducerModule.algorithm.markovChainOptions.initMode = cms.string("none")
        svFitProducerModule.algorithm.markovChainOptions.numIterBurnin = cms.uint32(50000)
        svFitProducerModule.algorithm.markovChainOptions.numIterSampling = cms.uint32(500000)
        svFitProducerModule.algorithm.markovChainOptions.numIterSimAnnealingPhase1 = cms.uint32(10000)
        svFitProducerModule.algorithm.markovChainOptions.numIterSimAnnealingPhase2 = cms.uint32(30000)
        svFitProducerModule.algorithm.markovChainOptions.alpha = cms.double(0.9999)
        svFitProducerModule.algorithm.markovChainOptions.numChains = cms.uint32(1)
        svFitProducerModule.algorithm.markovChainOptions.numBatches = cms.uint32(1)
        svFitProducerModule.algorithm.markovChainOptions.epsilon0 = cms.double(1.e-2)
        svFitProducerModule.algorithm.monitorMarkovChain = cms.bool(False)
        svFitProducerModule.algorithm.max_or_median = cms.string("max")
        svFitProducerModule.algorithm.pluginName = cms.string("svFitProducerByIntegration2WOtracksMaxNewMEtUncertainty")
        svFitProducerModule.algorithm.verbosity = cms.int32(1)
        svFitProducerModuleName = "svFitProducerByIntegration2WOtracksMaxNewMEtUncertainty"
        svFitAnalyzerModuleType = "SVfitEventHypothesisAnalyzer"
    else:
        raise ValueError("Option %i undefined !!" % option)

    if svFitProducerModule:
        svFitProducerModuleNames[option] = svFitProducerModuleName
        svFitAnalyzerModuleTypes[option] = svFitAnalyzerModuleType

        setattr(process, svFitProducerModuleName, svFitProducerModule)
        process.testSVfitTrackLikelihoodProductionSequence += svFitProducerModule
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# CV: run "stand-alone" version of SVfit for comparison

process.svFitStandaloneAnalyzer = cms.EDAnalyzer("SVfitStandaloneTestAnalyzer",
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
    dqmDirectory = cms.string("svFitStandaloneAnalyzer"),
    verbosity = cms.int32(0)                                                 
)                                          
process.testSVfitTrackLikelihoodProductionSequence += process.svFitStandaloneAnalyzer
#--------------------------------------------------------------------------------

process.testSVfitTrackLikelihoodAnalysisSequence = cms.Sequence()

#--------------------------------------------------------------------------------
# fill histograms of mass reconstructed by running SVfit with and without Track likelihoods

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
        numBinsSVfitMass = cms.int32(2500),
        svFitMassMax = cms.double(2500.),
        numBinsSVfitSigma = cms.int32(250),
        svFitSigmaMax = cms.double(250.),
        dqmDirectory = cms.string(svFitProducerModuleName)
    )
    svFitAnalyzerModuleName = svFitProducerModuleName.replace("Producer", "Analyzer")
    setattr(process, svFitAnalyzerModuleName, svFitAnalyzerModule)
    process.testSVfitTrackLikelihoodAnalysisSequence += svFitAnalyzerModule
#--------------------------------------------------------------------------------

process.testSVfitTrackLikelihoodSequence = cms.Sequence(process.testSVfitTrackLikelihoodProductionSequence + process.testSVfitTrackLikelihoodAnalysisSequence)

process.p = cms.Path(process.testSVfitTrackLikelihoodSequence)

#--------------------------------------------------------------------------------
# make plots of SVfit input variables, separetely for 
# events in which SVfit mass reconstructed using track information is better/worse
# than SVfit mass reconstructed without using track information
# in likelihood model

for option1_vs_2 in [ [ '1', '2' ], [ '1', '3' ], [ '1', '4' ], [ '5a', '6' ], [ '5a', '7' ], [ '5a', '8' ] ]:
        
    option1 = option1_vs_2[0]
    if not option1 in svFitProducerModuleNames:
        continue
    svFitProducerModuleName1 = svFitProducerModuleNames[option1]
    
    option2 = option1_vs_2[1]
    if not option2 in svFitProducerModuleNames:
        continue
    svFitProducerModuleName2 = svFitProducerModuleNames[option2]

    analyzerCorrelationOption1vs2 = cms.EDAnalyzer("SVfitEventHypothesisCorrelationAnalyzer",
        srcEventHypotheses1 = cms.InputTag(svFitProducerModuleName1),
        srcEventHypotheses2 = cms.InputTag(svFitProducerModuleName2),
        srcWeights = cms.VInputTag(),
        dqmDirectory = cms.string("correlation%sVs%s" % (svFitProducerModuleName1, svFitProducerModuleName2)) 
    )
    analyzerNameCorrelationOption1vs2 = "correlation%sVs%s" % (svFitProducerModuleName1, svFitProducerModuleName2)
    setattr(process, analyzerNameCorrelationOption1vs2, analyzerCorrelationOption1vs2)
    process.testSVfitTrackLikelihoodSequence += analyzerCorrelationOption1vs2

    filterTypeOption1betterThan2 = None
    if option1 == 1:
        filterTypeOption1betterThan2 = "BestMatchFilterCandidateToSVfitEventHypothesisByIntegration"
    elif option1 == 5:
        filterTypeOption1betterThan2 = "BestMatchFilterCandidateToSVfitEventHypothesis"
    else:
        raise ValueError("No BestMatchFilter type defined for option = %i" % option1)
    filterOption1betterThan2 = cms.EDFilter(filterTypeOption1betterThan2,
        srcRef = cms.InputTag('genAHdecayToTaus'),
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

#--------------------------------------------------------------------------------
# save plots

process.DQMStore = cms.Service("DQMStore")

process.savePlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string("/data1/veelken/tmp/svFitStudies/testSVfitTrackLikelihoods2_%s_%s_%s_2013Jan03.root" % (sample_type, channel, massPoint))
)

process.q = cms.EndPath(process.savePlots)
#--------------------------------------------------------------------------------

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

processDumpFile = open('testSVfitTrackLikelihoods2.dump' , 'w')
print >> processDumpFile, process.dumpPython()

