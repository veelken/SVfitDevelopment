import FWCore.ParameterSet.Config as cms

process = cms.Process("displaySVfitLikelihood")

import os
import re

import TauAnalysis.Configuration.tools.castor as castor
import TauAnalysis.Configuration.tools.eos as eos
from TauAnalysis.Skimming.EventContent_cff import *

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
#process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START53_V15::All')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(11000)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        ##'file:/data1/veelken/CMSSW_5_2_x/skims/genHtautauLeptonPairAcc/user/v/veelken/CMSSW_5_2_x/skims/genTauLeptonsPairAccSkim_ggHiggs125_mutau_5_1_sZL.root'
        ##'/store/user/veelken/CMSSW_5_3_x/skims/SVfitStudies/AccLowPtThresholds/ZplusJets_mutau/genTauLeptonPairSkim_ZplusJets_mutau_1_1_aWJ.root'
        'file:/data1/veelken/CMSSW_5_3_x/skims/selEvents_testSVfitTrackLikelihoods_Z_tautau_90_2013Aug08_AOD.root'                        
    ),
    eventsToProcess = cms.untracked.VEventRange(
        ##'1:32274:12899383',
        ##'1:43800:17505929',
        ##'1:44259:17689319',
        ##'1:102093:40804761',
        ##'1:113221:45252099',
        ##'1:137565:54982242'
        '1:102093:40804761',
        #'1:137071:54784964',
        #'1:154713:61835786',
        #'1:80821:32302660',
        #'1:172248:68844329',
        #'1:3369:1346645',
        #'1:89621:35819785',
        #'1:79157:31638382'
    )
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
            '/data2/veelken/CMSSW_5_3_x/skims/user/veelken/CMSSW_5_3_x/skims/SVfitStudies/AccLowPtThresholds/ZplusJets_mutau_highStats/'
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

##process.source.fileNames = cms.untracked.vstring(inputFileNames)
#--------------------------------------------------------------------------------

process.displaySVfitLikelihoodSequence = cms.Sequence()

#--------------------------------------------------------------------------------
# require event to contain generator level tau lepton pair,
# decaying to the specified decay mode
# matched to tau decay products on generator level

process.load("PhysicsTools/JetMCAlgos/TauGenJets_cfi")
process.displaySVfitLikelihoodSequence += process.tauGenJets
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
process.displaySVfitLikelihoodSequence += process.tauGenJetsSelectorAllHadrons

genElectronsFromTauDecays = None
genMuonsFromTauDecays = None
genTauJetsFromTauDecays = None
genTauPairs = None
genTaus = None
if sample_type == 'Z':
    process.load("TauAnalysis/GenSimTools/gen_decaysFromZs_cfi")
    ##process.genTauJetsFromZs.srcTauGenJets = cms.InputTag('tauGenJetsSelectorAllHadrons')
    process.displaySVfitLikelihoodSequence += process.produceGenDecayProductsFromZs
    genElectronsFromTauDecays = 'genElectronsFromZtautauDecays'
    genMuonsFromTauDecays = 'genMuonsFromZtautauDecays'
    genTauJetsFromTauDecays = 'genHadronsFromZtautauDecays'
    genTauPairs = 'genZdecayToTaus'
    genTaus = 'genTausFromZs'
elif sample_type == 'Higgs':
    process.load("TauAnalysis/GenSimTools/gen_decaysFromAHs_cfi")
    ##process.genTauJetsFromAHs.srcTauGenJets = cms.InputTag('tauGenJetsSelectorAllHadrons')
    process.displaySVfitLikelihoodSequence += process.produceGenDecayProductsFromAHs
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
process.displaySVfitLikelihoodSequence += process.genTauPairMassFilter

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
        src = cms.InputTag('genElectronsFromTauDecays'),
        minNumber = cms.uint32(numElectrons),
        maxNumber = cms.uint32(numElectrons)                      
    ) 
    process.displaySVfitLikelihoodSequence += process.genElectronFilter

if numMuons > 0:    
    process.genMuonFilter = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag('genMuonsFromTauDecays'),
        minNumber = cms.uint32(numMuons),
        maxNumber = cms.uint32(numMuons)                      
    )
    process.displaySVfitLikelihoodSequence += process.genMuonFilter

if numTauJets > 0:
    process.genTauFilter = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag('tauGenJetsSelectorAllHadrons'),
        minNumber = cms.uint32(numTauJets),
        maxNumber = cms.uint32(numTauJets)                      
    )
    process.displaySVfitLikelihoodSequence += process.genTauFilter
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
process.displaySVfitLikelihoodSequence += process.genTauMatchedPFJets

process.load("RecoTauTag/Configuration/RecoPFTauTag_cff")
process.displaySVfitLikelihoodSequence += process.PFTau

process.load("TauAnalysis/Skimming/goldenZmmSelectionVBTFnoMuonIsolation_cfi")
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
process.displaySVfitLikelihoodSequence += process.muonSelectionSequence

process.selectedMuons = cms.EDFilter("MuonAntiOverlapSelector",
    src = cms.InputTag('muons'),
    srcNotToBeFiltered = cms.VInputTag('goodIsoMuons'),
    dRmin = cms.double(0.3),
    invert = cms.bool(True),
    filter = cms.bool(False)                                                          
)
process.selectedMuons.filter = cms.bool(False)
process.displaySVfitLikelihoodSequence += process.selectedMuons

process.genMatchedMuons = cms.EDFilter("MuonAntiOverlapSelector",
    src = cms.InputTag('selectedMuons'),
    srcNotToBeFiltered = cms.VInputTag(genMuonsFromTauDecays),
    dRmin = cms.double(0.3),
    invert = cms.bool(True),
    filter = cms.bool(False)                                                          
)
process.genMatchedMuons.filter = cms.bool(False)
process.displaySVfitLikelihoodSequence += process.genMatchedMuons

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
process.displaySVfitLikelihoodSequence += process.selectedTaus

process.genMatchedTaus = cms.EDFilter("PFTauAntiOverlapSelector",
    src = cms.InputTag('selectedTaus'),
    srcNotToBeFiltered = cms.VInputTag(genTauJetsFromTauDecays),
    dRmin = cms.double(0.3),
    invert = cms.bool(True),
    filter = cms.bool(False)                                                          
)
process.displaySVfitLikelihoodSequence += process.genMatchedTaus
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# require event to contain reconstructed lepton pair,
# matched to tau decay products on generator level

if numElectrons > 0:
    process.recElectronFilter = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag('genMatchedElectrons'), # CV: not implemented yet
        minNumber = cms.uint32(numElectrons),
        maxNumber = cms.uint32(numElectrons)                      
    ) 
    process.displaySVfitLikelihoodSequence += process.recElectronFilter

if numMuons > 0:    
    process.recMuonFilter = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag('genMatchedMuons'),
        minNumber = cms.uint32(numMuons),
        maxNumber = cms.uint32(numMuons)                      
    )
    process.displaySVfitLikelihoodSequence += process.recMuonFilter

if numTauJets > 0:
    process.recTauFilter = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag('genMatchedTaus'),
        minNumber = cms.uint32(numTauJets),
        maxNumber = cms.uint32(numTauJets)                      
    )
    process.displaySVfitLikelihoodSequence += process.recTauFilter
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# produce genMET
process.load("RecoMET.Configuration.GenMETParticles_cff")
process.displaySVfitLikelihoodSequence += process.genParticlesForMETAllVisible

process.load("RecoMET.METProducers.genMetTrue_cfi")
process.genMetFromGenParticles = process.genMetTrue.clone(
    src = cms.InputTag('genParticlesForMETAllVisible'),
    alias = cms.string('genMetFromGenParticles')
)
process.displaySVfitLikelihoodSequence += process.genMetFromGenParticles
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# reconstructed Type 1 corrected PFMET and its uncertainty 

process.load("JetMETCorrections/Type1MET/pfMETCorrectionType0_cfi")
process.displaySVfitLikelihoodSequence += process.type0PFMEtCorrection

process.load("JetMETCorrections/Type1MET/pfMETsysShiftCorrections_cfi")
process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runAplusBvsNvtx_mc
process.displaySVfitLikelihoodSequence += process.pfMEtSysShiftCorrSequence

process.load("JetMETCorrections/Type1MET/pfMETCorrections_cff")
process.pfType1CorrectedMet.applyType0Corrections = cms.bool(False)
process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag(
    cms.InputTag('pfJetMETcorr', 'type1'),
    cms.InputTag('pfMETcorrType0'),
    cms.InputTag('pfMEtSysShiftCorr')
)
process.displaySVfitLikelihoodSequence += process.producePFMETCorrections

# CV: compute PFMET significance cov. matrix for uncorrected jets
#     in order to include pile-up jets
#    (which to a significant fraction get killed by L1Fastjet corrections)
process.ak5PFJetsNotOverlappingWithLeptons = cms.EDFilter("PFJetAntiOverlapSelector",
    src = cms.InputTag('ak5PFJets'),
    srcNotToBeFiltered = cms.VInputTag(
        ##'genMatchedElectrons', # CV: not implemented yet
        'genMatchedMuons',
        'genMatchedTaus'
    ),
    dRmin = cms.double(0.5),
    invert = cms.bool(False),
    filter = cms.bool(False)                                                          
)
process.displaySVfitLikelihoodSequence += process.ak5PFJetsNotOverlappingWithLeptons

from JetMETCorrections.Type1MET.pfMETCorrections_cff import pfCandsNotInJet
process.pfCandsNotInJetForPFMEtSignCovMatrix = pfCandsNotInJet.clone()
process.displaySVfitLikelihoodSequence += process.pfCandsNotInJetForPFMEtSignCovMatrix

from RecoMET.METProducers.METSigParams_cfi import *
process.pfMEtSignCovMatrix = cms.EDProducer("PFMEtSignCovMatrixProducer",
    METSignificance_params,                     
    src = cms.VInputTag(
        ##'genMatchedElectrons', # CV: not implemented yet
        'genMatchedMuons',
        'genMatchedTaus',
        'ak5PFJetsNotOverlappingWithLeptons',                                        
        'pfCandsNotInJetForPFMEtSignCovMatrix'
    ),
    addJERcorr = cms.PSet(
        inputFileName = cms.FileInPath('PhysicsTools/PatUtils/data/pfJetResolutionMCtoDataCorrLUT.root'),
        lutName = cms.string('pfJetResolutionMCtoDataCorrLUT')
    )
)
process.displaySVfitLikelihoodSequence += process.pfMEtSignCovMatrix
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# rerun primary event vertex reconstruction
process.load("RecoVertex/PrimaryVertexProducer/OfflinePrimaryVertices_cfi")
process.offlinePrimaryVertices.verbose = cms.untracked.bool(True)
##process.displaySVfitLikelihoodSequence += process.offlinePrimaryVertices
process.load("RecoVertex/PrimaryVertexProducer/OfflinePrimaryVerticesWithBS_cfi")
##process.displaySVfitLikelihoodSequence += process.offlinePrimaryVerticesWithBS
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# select primary event vertex
process.load("TauAnalysis/RecoTools/recoVertexSelectionByLeptonTracks_cff")
process.selectedPrimaryVertexQuality.src = cms.InputTag('offlinePrimaryVerticesWithBS')
process.selectedPrimaryVertexByLeptonMatch.srcLeptons = cms.VInputTag(
    ##'genMatchedElectrons', # CV: not implemented yet
    'genMatchedMuons',
    'genMatchedTaus'
)
process.displaySVfitLikelihoodSequence += process.selectPrimaryVertexByLeptonTracks

# require event to have exactly one vertex associated to tracks of tau decay products
process.recEventVertexFilter = cms.EDFilter("VertexCountFilter",
    src = cms.InputTag('selectedPrimaryVertexByLeptonMatch'),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(1)                                            
)
process.displaySVfitLikelihoodSequence += process.recEventVertexFilter
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# print list of generator level particles in the event
##
##process.printGenParticleList = cms.EDAnalyzer("ParticleListDrawer",
##    src = cms.InputTag('genParticles'),
##    maxEventsToPrint = cms.untracked.int32(100)
##)
##process.displaySVfitLikelihoodSequence += process.printGenParticleList
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# fill histograms

process.displaySVfitLikelihood = cms.EDAnalyzer("SVfitLikelihoodDisplay",
    srcGenParticles = cms.InputTag('genParticles'),
    srcElectrons = cms.InputTag(''),                                            
    ##srcMuons = cms.InputTag('genMatchedMuons'),
    srcMuons = cms.InputTag(''), 
    srcTaus = cms.InputTag('genMatchedTaus'),
    srcMEt = cms.InputTag('pfType1CorrectedMet'),
    ##srcMEt = cms.InputTag('genMetFromGenParticles'),                                            
    srcMEtCov = cms.InputTag('pfMEtSignCovMatrix'),                                            
    srcVertices = cms.InputTag('selectedPrimaryVertexByLeptonMatch'),
    srcBeamSpot = cms.InputTag('offlineBeamSpot'),
    applyBeamSpotConstraint = cms.bool(True),
    sfProdVertexCov = cms.double(1.0),
    sfDecayVertexCov = cms.double(1.0),
    makePlots_vs_X = cms.bool(True),
    makePlots_vs_GJangle = cms.bool(False),
    srcWeights = cms.VInputTag(),
    verbosity = cms.int32(0)                                                
)
process.displaySVfitLikelihoodSequence += process.displaySVfitLikelihood
#--------------------------------------------------------------------------------

process.p = cms.Path(process.displaySVfitLikelihoodSequence)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

processDumpFile = open('displaySVfitLikelihood.dump' , 'w')
print >> processDumpFile, process.dumpPython()
