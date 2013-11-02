#!/usr/bin/env python

import copy
import os
import shlex
import socket
import subprocess
import sys
import time

from TauAnalysis.CandidateTools.tools.buildConfigFilesSVfitPerformanceAnalysis import *
from TauAnalysis.TauIdEfficiency.tools.buildConfigFilesTauIdEffAnalysis import buildConfigFile_hadd
from TauAnalysis.Configuration.tools.jobtools import make_bsub_script
from TauAnalysis.Configuration.tools.harvestingLXBatch import make_harvest_scripts
from TauAnalysis.Configuration.tools.harvesting import castor_source
import TauAnalysis.Configuration.tools.castor as castor

version = '2012Mar13'

inputFilePath  = '/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/SVfitStudies/'
harvestingFilePath = '/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/harvesting/SVfitStudies/AHtautau_2012May01/'
outputFilePath = '/tmp/veelken/svFitStudies/'
#outputFilePath = '/data1/veelken/tmp/svFitStudies/AHtautau/2012Apr19/'

samplesToAnalyze = [
    'ggHiggs120',
    'vbfHiggs120',
    'ggHiggs130',
    'vbfHiggs130',
    'ggPhi160',
    'bbPhi160',
    'ggPhi200',
    'bbPhi200',
    'ggPhi300',
    'bbPhi300',
    'ggPhi450',
    'bbPhi450',
    'ZplusJets'
]

channelsToAnalyze = [
    ##'muTau',
    ##'eleTau',
    'eleMu',
    ##'diTau'
]

metResolutions = [
    None,
    5.,
    10.,
    15.,
    20.,
    25.
]    

#runSVfitEventHypothesisAnalyzer = True
runSVfitEventHypothesisAnalyzer = False
runLXBatchHarvesting = True
#runLXBatchHarvesting = False
#runMakePlots = True
runMakePlots = False

if runSVfitEventHypothesisAnalyzer:
    print "submitting analysis jobs..."
elif runLXBatchHarvesting:
    print "submitting harvesting jobs..."
elif runMakePlots:
    print "making final plots..."

maxEventsPerJob = 100 # max. number of events analyzed per job

def createFilePath_recursively(filePath):
    filePath_items = filePath.split('/')
    currentFilePath = "/"
    for filePath_item in filePath_items:
        currentFilePath = os.path.join(currentFilePath, filePath_item)
        if len(currentFilePath) <= 1:
            continue
        if not os.path.exists(currentFilePath):
            #sys.stdout.write("creating directory %s\n" % currentFilePath)
            os.mkdir(currentFilePath)

hostname = socket.gethostname()
print "hostname = %s" % hostname
if hostname == 'ucdavis.cern.ch':
    print "Running on %s" % hostname

try:
    castor.rfstat(harvestingFilePath)
except RuntimeError:
    # harvestingFilePath does not yet exist, create it
    print "harvestingFilePath does not yet exist, creating it."
    os.system("rfmkdir %s" % harvestingFilePath)
    os.system("rfchmod 777 %s" % harvestingFilePath)
    
outputFilePath = os.path.join(outputFilePath, version)
print "outputFilePath = %s" % outputFilePath
createFilePath_recursively(outputFilePath)
for channelToAnalyze in channelsToAnalyze:
    outputFilePath_channel = os.path.join(outputFilePath, channelToAnalyze)
    createFilePath_recursively(outputFilePath_channel)

configFilePath = os.path.join(os.getcwd(), "lxbatch")
print "configFilePath = %s" % configFilePath
createFilePath_recursively(configFilePath)
for channelToAnalyze in channelsToAnalyze:
    configFilePath_channel = os.path.join(configFilePath, channelToAnalyze)
    createFilePath_recursively(configFilePath_channel)

logFilePath = os.path.join(os.getcwd(), "lxbatch_log")
print "logFilePath = %s" % logFilePath
createFilePath_recursively(logFilePath)

execDir = "%s/bin/%s/" % (os.environ['CMSSW_BASE'], os.environ['SCRAM_ARCH'])

executable_cmsRun = 'cmsRun'
executable_bsub = 'bsub'
executable_waitForLXBatchJobs = 'python %s/src/TauAnalysis/Configuration/python/tools/waitForLXBatchJobs.py' % os.environ['CMSSW_BASE']
executable_rfcp = 'rfcp'
executable_rfrm = 'rfrm'
executable_stager = 'stager_get'
executable_hadd = 'hadd -f'
executable_shell = '/bin/csh'

bsubQueue = "1nd"

configFileName_SVfitEventHypothesisAnalyzer_template = "runSVfitPerformanceAnalysis_AHtautau_cfg.py"

def runCommand(commandLine):
    sys.stdout.write("%s\n" % commandLine)
    args = shlex.split(commandLine)
    retVal = subprocess.Popen(args, stdout = subprocess.PIPE)
    retVal.wait()
    return retVal

# find and delete "bad" files
##files = [ file_info for file_info in castor.nslsl(harvestingFilePath) ]
##for file in files:
##    if file['size'] < 1000:
##        runCommand("%s %s" % (executable_rfrm, file['path']))

def getMEtResolution_label(metResolution):
    retVal = "MEtResMC"
    if metResolution is not None:
        retVal = "MEtRes%1.0f" % metResolution
        retVal = retVal.replace(".", "_")
    return retVal

# CV: fill mapping of fileName to number of events contained in file into temporary cache
#     in order to reduce castor file I/O
print "initializing mapping of fileNames to number of events contained in each file..."
numEventsMap = {}
fileNamesToMap = []
for sampleToAnalyze in samplesToAnalyze:
    for channelToAnalyze in channelsToAnalyze:
        inputFilePath_channel = os.path.join(inputFilePath, version, channelToAnalyze)
        inputFileNames = None
        if inputFilePath_channel.find('/castor/') != -1:
            inputFileNames = [ file_info['path'] for file_info in castor.nslsl(inputFilePath_channel) ]
        else:
            inputFileNames = os.listdir(inputFilePath_channel)
        for inputFileName in inputFileNames:        
            if inputFileName.find("".join(['_', sampleToAnalyze, '_'])) != -1 or \
               inputFileName.find("".join(['/', sampleToAnalyze, '_'])) != -1:
                fileNamesToMap.append(inputFileName)
                # CV: request inputFiles located on castor to be prestaged
                #     in order to speed-up computation of numbers of events contained in each file
                #     by 'buildConfigFile_SVfitEventHypothesisAnalyzer' function later
                if inputFilePath_channel.find('/castor/') != -1:
                    commandLine = '%s -M %s -U myfiles' % (executable_stager, inputFileName)
                    runCommand(commandLine)
print " done."

#--------------------------------------------------------------------------------
#
# build config files for running FWLiteTauFakeRateAnalyzer macro on lxbatch
#
fileNames_SVfitEventHypothesisAnalyzer         = {}
bsubJobNames_SVfitEventHypothesisAnalyzer      = {}
bjobListFileNames_SVfitEventHypothesisAnalyzer = {}
for sampleToAnalyze in samplesToAnalyze:
    fileNames_SVfitEventHypothesisAnalyzer[sampleToAnalyze]         = {}
    bsubJobNames_SVfitEventHypothesisAnalyzer[sampleToAnalyze]      = {}
    bjobListFileNames_SVfitEventHypothesisAnalyzer[sampleToAnalyze] = {}
    for channelToAnalyze in channelsToAnalyze:
        fileNames_SVfitEventHypothesisAnalyzer[sampleToAnalyze][channelToAnalyze]         = {}
        bsubJobNames_SVfitEventHypothesisAnalyzer[sampleToAnalyze][channelToAnalyze]      = {}
        bjobListFileNames_SVfitEventHypothesisAnalyzer[sampleToAnalyze][channelToAnalyze] = {}
        configFilePath_channel = os.path.join(configFilePath, channelToAnalyze)
        for metResolution in metResolutions:
            metResolution_label = getMEtResolution_label(metResolution)
            print "building config file for sample = %s, channel = %s, metResolution = %s" % \
              (sampleToAnalyze, channelToAnalyze, metResolution_label)
            retVal_SVfitEventHypothesisAnalyzer = \
              buildConfigFile_SVfitEventHypothesisAnalyzer(sampleToAnalyze, channelToAnalyze, metResolution,
                                                           configFileName_SVfitEventHypothesisAnalyzer_template,
                                                           os.path.join(inputFilePath, version, channelToAnalyze), 1, maxEventsPerJob,
                                                           configFilePath_channel, logFilePath, harvestingFilePath, numEventsMap)        
            if retVal_SVfitEventHypothesisAnalyzer is not None:
                fileNames_SVfitEventHypothesisAnalyzer[sampleToAnalyze][channelToAnalyze][metResolution_label] = \
                  retVal_SVfitEventHypothesisAnalyzer
            else:
                fileNames_SVfitEventHypothesisAnalyzer[sampleToAnalyze][channelToAnalyze][metResolution_label] = {
                    'inputFileNames'  : [],
                    'configFileNames' : [],
                    'outputFileNames' : [],
                    'logFileNames'    : []
                }
            fileNameEntry = fileNames_SVfitEventHypothesisAnalyzer[sampleToAnalyze][channelToAnalyze][metResolution_label]
            fileNameEntry['bsubScriptFileNames'] = []
  
            bsubJobNames_SVfitEventHypothesisAnalyzer[sampleToAnalyze][channelToAnalyze][metResolution_label] = []
            bjobListFileNames_SVfitEventHypothesisAnalyzer[sampleToAnalyze][channelToAnalyze][metResolution_label] = None

            if retVal_SVfitEventHypothesisAnalyzer is None:
                continue

            for i in range(len(retVal_SVfitEventHypothesisAnalyzer['inputFileNames'])):

                # The None in the tuple indicates that batch job has no dependencies on other batch jobs
                input_files_and_jobs = \
                  [ (None, os.path.join(inputFilePath, version, channelToAnalyze, inputFileName)) \
                    for inputFileName in retVal_SVfitEventHypothesisAnalyzer['inputFileNames'][i] ]

                def log_file_maker(job_hash):
                    log_fileName = os.path.join(logFilePath, retVal_SVfitEventHypothesisAnalyzer['logFileNames'][i])
                    # CV: delete log-files from previous job submissions
                    os.system("rm -f %s" % log_fileName)
                    return log_fileName

                # Build script for batch job submission
                jobName, bsubScript = make_bsub_script(
                    os.path.join(harvestingFilePath, retVal_SVfitEventHypothesisAnalyzer['outputFileNames'][i]),
                    input_files_and_jobs,
                    log_file_maker,
                    "%s %s" % (executable_cmsRun,
                               os.path.join(configFilePath_channel, retVal_SVfitEventHypothesisAnalyzer['configFileNames'][i])))

                #print "configFilePath_channel = %s" % configFilePath_channel
                #print "retVal_SVfitEventHypothesisAnalyzer['logFileNames'][i] = %s" % \
                #  retVal_SVfitEventHypothesisAnalyzer['logFileNames'][i]
            
                bsubScriptFileName = \
                  os.path.join(configFilePath_channel, retVal_SVfitEventHypothesisAnalyzer['logFileNames'][i].replace(".log", ".sh"))
                bsubScriptFile = open(bsubScriptFileName, "w")
                bsubScriptFile.write(bsubScript)
                bsubScriptFile.close()
                time.sleep(0.100) # CV: wait for 100 milliseconds in order to avoid opening/closing too many files in too short time

                fileNameEntry['bsubScriptFileNames'].append(bsubScriptFileName)
                
                bsubJobName = "svFitPerfAna%s%s%s_%i" % (sampleToAnalyze, channelToAnalyze, metResolution_label, i)
                bsubJobNames_SVfitEventHypothesisAnalyzer[sampleToAnalyze][channelToAnalyze][metResolution_label].append(bsubJobName)

                bjobListFileName = \
                  os.path.join(configFilePath_channel, "batchJobs_SVfitEventHypothesisAnalyzer_%s_%s_pf%s.lst" %
                    (sampleToAnalyze, channelToAnalyze, metResolution_label))
                bjobListFile = open(bjobListFileName, "w")
                for bsubJobName in bsubJobNames_SVfitEventHypothesisAnalyzer[sampleToAnalyze][channelToAnalyze][metResolution_label]:
                    bjobListFile.write("%s\n" % bsubJobName)
                bjobListFile.close()
        
                bjobListFileNames_SVfitEventHypothesisAnalyzer[sampleToAnalyze][channelToAnalyze][metResolution_label] = bjobListFileName
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# build shell script for running 'hadd' in order to collect histograms
# for single sample and single channel into single .root file
#
bsubFileNames_harvesting    = {}
bsubJobNames_harvesting     = {}
bsubJobNames_harvesting_all = []
for sampleToAnalyze in samplesToAnalyze:
    bsubFileNames_harvesting[sampleToAnalyze] = {}
    bsubJobNames_harvesting[sampleToAnalyze]  = {}
    for channelToAnalyze in channelsToAnalyze:
        bsubFileNames_harvesting[sampleToAnalyze][channelToAnalyze] = {}
        bsubJobNames_harvesting[sampleToAnalyze][channelToAnalyze]  = {}
        configFilePath_channel = os.path.join(configFilePath, channelToAnalyze)
        outputFilePath_channel = os.path.join(outputFilePath, channelToAnalyze)
        for metResolution in metResolutions:
            metResolution_label = getMEtResolution_label(metResolution)
    
            plot_regex = r"[a-zA-Z0-9._]+"
            skim_regex = r"dont match anything"
        
            def local_copy_mapper(sample):
                return os.path.join(
                    outputFilePath_channel,
                    'svFitPerformanceAnalysisPlots_%s_%s_%s_%s_harvested.root' %
                      (channelToAnalyze, sampleToAnalyze, metResolution_label, version))

            fileNameEntry = fileNames_SVfitEventHypothesisAnalyzer[sampleToAnalyze][channelToAnalyze][metResolution_label]

            inputFileInfos = []        
            for inputFileName in fileNameEntry['outputFileNames']:
                inputFileInfo = {
                    'path'        : os.path.join(harvestingFilePath, inputFileName),
                    'size'        : 1,           # dummy
                    'time'        : time.localtime(),
                    'file'        : inputFileName,
                    'permissions' : 'mrw-r--r--' # "ordinary" file access permissions
                }
                #print "inputFileInfo = %s" % inputFileInfo
                inputFileInfos.append(inputFileInfo)

            retVal_make_harvest_scripts = make_harvest_scripts(
                plot_regex,
                skim_regex,
                channel = channelToAnalyze,
                sampleToAnalyze = sampleToAnalyze,
                job_id = "_".join([ metResolution_label, version ]),
                input_files_info = inputFileInfos,
                harvester_command = executable_hadd,
                abort_on_rfcp_error = False,
                castor_output_directory = harvestingFilePath,
                script_directory = configFilePath_channel,
                merge_script_name = \
                os.path.join(configFilePath_channel, "_".join([ 'submit', sampleToAnalyze, channelToAnalyze,
                                                                metResolution_label, 'merge' ]) + '.sh'),
                local_copy_mapper = local_copy_mapper,
                chunk_size = 2.e+9, # 2 GB
                check_old_files = False,
                max_bsub_concurrent_file_access = 250,
                verbosity = 0
            )

            bsubFileNames_harvesting[sampleToAnalyze][channelToAnalyze][metResolution_label] = retVal_make_harvest_scripts

            bsubJobName = "harvest%s%s%s" % (sampleToAnalyze, channelToAnalyze, metResolution_label)
            bsubJobNames_harvesting[sampleToAnalyze][channelToAnalyze][metResolution_label] = bsubJobName

            if len(retVal_make_harvest_scripts['final_harvest_files']) > 0:
                bsubJobNames_harvesting_all.append(bsubJobName)

bjobListFileNames_harvesting = {}
for channelToAnalyze in channelsToAnalyze:
    configFilePath_channel = os.path.join(configFilePath, channelToAnalyze)
    bjobListFileNames_harvesting[channelToAnalyze] = \
      os.path.join(configFilePath_channel, "batchJobs_harvesting_%s.lst" % channelToAnalyze)
    bjobListFile_harvesting = open(bjobListFileNames_harvesting[channelToAnalyze], "w")
    for sampleToAnalyze in samplesToAnalyze:
        for metResolution in metResolutions:
            metResolution_label = getMEtResolution_label(metResolution)
            for bsubJobName in bsubFileNames_harvesting[sampleToAnalyze][channelToAnalyze][metResolution_label]['bsub_job_names']:        
                bjobListFile_harvesting.write("%s\n" % bsubJobName)
    bjobListFile_harvesting.close()
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# build shell script for running 'hadd' in order to collect histograms
# for all samples and event selections into single .root file
#
final_haddInputFileNames  = {}
final_haddShellFileNames  = {}
final_haddOutputFileNames = {}
final_haddJobNames        = {}
final_haddLogFileNames    = {}
for channelToAnalyze in channelsToAnalyze:
    final_haddInputFileNames[channelToAnalyze] = []
    outputFilePath_channel = os.path.join(outputFilePath, channelToAnalyze)
    for sampleToAnalyze in samplesToAnalyze:    
        for metResolution in metResolutions:
            metResolution_label = getMEtResolution_label(metResolution)
            for final_harvest_file in \
              bsubFileNames_harvesting[sampleToAnalyze][channelToAnalyze][metResolution_label]['final_harvest_files']:
                # CV:
                #    (1) file name of final harvesting output file is stored at index[1] in final_harvest_file-tuple
                #       (cf. TauAnalysis/Configuration/python/tools/harvestingLXBatch.py)
                #    (2) assume that .root files containing histograms for single sample and single channel
                #        are copied to local disk via rfcp prior to running 'hadd'
                final_haddInputFileNames[channelToAnalyze].append(
                  os.path.join(outputFilePath_channel, os.path.basename(final_harvest_file[1])))
    final_haddShellFileNames[channelToAnalyze] = \
      os.path.join(configFilePath, 'harvestSVfitPerformanceHistograms_%s.csh' % version)
    final_haddOutputFileNames[channelToAnalyze] = \
      os.path.join(outputFilePath_channel, 'svFitPerformanceAnalysisPlots_all_%s.root' % version)
    retVal_final_hadd = \
      buildConfigFile_hadd(executable_hadd,
                           final_haddShellFileNames[channelToAnalyze],
                           final_haddInputFileNames[channelToAnalyze],
                           final_haddOutputFileNames[channelToAnalyze])
    final_haddJobNames[channelToAnalyze] = 'final_hadd_%s' % channelToAnalyze
    final_haddLogFileNames[channelToAnalyze] = retVal_final_hadd['logFileName']
#--------------------------------------------------------------------------------

def make_MakeFile_vstring(list_of_strings):
    retVal = ""
    for i, string_i in enumerate(list_of_strings):
        if i > 0:
            retVal += " "
        retVal += string_i
    return retVal

# done building config files, now build Makefile...
makeFileName = "Makefile_SVfitPerformanceAnalysis_AHtautau"
makeFileName = makeFileName + "_" + "_".join(channelsToAnalyze)
makeFile = open(makeFileName, "w")
makeFile.write("\n")
outputFileNames_make = []
if runSVfitEventHypothesisAnalyzer:
    for sampleToAnalyze in samplesToAnalyze:
        for channelToAnalyze in channelsToAnalyze:
            for metResolution in metResolutions:
                metResolution_label = getMEtResolution_label(metResolution)
                fileNameEntry = fileNames_SVfitEventHypothesisAnalyzer[sampleToAnalyze][channelToAnalyze][metResolution_label]
                if fileNameEntry is None or len(fileNameEntry['inputFileNames']) == 0:
                    continue
                for i in range(len(fileNameEntry['inputFileNames'])):
                    outputFileNames_make.append(fileNameEntry['outputFileNames'][i])
if runLXBatchHarvesting:
    for sampleToAnalyze in samplesToAnalyze:
        for channelToAnalyze in channelsToAnalyze:
            for metResolution in metResolutions:
                metResolution_label = getMEtResolution_label(metResolution)
                outputFileNames_make.append(bsubJobNames_harvesting[sampleToAnalyze][channelToAnalyze][metResolution_label])
if runMakePlots:
    for channelToAnalyze in channelsToAnalyze:
        outputFileNames_make.append(final_haddJobNames[channelToAnalyze])
makeFile.write("all: %s\n" %
  (make_MakeFile_vstring(outputFileNames_make)))
makeFile.write("\techo 'Finished running SVfitPerformanceAnalysis.'\n")
makeFile.write("\n")
if runSVfitEventHypothesisAnalyzer:
    for sampleToAnalyze in samplesToAnalyze:
        for channelToAnalyze in channelsToAnalyze:
            for metResolution in metResolutions:
                metResolution_label = getMEtResolution_label(metResolution)
                fileNameEntry = fileNames_SVfitEventHypothesisAnalyzer[sampleToAnalyze][channelToAnalyze][metResolution_label]
                if fileNameEntry is None or len(fileNameEntry['inputFileNames']) == 0:
                    continue
                bsubJobEntry = bsubJobNames_SVfitEventHypothesisAnalyzer[sampleToAnalyze][channelToAnalyze][metResolution_label]
                for i in range(len(fileNameEntry['inputFileNames'])):
                    makeFile.write("%s:\n" %
                      (fileNameEntry['outputFileNames'][i]))
                    makeFile.write("\t%s -q %s -J %s < %s\n" %
                      (executable_bsub,
                       bsubQueue,
                       bsubJobEntry[i],
                       fileNameEntry['bsubScriptFileNames'][i]))
makeFile.write("\n")
if runLXBatchHarvesting:     
    for sampleToAnalyze in samplesToAnalyze:
        for channelToAnalyze in channelsToAnalyze:
            for metResolution in metResolutions:
                metResolution_label = getMEtResolution_label(metResolution)
                fileNameEntry = fileNames_SVfitEventHypothesisAnalyzer[sampleToAnalyze][channelToAnalyze][metResolution_label]
                if runSVfitEventHypothesisAnalyzer:
                    makeFile.write("%s: %s\n" %
                      (bsubJobNames_harvesting[sampleToAnalyze][channelToAnalyze],
                       make_MakeFile_vstring(fileNameEntry['outputFileNames'])))
                    makeFile.write("\t%s %s\n" %
                      (executable_waitForLXBatchJobs,
                       bjobListFileNames_SVfitEventHypothesisAnalyzer[sampleToAnalyze][channelToAnalyze][metResolution_label]))
                else:
                    makeFile.write("%s:\n" %
                      (bsubJobNames_harvesting[sampleToAnalyze][channelToAnalyze][metResolution_label]))
                makeFile.write("\t%s %s\n" %
                  (executable_shell,
                   bsubFileNames_harvesting[sampleToAnalyze][channelToAnalyze][metResolution_label]['harvest_script_name']))
if runMakePlots:
    for channelToAnalyze in channelsToAnalyze:
        makeFile.write("%s:\n" %
          (final_haddJobNames[channelToAnalyze]))
        if runLXBatchHarvesting:  
            makeFile.write("\t%s %s\n" %
              (executable_waitForLXBatchJobs,
               bjobListFileNames_harvesting[channelToAnalyze]))
        for final_haddInputFileName in final_haddInputFileNames[channelToAnalyze]:
            makeFile.write("\t%s %s %s\n" %
             (executable_rfcp,
               os.path.join(harvestingFilePath, os.path.basename(final_haddInputFileName)),
               final_haddInputFileName))
        makeFile.write("\t%s %s\n" %
          (executable_shell,
           final_haddShellFileNames[channelToAnalyze]))
makeFile.write("\n")
makeFile.write(".PHONY: clean\n")
makeFile.write("clean:\n")
for sampleToAnalyze in samplesToAnalyze:
    for channelToAnalyze in channelsToAnalyze:
        for metResolution in metResolutions:
            metResolution_label = getMEtResolution_label(metResolution)
            fileNameEntry = fileNames_SVfitEventHypothesisAnalyzer[sampleToAnalyze][channelToAnalyze][metResolution_label]
            if runSVfitEventHypothesisAnalyzer:
                if fileNameEntry is None:
                    continue
                for outputFileName in fileNameEntry['outputFileNames']:
                    makeFile.write("\t- %s %s\n" %
                      (executable_rfrm,
                       os.path.join(harvestingFilePath, outputFileName)))
            if runLXBatchHarvesting: 
                for final_harvest_file in bsubFileNames_harvesting[sampleToAnalyze][channelToAnalyze][metResolution_label]['final_harvest_files']:
                    # CV: file name of final harvesting output file is stored at index[1] in final_harvest_file-tuple
                    #    (cf. TauAnalysis/Configuration/python/tools/harvestingLXBatch.py)    
                    makeFile.write("\t- %s %s\n" %
                      (executable_rfrm,
                       os.path.join(harvestingFilePath, final_harvest_file[1])))
for channelToAnalyze in channelsToAnalyze:                    
    makeFile.write("\trm -f %s\n" % make_MakeFile_vstring(final_haddInputFileNames[channelToAnalyze]))            
    makeFile.write("\trm -f %s\n" % final_haddShellFileNames[channelToAnalyze])
    makeFile.write("\trm -f %s\n" % final_haddOutputFileNames[channelToAnalyze])
makeFile.write("\techo 'Finished deleting old files.'\n")
makeFile.write("\n")
makeFile.close()

print("Finished building Makefile. Now execute 'make -j 8 -f %s'." % makeFileName)
