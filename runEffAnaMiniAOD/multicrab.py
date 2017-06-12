#!/usr/bin/env python
"""
This is a small script that does the equivalent of multicrab.
"""
import os
from optparse import OptionParser

from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException


def getOptions():
    """
    Parse and return the arguments provided by the user.
    """
    usage = ("Usage: %prog --crabCmd CMD [--workArea WAD --crabCmdOpts OPTS]"
             "\nThe multicrab command executes 'crab CMD OPTS' for each project directory contained in WAD"
             "\nUse multicrab -h for help")

    parser = OptionParser(usage=usage)

    parser.add_option('-c', '--crabCmd',
                      dest = 'crabCmd',
                      default = '',
                      help = "crab command",
                      metavar = 'CMD')

    parser.add_option('-w', '--workArea',
                      dest = 'workArea',
                      default = '',
                      help = "work area directory (only if CMD != 'submit')",
                      metavar = 'WAD')

    parser.add_option('-o', '--crabCmdOpts',
                      dest = 'crabCmdOpts',
                      default = '',
                      help = "options for crab command CMD",
                      metavar = 'OPTS')

    (options, arguments) = parser.parse_args()

    if arguments:
        parser.error("Found positional argument(s): %s." % (arguments))
    if not options.crabCmd:
        parser.error("(-c CMD, --crabCmd=CMD) option not provided.")
    if options.crabCmd != 'submit':
        if not options.workArea:
            parser.error("(-w WAR, --workArea=WAR) option not provided.")
        if not os.path.isdir(options.workArea):
            parser.error("'%s' is not a valid directory." % (options.workArea))

    return options


def main():

    options = getOptions()

    # The submit command needs special treatment.
    if options.crabCmd == 'submit':

        #--------------------------------------------------------
        # This is the base config:
        #--------------------------------------------------------
        from CRABClient.UserUtilities import config
        config = config()

        config.General.requestName = None
        config.General.workArea = "Full2016"

        config.JobType.pluginName = 'Analysis'
        config.JobType.psetName = "runEffAnaMiniAOD_cfg.py"

        config.Data.inputDataset = None
        config.Data.splitting = 'LumiBased'
        config.Data.unitsPerJob = 10
        # config.Data.totalUnits = 30
        config.Data.lumiMask = "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"
        config.Data.outLFNDirBase = "/store/group/cafcern/pyu/EleEfficiencyHWW"

        config.Site.storageSite = "T2_CH_CERN"
        #--------------------------------------------------------

        # Will submit one task for each of these input datasets.
        inputDatasets = [
                         # '/SingleElectron/Run2016B-03Feb2017_ver1-v1/MINIAOD',
                         # '/SingleElectron/Run2016B-03Feb2017_ver2-v2/MINIAOD',
                         # '/SingleElectron/Run2016C-03Feb2017-v1/MINIAOD',
                         # '/SingleElectron/Run2016D-03Feb2017-v1/MINIAOD',
                         # '/SingleElectron/Run2016E-03Feb2017-v1/MINIAOD',
                         # '/SingleElectron/Run2016F-03Feb2017-v1/MINIAOD',
                         '/SingleElectron/Run2016G-03Feb2017-v1/MINIAOD',
                         '/SingleElectron/Run2016H-03Feb2017_ver2-v1/MINIAOD',
                         '/SingleElectron/Run2016H-03Feb2017_ver3-v1/MINIAOD',
                        ]

        for inDS in inputDatasets:
            # inDS is of the form /A/B/C. Since B is unique for each inDS, use this in the CRAB request name.
            config.General.requestName = inDS.split('/')[2]
            config.Data.inputDataset = inDS
            config.Data.outputDatasetTag = '%s_%s' % (config.General.workArea, config.General.requestName)
            # Submit.
            try:
                print "Submitting for input dataset %s" % (inDS)
                crabCommand(options.crabCmd, config = config, *options.crabCmdOpts.split())
            except HTTPException as hte:
                print "Submission for input dataset %s failed: %s" % (inDS, hte.headers)
            except ClientException as cle:
                print "Submission for input dataset %s failed: %s" % (inDS, cle)

    # All other commands can be simply executed.
    elif options.workArea:

        for dir in os.listdir(options.workArea):
            projDir = os.path.join(options.workArea, dir)
            if not os.path.isdir(projDir):
                continue
            # Execute the crab command.
            msg = "Executing (the equivalent of): crab %s --dir %s %s" % (options.crabCmd, projDir, options.crabCmdOpts)
            print "-"*len(msg)
            print msg
            print "-"*len(msg)
            try:
                crabCommand(options.crabCmd, dir = projDir, *options.crabCmdOpts.split())
            except HTTPException as hte:
                print "Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, hte.headers)
            except ClientException as cle:
                print "Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, cle)


if __name__ == '__main__':
    main()
