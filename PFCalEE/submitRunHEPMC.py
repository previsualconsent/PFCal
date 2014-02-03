#!/usr/bin/env python

import os,sys
import optparse
import commands
import math

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-q', '--queue'      ,    dest='queue'              , help='batch queue'                  , default='2nw')
parser.add_option('-v', '--version'    ,    dest='version'            , help='detector version'             , default=0,      type=int)
parser.add_option('-f', '--datafile'   ,    dest='datafile'           , help='HepMC input file'             , default='example_MyPythia.dat')
parser.add_option('-d', '--datatype'   ,    dest='datatype'           , help='data type'                    , default='PythiaTest')
parser.add_option('-n', '--nevts'      ,    dest='nevts'              , help='number of events to generate' , default=1000,    type=int)
parser.add_option('-o', '--out'        ,    dest='out'                , help='output directory'             , default=os.getcwd() )
parser.add_option('-S', '--no-submit'  ,    action="store_true",  dest='nosubmit'           , help='Do not submit batch job.')
(opt, args) = parser.parse_args()

nevents=opt.nevts

outDir='%s/version_%d/%s/'%(opt.out,opt.version,opt.datatype)
os.system('mkdir -p %s'%outDir)

#wrapper
scriptFile = open('%s/runJob.sh'%(outDir), 'w')
scriptFile.write('#!/bin/bash\n')
scriptFile.write('source %s/g4env.sh\n'%(os.getcwd()))
scriptFile.write('cd %s\n'%(outDir))
scriptFile.write('PFCalEE g4steer.mac %d\n'%opt.version)
scriptFile.write('echo "All done"\n')
scriptFile.close()

#write geant 4 macro
g4Macro = open('%s/g4steer.mac'%(outDir), 'w')
g4Macro.write('/control/verbose 0\n')
g4Macro.write('/control/saveHistory\n')
g4Macro.write('/run/verbose 0\n')
g4Macro.write('/event/verbose 0\n')
g4Macro.write('/tracking/verbose 0\n')

g4Macro.write('/generator/select hepmcAscii\n')
g4Macro.write('/generator/hepmcAscii/open %s/data/%s\n'%(os.getcwd(),opt.datafile))
g4Macro.write('/generator/hepmcAscii/verbose 1\n')

g4Macro.write('/run/beamOn %d\n'%(nevents))
g4Macro.close()

#submit
os.system('chmod u+rwx %s/runJob.sh'%outDir)
if opt.nosubmit : os.system('echo bsub -q %s %s/runJob.sh'%(opt.queue,outDir)) 
else: os.system("bsub -q %s \'%s/runJob.sh\'"%(opt.queue,outDir))

