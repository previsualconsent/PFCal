#!/usr/bin/env python

import os,sys
import optparse
import commands
import math

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-q', '--queue'      ,    dest='queue'              , help='batch queue'                  , default='2nw')
parser.add_option('-v', '--version'    ,    dest='version'            , help='detector version'             , default=0,      type=int)
parser.add_option('-a', '--alpha'      ,    dest='alpha'              , help='incidence angle'              , default=0,      type=float)
parser.add_option('-g', '--gun'        ,    dest='gun'                , help='particle to shoot'            , default='e-')
parser.add_option('-n', '--nevts'      ,    dest='nevts'              , help='number of events to generate' , default=1000,    type=int)
parser.add_option('-o', '--out'        ,    dest='out'                , help='output directory'             , default=os.getcwd() )
parser.add_option('-e', '--eos'        ,    dest='eos'                , help='save root file to EOS',         default='')
parser.add_option('-S', '--no-submit'  ,    action="store_true",  dest='nosubmit'           , help='Do not submit batch job.')
(opt, args) = parser.parse_args()

nevents=opt.nevts

for en in [5,10,25,50,75,100,150,200,300]: #,500]:

    outDir='%s/version_%d/%s/e_%d'%(opt.out,opt.version,opt.gun,en)
    if opt.alpha>0 : outDir='%s_%3.3f'%(outDir,opt.alpha) 
    os.system('mkdir -p %s'%outDir)

    #wrapper
    scriptFile = open('%s/runJob.sh'%(outDir), 'w')
    scriptFile.write('#!/bin/bash\n')
    scriptFile.write('source %s/g4env.sh\n'%(os.getcwd()))
    scriptFile.write('cd %s\n'%(outDir))
    scriptFile.write('PFCalEE g4steer.mac %d\n'%opt.version)
    if len(opt.eos)>0:
        outTag='version_%d_e%d'%(opt.version,en)
        scriptFile.write('cmsStage PFcal.root %s/HGcal_%s.root'%(opt.eos,outTag))
    scriptFile.write('echo "All done"\n')
    scriptFile.close()

    #write geant 4 macro
    g4Macro = open('%s/g4steer.mac'%(outDir), 'w')
    g4Macro.write('/control/verbose 0\n')
    g4Macro.write('/control/saveHistory\n')
    g4Macro.write('/run/verbose 0\n')
    g4Macro.write('/event/verbose 0\n')
    g4Macro.write('/tracking/verbose 0\n')
    g4Macro.write('/generator/select particleGun\n')
    g4Macro.write('/gun/particle %s\n'%(opt.gun))    
    g4Macro.write('/gun/energy %f GeV\n'%(en))
    g4Macro.write('/gun/direction %f %f %f\n'%(math.cos(opt.alpha),math.sin(opt.alpha),0.))
    g4Macro.write('/run/beamOn %d\n'%(nevents))
    g4Macro.close()

    #submit
    os.system('chmod u+rwx %s/runJob.sh'%outDir)
    if opt.nosubmit : os.system('echo bsub -q %s %s/runJob.sh'%(opt.queue,outDir)) 
    else: os.system("bsub -q %s \'%s/runJob.sh\'"%(opt.queue,outDir))

