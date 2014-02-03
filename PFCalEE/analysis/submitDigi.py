#!/usr/bin/env python

import os,sys
import optparse
import commands
import math

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-q', '--queue'      ,    dest='queue'              , help='batch queue'                  , default='1nh')
parser.add_option('-g', '--granularity',    dest='granularity'        , help='granularities \"layer_i-layer_j:factor,layer:factor,...\"'             , default='0-29:4')
parser.add_option('-n', '--noise'      ,    dest='noise'              , help='noise (in Mips) \"layer_i-layer_j:factor,layer:factor,...\"'           , default='0-29:0.1')
parser.add_option('-t', '--thresh'     ,    dest='threshold'          , help='threshold (in ADC) \"layer_i-layer_j:factor,layer:factor,...\"'            , default='0-29:25')
parser.add_option('-N', '--Nevts'      ,    dest='Nevts'              , help='number of events to process' , default=0,    type=int)
parser.add_option('-m', '--miptoadc'   ,    dest='MipToADC'           , help='adc counts per MIP' , default=50,    type=int)
parser.add_option('-r', '--randomseed' ,    dest='seed'               , help='random seed' , default=0,    type=int)
parser.add_option('-d', '--debug'      ,    dest='debug'              , help='debug output' , default=0,    type=int)
parser.add_option('-s', '--scenario'   ,    dest='scenario'           , help='Integer describing scenario.' , default=0,    type=int)
parser.add_option('-o', '--out'        ,    dest='out'                , help='output directory'             , default=os.getcwd() )
parser.add_option('-S', '--no-submit'  ,    action="store_true",  dest='nosubmit'           , help='Do not submit batch job.')
(opt, args) = parser.parse_args()

nevents=opt.Nevts

for version in [20]:
    
    for en in [5,10,25,50,75,100,150,200,300,500]:
    #for en in [5,10,25,50,75,100]:
    #for en in [150,200,300,500]:
        
        inDir='%s/../version_%d/e-/e_%d/'%(os.getcwd(),version,en)
        outDir='%s/version_%d/scenario_%d/e-/e_%d'%(opt.out,version,opt.scenario,en)
        outlog='%s/digitizer.log'%(outDir)
        os.system('mkdir -p %s'%outDir)

        #wrapper
        scriptFile = open('%s/runJob.sh'%(outDir), 'w')
        scriptFile.write('#!/bin/bash\n')
        scriptFile.write('cd  %s/../\n'%(os.getcwd()))
        scriptFile.write('source g4env.sh\n')
        scriptFile.write('cd %s\n'%(outDir))
        scriptFile.write('%s/bin/digitizer %d %s %s %s %s %s %d %d %d | tee %s\n'%(os.getcwd(),opt.Nevts,inDir,outDir,opt.granularity,opt.noise,opt.threshold,opt.MipToADC,opt.seed,opt.debug,outlog))
        scriptFile.write('echo "All done"\n')
        scriptFile.close()
        
        #submit
        os.system('chmod u+rwx %s/runJob.sh'%outDir)
        #
        if opt.nosubmit : os.system('echo bsub -q %s %s/runJob.sh'%(opt.queue,outDir)) 
        else: os.system("bsub -q %s \'%s/runJob.sh\'"%(opt.queue,outDir))
