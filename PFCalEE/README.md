# PFCalEE

Geant4 simulation of a Si-base sampling calorimeter

Check https://twiki.cern.ch/twiki/bin/view/CMS/HGCalPerformanceStudiesWithG4

Geometry implementation is instantiated by detector versions in an enum - cf. src/DetectorConstruction.cc and src/SamplingSection.cc

A small ntuple is stored with the energy deposits - cf. src/EventAction.cc 

When changing the ttree content, adding homemade classes rather than
simple objects, the classes should be added to userlib. The dictionary
for root to understand the classes also need to be remade. Use "make
dictionary" before make, inside of userlib/. Follow instructions in
https://twiki.cern.ch/twiki/bin/view/Sandbox/AnnemarieMagnanSandbox
for what needs to be put in the homemade classes for root to
understand them (see also example in class userlib/include/HGCSSSimHit.hh).

## Setup the environment (SLC6)

source g4env.sh

## Compile
mkdir userlib/lib
mkdir userlib/obj
make


## Submit in parallel the runs submitRun.py. Some examples are:

Absorber width scan

for i in `seq 0 6`; do python submitRun.py -v ${i} -g mu-; done
for i in `seq 0 6`; do python submitRun.py -v ${i}; done

Si width scan Uniform 1xX0

for i in `seq 7 12`; do python submitRun.py -v ${i}; done

VJ version

python submitRun.py -v 13

Si width scan on Calice_Pb

for i in `seq 14 19`; do python submitRun.py -v ${i}; done

HGCALEE

for i in `seq 20 21`; do python submitRun.py -v ${i}; done

python submitRun.py -v 20 -g mu-
python submitRun.py -v 20 -g e+
python submitRun.py -v 20 -g pi-

for i in `seq 100 105`; do python submitRun.py -v ${i} -e /store/cmst3/group/hgcal/Geant4 -o ~/work/ntuples -g e-; done
