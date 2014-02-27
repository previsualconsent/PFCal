#!/bin/bash
source /afs/cern.ch/user/p/phansen/public/PFCal/PFCalEE/g4env.sh
cd /afs/cern.ch/user/p/phansen/public/PFCal/PFCalEE/
analysis/bin/vectorfit 0 15 version_20/e-/ 0.250
#                  <nevents> <radius> <path to data> <angle in rads>
echo "All done"
