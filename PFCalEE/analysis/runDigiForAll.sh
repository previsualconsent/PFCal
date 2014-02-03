#!/bin/sh

##default baseline
python submitDigi.py -s 0 -o /afs/cern.ch/work/a/amagnan/public/HGCalEEDigi
##scan adc range
python submitDigi.py -s 1 -m 1 -t 0-29:1 -o /afs/cern.ch/work/a/amagnan/public/HGCalEEDigi
python submitDigi.py -s 2 -m 10 -t 0-29:5 -o /afs/cern.ch/work/a/amagnan/public/HGCalEEDigi
python submitDigi.py -s 3 -m 25 -t 0-29:13 -o /afs/cern.ch/work/a/amagnan/public/HGCalEEDigi

##scan granularity and noise
python submitDigi.py -s 4 -g 0-29:4 -n 0-29:0.1 -t 0-29:25 -o /afs/cern.ch/work/a/amagnan/public/HGCalEEDigi
python submitDigi.py -s 5 -g 0-29:2 -n 0-29:0.02 -t 0-29:6 -o /afs/cern.ch/work/a/amagnan/public/HGCalEEDigi
python submitDigi.py -s 6 -g 0-10:4,11-20:2,21-29:4 -n 0-10:0.1,11-20:0.02,21-29:0.1 -t 0-10:25,11-20:6,21-29:25 -o /afs/cern.ch/work/a/amagnan/public/HGCalEEDigi

