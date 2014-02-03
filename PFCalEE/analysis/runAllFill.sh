#!/bin/sh

for v in 20;
do

    for s in 0 1 2 3 4 5 6;
    do
	mkdir -p PLOTS/version_$v/scenario_$s
	./bin/fillHistos 0 /afs/cern.ch/work/a/amagnan/public/HGCalEEDigi/version_$v/scenario_$s/ DigiPFcal.root

    done

done
