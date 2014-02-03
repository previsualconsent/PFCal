##################################
## Example analysis

#to compile
mkdir lib
mkdir bin
mkdir obj

make

#any .cc and .hh files will be compiled into obj files and into a library in lib/
#any .cpp file in test/ will be compiled as executable, in bin/

./bin/studyTransverseProperties -i ../ -v 20 -e 5
./bin/compareCaloStackPerformances -i ../ -v 20 -e 5

###################################
Info on existing Executables:
###################################

######################
## digitizer.cpp
# Submit script is submitDigi.py. Example:
##default baseline
python submitDigi.py -s 0 -o /afs/cern.ch/work/a/amagnan/public/HGCalEEDigi

#Note: the scenario number is decided when submitting the script, with the option -s, taking an integer:
#this is to avoid overwritting the output file ! A directory "scenario_i" is created.

# Example command lines for different scenarios in:
./runDigiForAll.sh


#######################
## fillHistos.cpp
# Run on G4 or on digitized files, to fill histograms, e.g. energy per layer, total energy, ...
# This is practical to then have just a plot macro to execute interactively in root to make nice plots.
#Example command line:
## on digi file
./bin/fillHistos 0 /afs/cern.ch/work/a/amagnan/public/HGCalEEDigi/version_$v/scenario_$s/ DigiPFcal.root
## on G4 file:
./bin/fillHistos 0 ../version_20/ PFcal.root

# Example submit script in:
./runAllFill.sh

# Example plotting macros in macros/plotE.C, plotXY.C, etc...
cd macros
root plotE.C++