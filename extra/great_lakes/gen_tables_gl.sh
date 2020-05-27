#! /bin/bash
#
#setup_gl.sh
#Description:
#	This shell script does the following
#   - Load the desired modules

###############
## Constants ##
###############

module_warning=$'  (The "module" keyword is only available on the Great Lakes Computing Cluster\n  and thus will fail if called on any other computer.)';

###############
## Algorithm ##
###############

echo " "
echo "gen_tables_gl.sh"
echo " "
echo "Warning: This shell script is designed for the Great Lakes Unix System."
echo "Results not guaranteed for other operating systems."
echo " "

## Setup the Great Lakes Environment by Loading Modules
source load_modules_gl.sh
module list

echo "Modules successfully loaded."

## Run The MATLAB Script for ACC
cd matlab_files
matlab -batch "acc_generate_tables_gl"
