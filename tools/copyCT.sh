#! /bin/bash

# This script copies CT from Results/r01forceHist*.txt to clipboard

forceFile="r01ForceNonDim.dat"

# Check for input arguments
if [[ $# -eq 0 ]]; then
  resultsDir="Results"
else
  resultsDir=$1
fi

tail -n 1 $resultsDir/$forceFile | head -n 1 | awk '{print $2}' | xsel -ib
tail -n 1 $resultsDir/$forceFile | head -n 1
