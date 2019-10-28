#! /bin/bash

# This script copies CT from Results/r01forceHist*.txt to clipboard

# Check for input arguments
if [[ $# -eq 0 ]]; then
  lineFromLast=1
else
  lineFromLast=$1
fi

tail -n $lineFromLast Results/r01forceHist.txt | head -n 1 | awk '{print $2}' | xsel -ib
tail -n $lineFromLast Results/r01forceHist.txt | head -n 1
