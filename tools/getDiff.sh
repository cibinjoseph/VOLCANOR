#! /bin/bash

# This script checks the difference between two files and
# reports if there is a mismatch. 
# This may be used during debugging to ensure changes do not cause
# differences in results.

# Check for input arguments
if [[ $# -eq 0 ]]; then
  testFile=wingNwake00020.plt
  inputFile=Results/wingNwake00020.plt
else
  testFile=$1
  inputFile=$2
fi

diff -q $testFile $inputFile && echo "NO CHANGES" || diff $testFile $inputFile
