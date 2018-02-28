#!/bin/bash
# Code to save parameters from current case

casefile='casefile.txt'

# Parameters from library.f90
head -n10 library.f90 | tail -n3 > dummyfile1

# Parameters from input file and switches.f90
cat <(echo '==> library.f90 <==') dummyfile1 <(echo)\
  <(echo '==> inputfile <==') inputfile <(echo)\
  <(echo '==> switches.f90  <==') switches.f90 > $casefile

# Cleanup
rm dummyfile1

# Delete comments
sed -i '/^!/ d' $casefile
sed -i 's/integer, parameter :://g' $casefile
