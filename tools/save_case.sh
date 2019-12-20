#!/bin/bash
# Code to save inputs and output from current case

set -e

if [ $# -eq 0 ]; then
  echo "Error: No arguments supplied"
  exit 1
fi

# Make folder provided as argument
recordDir=$1
mkdir -p $recordDir

# Copy input files
echo "Copying input files ..."
rsync -P *.in $recordDir

# Copy output files
echo "Copying output files ..."
rsync -rP Results $recordDir

# Get airfoil used
airfoilName="$(awk 'END {print $2}' rotor01.in)"
ls airfoils/$airfoilName > /dev/null 2>&1
if [ $? -eq 0 ]; then
  echo "Copying airfoil files ... $airfoilName"
  rsync -P airfoils/$airfoilName $recordDir
fi

# Get geometry used
geometryName="$(head -n3 rotor01.in | awk 'END {print $2}')"
ls geometry/$geometryName > /dev/null 2>&1
if [ $? -eq 0 ]; then
  echo "Copying geometry files ... $geometryName"
  rsync -P geometry/$geometryName $recordDir
fi

# Check for multiple rotor files
numOfrotorFiles="$(ls rotor??.in | wc -l)"
if [ $numOfrotorFiles -gt 1 ]; then
  echo "Warning: Multiple rotor files found. Copy others separately"
fi
