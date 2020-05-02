#!/bin/bash
# Sets up files and links for a new case

set -e

# Make a new case directory
echo "1. Creating directories"
if [ $# -eq 0 ]; then
  caseName=caseNew
else
  caseName=$1
fi
mkdir -v $caseName
mkdir -v "$caseName/Results"

# Copy template files
echo "2. Copying template files"
cp -v tools/caseTemplate/* $caseName/

# Create links
echo "3. Creating links to executables"
if [ -f bin/volcanor ]; then
  ln -v -s "$PWD/bin/volcanor" "$caseName/volcanor"
else
  echo -e "\e[33mWarning\e[0m" ": volcanor executable not found in bin/
  \tBuild target 'run' in build/"
  ln -v -s "$PWD/bin/volcanor" "$caseName/volcanor"
fi

if [ -f bin/gridgen ]; then
  ln -v -s "$PWD/bin/gridgen" "$caseName/gridgen"
fi
ln -v -s "$PWD/tools/plotit.py" "$caseName/plotit.py"
