#!/bin/bash
# Sets up files and links for a new case

set -e

# Make a new case directory
if [ $# -eq 0 ]; then
  caseName="new.case"
elif [ $1 == '-h' ]; then
  echo 'Usage: newcase.sh [caseName]'
  exit
else
  caseName=$1
fi

echo "1. Creating directories"
mkdir -v $caseName
mkdir -v "$caseName/Results"
mkdir -v "$caseName/Restart"

# Copy template files
echo "2. Copying template files"
cp -v tools/template.case/* $caseName/

# Create links
echo "3. Creating links to executables"
ln -v -s "$PWD/bin" "$caseName/bin"
ln -v -s "$PWD/tools" "$caseName/tools"
ln -v -s "$PWD/tools/Makefile" "$caseName/Makefile"

