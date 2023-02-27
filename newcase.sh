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

REPO_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

# Copy template files
echo "2. Copying template files"
cp -v ${REPO_DIR}/tools/template.case/*.nml $caseName/

# Create links
echo "3. Creating links to executables"
ln -v -s "${REPO_DIR}/bin" "$caseName/bin"
ln -v -s "${REPO_DIR}/tools" "$caseName/tools"
ln -v -s "${REPO_DIR}/tools/Makefile" "$caseName/Makefile"
