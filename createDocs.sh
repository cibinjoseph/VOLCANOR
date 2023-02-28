#!/bin/bash

##############################################
# Run this in the topmost folder of the repo #
#       to generate the documentation        #
##############################################

rm -rf docs/src
if [[ $* == '--help' ]] || [[ $* == '-h' ]]; then
  echo 'Usage createDocs.sh [OPTION]'
  echo 'Generate documentation using FORD.'
  echo
  echo '  -h, --help     Display help' 
  echo '  -d, --debug    Generate docs without search' 
elif [[ $* == '--debug' ]] || [[ $* == '-d' ]]; then
  ford -p ford_input --no-search ford_input.md
else
  ford -p ford_input ford_input.md
fi
