#!/bin/bash
# Script to archive a current case

set -e

# Compress case for archiving
tar -I "pigz -9 -k " -cf $1.tar.gz $1

# Alternative compression using zip
# zip -9 -r $1.zip $1

# Uncompress archived case
# unpigz $1
