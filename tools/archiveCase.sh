#!/bin/bash
# Script to archive a current case

set -e

# zip -9 -r $1.zip $1
tar --use-compress-program="pigz -9 -k " $1.tar.gz $1
