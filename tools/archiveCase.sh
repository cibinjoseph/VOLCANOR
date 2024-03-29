#!/bin/bash
# Script to archive a current case

set -e

# Compress case for archiving and display progress bar
# tar -I "pigz -9 -k " -cf $1.tar.gz $1 && pigz -tl $1.tar.gz
tar cf - $1 -P | pv -s $(du -sb $1 | awk '{print $1}') | pigz > $1.tar.gz && pigz -tl $1.tar.gz

# Alternative compression using zip
# zip -9 -r $1.zip $1

# Uncompress archived case
# unpigz $1
