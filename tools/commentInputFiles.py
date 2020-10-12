#!/usr/bin/python3
""" Changes description lines to comments """
""" in all input files """

import subprocess
import os
import glob

def fileLength(filename):
    """ Returns number of lines in file """
    with open(filename) as f:
        for i, l in enumerate(f, 1):
            pass
    return i

def commentFile(filename, excludeLines):
    """ Comments file excluding excludeLines """
    for lineNum in range(1, fileLength(filename)+1):
        if lineNum not in excludeLines:
            # Comment non-empty line
            command = "sed -i '" + str(lineNum) + "s/./# &/' " + filename
            subprocess.call([command], shell=True)
            # Comment empty line
            command = "sed -i '" + str(lineNum) + "s/^$/# &/' " + filename
            subprocess.call([command], shell=True)

# config.in
excludeLines = [3, 7, 11, 16, 21, 27, 33, 39, 44, 49, 54]
commentFile('config.in', excludeLines)

# rotor.in
rotorFiles = glob.glob('rotor??.in')
excludeLines = [5, 9, 14, 18, 22, 26, 31, 36, \
                40, 45, 51, 57, 63, 67, 72, 76,
                82, 85, 88, 97, 104, 110, 115]
for filename in rotorFiles:
    commentFile(filename, excludeLines)

# gridconfig.in
if os.path.exists('./gridconfig.in'):
    excludeLines = [4, 9, 13, 18, 24]
    commentFile('gridconfig.in', excludeLines)
