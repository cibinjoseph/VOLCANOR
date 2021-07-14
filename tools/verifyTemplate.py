#!/usr/bin/python3
""" Checks input formats against files in template directory """

import sys
import time
from pathlib import Path


def getLineStats(filename):
    """ Returns lineStats list with number of inputs in each line
    0 implies a comment
    """
    with open(filename) as fh:
        lines = fh.readlines()

    lineStats = []
    versionNum = None
    for line in lines:
        if line[0] == '#':
            lineStats.append(0)
        else:
            if not versionNum:
                # versionNum is NONE on the first encounter
                versionNum = line.strip()
            numInputs = len(line.split())
            lineStats.append(numInputs)

    return lineStats, versionNum

def isValid(filename, refInput):
    """ Checks if file is of the same template as refInput """
    returnVal = True

    refLineStats, refVersion = getLineStats(refInput)
    lineStats, versionNum = getLineStats(filename)

    # Check version numbers
    if versionNum != refVersion:
        print('Version mismath', end=': ')
        print('Expected ' + str(refVersion) + ', found ' + str(versionNum))

    # Check number of lines (neglecting comments)
    if (len(refLineStats)-refLineStats.count(0)) != \
    (len(lineStats)-lineStats.count(0)):
        print('No. of lines mismatch')
        returnVal = False

    # Input file stats without comments
    refLineStatsInpts = [x for x in refLineStats if x != 0]

    refIndex = -1
    for i, numInputs in enumerate(lineStats):
        if numInputs != 0:  # Non comment
            refIndex += 1
            if numInputs != refLineStatsInpts[refIndex]:
                print('Mismatch on line ' + str(i+1), end=': ')
                print('Expected ' + str(refLineStatsInpts[refIndex]) + \
                      ', found ' + str(numInputs))
                returnVal = False

    return returnVal


def main():
    scriptPath = str(Path(__file__).parent)
    args = sys.argv[1:]

    for arg in args:
        if arg[-3:] == '.in':
            if 'gridconfig' in arg:
                pass
            elif 'config.in' in arg:
                templateFile = scriptPath + '/template.case/config.in'
                print()
                print(arg)
                if isValid(arg, templateFile):
                    print('MATCH')
            elif 'geom' in arg:
                templateFile = scriptPath + '/template.case/geom01.in'
                print()
                print(arg)
                if isValid(arg, templateFile):
                    print('MATCH')

if __name__ == "__main__":
    main()
