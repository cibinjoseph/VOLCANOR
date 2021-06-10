#!/usr/bin/python3
""" Checks input formats against files in template directory """

import sys
import time

# Version 0.8
geomInput = [1, 4, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 2, \
             6, 4, 4, 2, 1, 2, 3, 3, 3, 3, 4, 1, 2, 3]

# Version 0.1
configInput = [1, 3, 2, 2, 1, 2, 4, 3, 2, 1, 1, 1]

def isValid(filename, refInput):
    returnVal = True
    with open(filename) as fh:
        lines = fh.readlines()

    l = []
    isComment = []
    for line in lines:
        l.append(len(line.split()))
        if line[0] == '#':
            isComment.append(True)
        else:
            isComment.append(False)

    if isComment.count(False) != len(refInput):
        print('No. of lines mismatch')
        returnVal = False

    refIndex = 0
    for i in range(len(l)):
        if isComment[i]:
            pass
        else:
            if l[i] != refInput[refIndex]:
                print('Mismatch on line ' + str(i+1), end='. ')
                print('Only ' + str(l[i]) + ' of ' + \
                      str(refInput[refIndex]) + ' inputs found')
                returnVal = False
            refIndex += 1

    return returnVal


args = sys.argv[1:]
for arg in args:
    if arg[-3:] == '.in':
        if 'gridconfig' in arg:
            pass
        elif 'config.in' in arg:
            print(arg)
            if isValid(arg, configInput):
                print('MATCH')
        elif 'geom' in arg:
            print(arg)
            if isValid(arg, geomInput):
                print('MATCH')

