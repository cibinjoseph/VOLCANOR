#!/usr/bin/python3

import glob, os, sys

dirList = sys.argv[1:]

cwd = os.getcwd()

for dirName in dirList:
    print(dirName + ' :')
    os.chdir(dirName)
    for fileName in glob.glob('rotor??.in'):
        # print(fileName)
        rotorNum = fileName[5:7]
        os.rename(fileName, 'geom'+rotorNum+'.in')
    os.chdir(cwd)

