""" Extract parameters from parameter file """
import numpy as np
from glob import glob

def getParams(paramsFile='Results/r01Params.dat'):
    if paramsFile == None:
        paramsFile = 'Results/r01Params.dat'

    with open(paramsFile, 'r') as fh:
        lines = fh.readlines()

    params = {}
    for line in lines:
        cols = line.split()
        params[cols[0]] = cols[1]

    return paramsFile, params

def getForceDist(forceDistFile=None):
    if forceDistFile == None:
        # Parse Results/ and get latest ForceDist file
        fileNumPrev = 0
        for file in glob('Results/r01ForceDist*'):
            fileNum = int(file[20:25])
            if fileNum > fileNumPrev:
                forceDistFile = file
                fileNumPrev = fileNum

    with open(forceDistFile, 'r') as fh:
        lines = fh.readlines()

    secSpan = []
    secCL = []
    secCD = []
    secArea = []
    secChord = []
    for line in lines[2:]:
        cols = line.split()
        secSpan.append(float(cols[1]))
        secCL.append(float(cols[2]))
        secCD.append(float(cols[3]))
        secArea.append(float(cols[4]))
        secChord.append(float(cols[5]))

    secSpan = np.array(secSpan)
    secCL = np.array(secCL)
    secCD = np.array(secCD)
    secArea = np.array(secArea)
    secChord = np.array(secChord)

    return forceDistFile, secSpan, secCL, secCD, secArea, secChord
