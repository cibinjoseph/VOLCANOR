""" Extract parameters from parameter file """
import numpy as np

def getParams(paramsFile):
    with open(paramsFile, 'r') as fh:
        lines = fh.readlines()

    params = {}
    for line in lines:
        cols = line.split()
        params[cols[0]] = cols[1]

    return params

def getForceDist(forceDistFile):
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

    return secSpan, secCL, secCD, secArea, secChord
