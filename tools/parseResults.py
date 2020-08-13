""" Extract parameters from parameter file """
import numpy as np
from glob import glob

ResultsDir = 'Results/'

def _getHeader(file):
    ''' Returns variables in header as list of strings '''
    with open(file, 'r') as fh:
        line = fh.readline()
    return line.split()

def _getDataDict(file):
    """ Parse File and extract header and dataDict """
    mat = np.loadtxt(file, skiprows=1)
    dataDict = {}
    for col, var in enumerate(_getHeader(file)):
        dataDict[var] = mat[:, col]
    return dataDict

def _getLatestFile(filename):
    """ Return filename of latest file with starting string """
    fileNumPrev = 0
    filenameLen = len(filename)
    for currentFile in glob(filename+'*'):
        fileNum = int(currentFile[filenameLen:filenameLen+5])
        if fileNum > fileNumPrev:
            latestFile = currentFile
            fileNumPrev = fileNum
    return latestFile

def getParams(file=ResultsDir + 'r01Params.dat'):
    """ Extract parameters from params file """
    with open(file, 'r') as fh:
        lines = fh.readlines()

    params = {}
    for line in lines:
        cols = line.split()
        params[cols[0]] = cols[1]
    return params

def getForceDist(file=None):
    """ Extract dataDict from forcedist file """
    if file == None:
        # Parse ResultsDir and get latest ForceDist file
        file = _getLatestFile(ResultsDir + 'r01b01ForceDist')
    return _getDataDict(file)


def getForceDim(file=ResultsDir + 'r01ForceDim.dat'):
    """ Extract dataDict from ForceDim file """
    return _getDataDict(file)

def getForceNonDim(file=ResultsDir + 'r01ForceNonDim.dat'):
    """ Extract dataDict from ForceDim file """
    return _getDataDict(file)
