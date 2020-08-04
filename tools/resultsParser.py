""" Extract parameters from parameter file """
import numpy as np
from glob import glob

ResultsDir = 'Results/'

def getHeader(file):
    ''' Returns variables in header as list of strings '''
    with open(file, 'r') as fh:
        line = fh.readline()
    return line.split()

def getParams(file=ResultsDir + 'r01Params.dat'):
    """ Extract parameters from params file """
    with open(file, 'r') as fh:
        lines = fh.readlines()

    params = {}
    for line in lines:
        cols = line.split()
        params[cols[0]] = cols[1]
    return params

def getData(file):
    """ Parse File and extract header and data """
    mat = np.loadtxt(file, skiprows=1)
    data = {}
    for col, var in enumerate(getHeader(file)):
        data[var] = mat[:, col]
    return data

def getLatestFile(filename):
    """ Return filename of latest file with starting string """
    fileNumPrev = 0
    filenameLen = len(filename)
    for currentFile in glob(filename+'*'):
        fileNum = int(currentFile[filenameLen:filenameLen+5])
        if fileNum > fileNumPrev:
            latestFile = currentFile
            fileNumPrev = fileNum
    return latestFile

def getForceDist(file=None):
    """ Extract data from forcedist file """
    if file == None:
        # Parse ResultsDir and get latest ForceDist file
        file = getLatestFile(ResultsDir + 'r01b01ForceDist')
    return getData(file)


def getForceDim(file=ResultsDir + 'r01ForceDim.dat'):
    """ Extract data from ForceDim file """
    return getData(file)

def getForceNonDim(file=ResultsDir + 'r01ForceNonDim.dat'):
    """ Extract data from ForceDim file """
    return getData(file)
