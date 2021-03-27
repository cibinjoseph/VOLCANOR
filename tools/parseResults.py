""" Extract parameters from parameter file """
import numpy as np
from glob import glob
import json

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
        try:
            dataDict[var] = np.array(mat[:, col], dtype='float64')
        except:
            dataDict[var] = mat[:, col]
    return dataDict

def _getLatestFile(filename):
    """ Return filename of latest file with starting string """
    fileNumPrev = 0
    filenameLen = len(filename)
    filelist = glob(filename + '*')
    if len(filelist) > 1:
        for currentFile in filelist:
            fileNum = int(currentFile[filenameLen:filenameLen+5])
            if fileNum > fileNumPrev:
                latestFile = currentFile
                fileNumPrev = fileNum
    else:
        latestFile = filelist[0]
    return latestFile

def getParamsDat(file=None):
    """ Extract parameters from params dat file """
    if file == None:
        file = ResultsDir + 'r01Params.dat'

    with open(file, 'r') as fh:
        lines = fh.readlines()

    params = {}
    for line in lines:
        cols = line.split()
        try:
            params[cols[0]] = float(cols[1])
        except ValueError:
            params[cols[0]] = cols[1]
    return params

def getParams(file=None):
    """ Extract parameters from params json file """
    if file == None:
        file = ResultsDir + 'r01Params.json'

    with open(file, 'r') as fh:
        params = json.load(fh)
    return params

def getForceDist(file=None):
    """ Extract dataDict from forcedist file """
    if file == None:
        file = _getLatestFile(ResultsDir + 'r01b01ForceDist')
    return _getDataDict(file)

def getForceDim(file=None):
    """ Extract dataDict from ForceDim file """
    if file == None:
        file = _getLatestFile(ResultsDir + 'r01ForceDim')
    return _getDataDict(file)

def getForceNonDim(file=None):
    """ Extract dataDict from ForceDim file """
    if file == None:
        file = _getLatestFile(ResultsDir + 'r01ForceNonDim')
    return _getDataDict(file)
