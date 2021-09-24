""" Extract parameters from parameter file """
import numpy as np
import pandas as pd
from glob import glob
import json

ResultsDir = 'Results/'
rotorNum = '01'
bladeNum = '01'
iterNum = ''

def _getHeader(file):
    ''' Returns variables in header as list of strings '''
    with open(file, 'r') as fh:
        line = fh.readline()
    return line.split()

def _getDataDict(file):
    """ Parse File and extract header and dataDict """
    dataDict = pd.read_csv(file, delim_whitespace=True).to_dict(orient='list')
    for col in dataDict.keys():
        dataDict[col] = np.array(dataDict[col])
    return dataDict

def _getLatestFile(filename):
    """ Return filename of latest file with starting string """
    fileNumPrev = 0
    filenameLen = len(filename)
    filelist = glob(filename + '*')
    if len(filelist) > 1:
        for currentFile in filelist:
            try:
                fileNum = int(currentFile[filenameLen:filenameLen+5])
            except ValueError:
                fileNum = fileNumPrev
            if fileNum > fileNumPrev:
                latestFile = currentFile
                fileNumPrev = fileNum
    else:
        latestFile = filelist[0]
    return latestFile

def getParamsDat(file=None):
    """ Extract parameters from params dat file """
    if file == None:
        file = ResultsDir + 'r' + rotorNum + 'Params.dat'

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
        file = ResultsDir + 'r' + rotorNum + 'Params.json'

    with open(file, 'r') as fh:
        params = json.load(fh)
    return params

def getForceDist(file=None):
    """ Extract dataDict from forcedist file """
    if file == None:
        file = _getLatestFile(ResultsDir + \
                              'r' + rotorNum + \
                              'b' + bladeNum + 'ForceDist' + \
                              iterNum)
    return _getDataDict(file), file

def getForceDim(file=None):
    """ Extract dataDict from ForceDim file """
    if file == None:
        file = _getLatestFile(ResultsDir + \
                              'r' + rotorNum + 'ForceDim')
    return _getDataDict(file)

def getForceNonDim(file=None):
    """ Extract dataDict from ForceDim file """
    if file == None:
        file = _getLatestFile(ResultsDir + \
                              'r' + rotorNum + 'ForceNonDim')
    return _getDataDict(file)
