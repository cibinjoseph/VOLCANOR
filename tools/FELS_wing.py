#!/usr/bin/python3
""" Computes FELS CL and CD from results file """
import sys
import c81utils
from resultsParser import getParams, getForceDist
import numpy as np

# Get filenames as arguments
args = sys.argv
if len(args) == 1:
    print('Usage: wingFELS.py a.C81 [CLline.dat] [r01ForceDistxx.dat] [r01Params.dat]')
    quit()
else:
    c81File = None
    paramsFile = None
    forceFile = None
    CLlineFile = None
    printSectional = False
    for arg in args[1:]:
        if '.C81' in arg:
            c81File = arg
        if 'Params' in arg:
            paramsFile = arg
        if 'ForceDist' in arg:
            forceFile = arg
        if 'line' in arg:
            CLlineFile = arg
        if '-s' in arg:
            printSectional = True

if c81File == None:
    print('Usage: wingFELS.py a.C81 [CLline.dat] [r01ForceDistxx.dat] [r01Params.dat]')
    quit()

paramsFile, params = getParams(paramsFile)
forceDistFile, secSpan, secCL, secCD, secArea, secChord = getForceDist(forceFile)

if (CLlineFile == None):
    CL0 = float(params['CL0'])
    CLa = float(params['CLa'])
else:
    with open(CLlineFile, 'r') as fh:
        CL0 = float(fh.readline())
        CLa = float(fh.readline())

refArea = float(params['radius'])*float(params['chord'])

with open(c81File, 'r') as fh:
    c81Airfoil = c81utils.load(fh)

secAlpha = []
secCL_nonLin = []
secCD_nonLin = []
for CL_Lin in secCL:
    alpha = (180.0/np.pi)*(CL_Lin - CL0)/CLa
    CL_nonLin = c81Airfoil.getCL(alpha, 0.1)
    CD_nonLin = c81Airfoil.getCD(alpha, 0.1)
    # Record for printing
    secAlpha.append(alpha)
    secCL_nonLin.append(CL_nonLin)
    secCD_nonLin.append(CD_nonLin)

# Add induced drag
secCD_nonLin += secCD

wingCL_Lin    = np.dot(np.array(secCL), secArea) / refArea
wingCD_Lin    = np.dot(np.array(secCD), secArea) / refArea

wingCL_nonLin = np.dot(np.array(secCL_nonLin), secArea) / refArea
wingCD_nonLin = np.dot(np.array(secCD_nonLin), secArea) / refArea

print('      c81File = ' + c81File)
print('   ParamsFile = ' + paramsFile)
print('ForceDistFile = ' + forceDistFile)
if (CLlineFile is not None):
    print('   CLlineFile = ' + CLlineFile)
print()
print('    theta0 = ' + str(params['theta0']))
print('    Lin CL = ' + str(wingCL_Lin))
print('    Lin CD = ' + str(wingCD_Lin))
print()
print('Non-lin CL = ' + str(wingCL_nonLin))
print('Non-lin CD = ' + str(wingCD_nonLin))
print()
if printSectional:
    for indx, alpha in enumerate(secAlpha):
        print((alpha, secCL_nonLin[indx], secCD_nonLin[indx]))
