#!/usr/bin/python3
""" Extract pitch input and unsteady CT from case dir """

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import parseResults as pr
from scipy import integrate


averageOver = 24*4

def normalize(l, npoints=averageOver):
    return np.array(l)*npoints/np.sum(l[-npoints:])

try:
    dirName = sys.argv[1] + '/'
except:
    dirName = './'

# Extract params
pr.ResultsDir = dirName + pr.ResultsDir
params = pr.getParams()

dt = params['dt']
nt = params['nt']
nb = params['nb']
omega = params['Omega']
dpitch = params['dpitch']
slopeCorrection = 5.73/(2*np.pi)

# Set whether flap is present
flapPresent = True
if params['bladeDynamicsSwitch'] == 0:
    flapPresent = False

# Extract pitch input from dynamics file
pitchData = pr._getDataDict(pr.ResultsDir + 'r01bladedynamics.csv')
if flapPresent:
    pitch = pitchData['pitch']
    flap = pitchData['flap']

# Extract Ct data
ctdata = pr.getForceNonDim()
ct = ctdata['CL/CT']
# print('CT_max/CT_ss = ', end='')
# print(max(ct)/ct[-1])

# Mean inlow 
viMean = []
for it in range(nt+1):
    pr.iterNum = f'{it:05d}'
    data, _ = pr.getForceDist()
    phi = data['secPhi']*np.pi/180.0
    vi = data['secSpan']*params['Omega']*np.tan(phi)
    viMeanVal = 2*integrate.simps(vi*data['secSpan'], data['secSpan']) \
            /((1.0-(params['root_cut'])**2)*(params['radius'])**2.0)
    viMean.append(viMeanVal)

pr.iterNum = ''


t = np.arange(nt+1)*dt

# Plots
generatePlot = False
if generatePlot == True:
    fig, ax = plt.subplots(3, 1)
    plt.rcParams['axes.grid'] = True

    ax[0].plot(t, normalize(viMean))
    ax[0].plot(t, np.ones(t.shape))
    ax[0].set_ylabel('inflow')
    ax[0].set_title('Ramp = ' + str(dpitch) + ' deg/s')
    try:
        expt = pr._getDataDict(dirName + 'exptVi.csv')
        ax[0].plot(expt['t'], normalize(expt['vi']), 'ro')
    except:
        print('No exptVi.csv data available', file=sys.stderr)
        pass

    ax[1].plot(t, normalize(flap))
    ax[1].plot(t, np.ones(t.shape))
    ax[1].set_ylabel('flap')
    try:
        expt = pr._getDataDict(dirName + 'exptFlap.csv')
        ax[1].plot(expt['t'], normalize(expt['flap']), 'ro')
    except:
        print('No exptFlap.csv data available', file=sys.stderr)
        pass

    ax[2].plot(t, normalize(ct))
    ax[2].plot(t, np.ones(t.shape))
    ax[2].set_ylabel('CT')
    try:
        expt = pr._getDataDict(dirName + 'exptCT.csv')
        ax[2].plot(expt['t'], normalize(expt['CT']), 'ro')
    except:
        print('No exptCT.csv data available', file=sys.stderr)
        pass

    fig.tight_layout()
    plt.show()

# Write to file
writeToFile = True
dataOut = {}
dataOut['t'] = t
dataOut['viMean'] = normalize(viMean)
dataOut['CT'] = normalize(ct)
if flapPresent:
    dataOut['pitch'] = normalize(pitch)
    dataOut['flap'] = normalize(flap)

if writeToFile:
    dataPd = pd.DataFrame.from_dict(dataOut)
    dataPd.to_csv(sys.stdout, sep='\t', index=False)
