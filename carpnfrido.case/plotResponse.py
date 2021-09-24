#!/usr/bin/python3
""" Extract pitch input and unsteady CT from case dir """

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import parseResults as pr


def normalized(l):
    return l/(l[-1])

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
dpitch = params['dpitch']
slopeCorrection = 5.73/(2*np.pi)

# Extract pitch input from dynamics file
pitchData = pr._getDataDict(pr.ResultsDir + 'r01bladedynamics.csv')
pitch = pitchData['pitch']
flap = pitchData['flap']

# Extract Ct data
ctdata = pr.getForceNonDim()
ct = ctdata['CL/CT']

t = np.arange(nt+1)*dt

fig, ax = plt.subplots(3, 1)
plt.title('Ramp = ' + str(dpitch) + ' deg/s')

ax[0].plot(t, pitch)
ax[0].grid()

ax[1].plot(t, normalized(flap))
ax[1].grid()
try:
    expt = pr._getDataDict(dirName + 'exptFlap.csv')
    ax[1].plot(expt['t'], normalized(expt['flap']), 'ro')
except:
    print('No exptFlap.csv data available')
    pass

ax[2].plot(t, normalized(ct))
ax[2].grid()
try:
    expt = pr._getDataDict(dirName + 'exptCT.csv')
    ax[2].plot(expt['t'], normalized(expt['CT']), 'ro')
except:
    print('No exptCT.csv data available')
    pass

fig.tight_layout()
plt.show()

# Write to file
print('t     pitch     flap     CT')
for ti, pi, fl, cti in zip(t, pitch, flap, ct):
    print(str(ti) + '\t' + str(pi) +'\t' + str(fl) + '\t' + str(cti))


