#!/usr/bin/python3
""" Computes thrust averaged over last 1 rev """

import parseResults as pr
import numpy as np
import sys

try:
    pr.ResultsDir = sys.argv[1]
except:
    pass

params = pr.getParams()
data = pr.getForceNonDim()
locals().update(params)

ntPerRev = np.floor(2*np.pi/(Omega*dt))
CT = data['CL/CT']
CTAvg = np.zeros(CT.shape)

for i, _ in enumerate(CT):
    if i >= ntPerRev:
        CTAvg[i] = sum(CT[int(i-ntPerRev):i])/(ntPerRev+1)
    else:
        CTAvg[i] = CT[i]
    rev = i/ntPerRev
    print(str(rev) + '  ' + str(CTAvg[i]))

