#!/usr/bin/python3

import numpy as np
import parseResults as pr
import matplotlib.pyplot as plt

# read alphadist
data = pr._getDataDict(pr.ResultsDir + 'r01b01alphaDist00300.curve')
secSpan = np.array(data['secSpan'])
secAlpha = np.array(data['alpha'])
plt.plot(secSpan, secAlpha, 'b-.', label='IndVel')

# read forcedist
data = pr.getForceDist()
secSpan = np.array(data['secSpan'])
secCL = np.array(data['secCL'])
CLa = 6.2832  # per rad
CL0 = 0.7106
secAlphaCalc = (secCL-CL0)/CLa*(180.0/np.pi)
plt.plot(secSpan, secAlphaCalc, 'r-.', label='FELS')
print('Approx. alpha = ' + str(secAlphaCalc[int(len(secAlphaCalc)/2)]))

# read vsp forcedist
vspMat = np.loadtxt('vsp.dat')
yavg = vspMat[:, 0]
CL = vspMat[:, 1]
secAlphaVSP = (CL-CL0)/CLa*(180.0/np.pi)
# plt.plot(yavg, secAlphaVSP, 'g.', label='VSP')

plt.ylabel('sectional alpha (deg)')
plt.grid(True)
plt.legend()
plt.show()
