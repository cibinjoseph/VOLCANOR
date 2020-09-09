#!/usr/bin/python3

import numpy as np
import parseResults as pr
import matplotlib.pyplot as plt

alpha0_deg = -6.480218

CLa = 2.0*np.pi  # per rad
alpha0 = alpha0_deg * np.pi/180.0

# read forcedist file
data = pr.getForceDist()
secSpan = np.array(data['secSpan'])
secCL = np.array(data['secCL'])
secAlpha = np.array(data['secAlpha'])

secAlphaEst = (secCL/CLa + alpha0)*(180.0/np.pi)

plt.plot(secSpan, secAlpha, 'b-.', label='IndVel')
plt.plot(secSpan, secAlphaEst, 'r-.', label='FELS')
print('Approx. alpha = ' + str(secAlphaEst[int(len(secAlphaEst)/2)]))

# read vsp forcedist
vspMat = np.loadtxt('vsp.dat')
yavg = vspMat[:, 0]
CL = vspMat[:, 1]
# secAlphaVSP = (CL-CL0)/CLa*(180.0/np.pi)
# plt.plot(yavg, secAlphaVSP, 'g.', label='VSP')

plt.ylabel('sectional alpha (deg)')
plt.grid(True)
plt.legend()
plt.show()
