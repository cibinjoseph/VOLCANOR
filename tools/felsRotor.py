#!/usr/bin/python3

import numpy as np
import c81utils as c81
import parseResults as pr
import matplotlib.pyplot as plt

try:
    pr.ResultsDir = sys.argv[1]
except:
    pass

params = pr.getParams()
data = pr.getForceDist()

# Get variables to global scope
locals().update(params)
locals().update(data)

CLa_lin = 2.0*np.pi
alf0_deg = 0.0

alf0 = alf0_deg*np.pi/180.0

densityFELS = density
vTip = radius*Omega
dx = secArea/secChord
vinf = secSpan*Omega

alphalist = (180.0/np.pi)*(secCL/CLa_lin + alf0)
machlist = secVel/velSound
# print(alphalist)

fig, ax = plt.subplots(2)
ax[0].plot(secSpan/radius, alphalist, 'b*-', label='FELS')
ax[0].plot(secSpan/radius, secAlpha, 'r*-', label='Ind. vel.')
ax[0].legend()
ax[0].set_ylabel('Alpha (deg)')
ax[0].grid(True)

c81File = "NACA0012.C81"
with open(c81File, 'r') as fh:
    c81Airfoil = c81.load(fh)

CL_nonlin = []
for i, alpha in enumerate(alphalist):
    CL_nonlin.append(c81Airfoil.getCL(alpha, machlist[i]))

secLift = secCL*(0.5*densityFELS*secArea*vinf*vinf)
secLift_nonlin = CL_nonlin*(0.5*densityFELS*secArea*vinf*vinf)

ThrustFELS = nb*np.sum(secLift)
ThrustFELS_nonlin = nb*np.sum(secLift_nonlin)

CT = ThrustFELS / (densityFELS*np.pi*radius*radius*(vTip)**2.0)
CTFELS_nonlin = ThrustFELS_nonlin / (densityFELS*np.pi*radius*radius*(vTip)**2.0)

print('Min/Max alpha (deg) = ' + str(np.min(alphalist)) +' / ' + str(np.max(alphalist)))
print('Thrust = ' + str(ThrustFELS))
print('Thrust nonlin = ' + str(ThrustFELS_nonlin))

print('CT = ' + str(CTFELS))
print('CT nonlin = ' + str(CTFELS_nonlin))

# Write distribution to file
outMat = np.column_stack((secSpan/radius, alphalist, secLift/dx))
np.savetxt('loadVLM.dat', outMat, delimiter='  ')

# ax[1].plot(secSpan/radius, CL_nonlin, label='FELS')
# ax[1].plot(secSpan/radius, secLift/dx, label='FELS')
# ax[1].set_ylabel('Lift per unit span')

denom = density*(np.pi*radius*radius)*vTip*vTip
ax[1].plot(secSpan/radius, CL_nonlin, label='FELS')
ax[1].set_ylabel('CL (nonlinear)')
ax[1].grid(True)
ax[1].legend()
plt.xlabel('sec. span (r/R)')

plt.show()
