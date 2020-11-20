#!/usr/bin/python3

import numpy as np
import c81utils as c81
import parseResults as pr
import matplotlib.pyplot as plt

params = pr.getParams()
data = pr.getForceDist()

CLa_lin = 2.0*np.pi
alf0_deg = 0.0

alf0 = alf0_deg*np.pi/180.0
omega = params['Omega']
rho = params['density']
Rad = params['radius']
nb = params['nb']
try:
    velSound = data['velSound']
except KeyError:
    velSound = 330.0

rhoMars = rho
vTip = Rad*omega

secSpan = data['secSpan']
secCL = data['secCL']
secArea = data['secArea']
secAlpha = data['secAlpha']
secVel = data['secVel']
dx = secArea/data['secChord']

vinf = secSpan*omega

alphalist = (180.0/np.pi)*(secCL/CLa_lin + alf0)
machlist = secVel/velSound
# print(alphalist)

fig, ax = plt.subplots(2)
ax[0].plot(secSpan/Rad, alphalist, 'b*-', label='FELS')
ax[0].plot(secSpan/Rad, secAlpha, 'r*-', label='Ind. vel.')
ax[0].legend()
ax[0].set_ylabel('Alpha (deg)')
ax[0].grid(True)

c81File = "NACA0012.C81"
with open(c81File, 'r') as fh:
    c81Airfoil = c81.load(fh)

CL_nonlin = []
for i, alpha in enumerate(alphalist):
    CL_nonlin.append(c81Airfoil.getCL(alpha, machlist[i]))

secLift = CL_nonlin*(0.5*rhoMars*secArea*vinf*vinf)
ThrustMars = nb*np.sum(secLift)
CTMars = ThrustMars / (rhoMars*np.pi*Rad*Rad*(vTip)**2.0)
print('Min/Max alpha (deg) = ' + str(np.min(alphalist)) +' / ' + str(np.max(alphalist)))
print('Thrust = ' + str(ThrustMars))
print('CT = ' + str(CTMars))

# Write distribution to file
outMat = np.column_stack((secSpan/Rad, alphalist, secLift/dx))
np.savetxt('loadVLM.dat', outMat, delimiter='  ')

ax[1].plot(secSpan/Rad, secLift/dx, label='FELS')
ax[1].set_ylabel('Lift per unit span')
ax[1].grid(True)
ax[1].legend()
plt.xlabel('sec. span (r/R)')

plt.show()
