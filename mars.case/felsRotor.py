#!/usr/bin/python3

import numpy as np
import c81utils as c81
import parseResults as pr
import matplotlib.pyplot as plt

params = pr.getParams()
data = pr.getForceDist()

CLa_lin = 2.0*np.pi
# alf0_deg = -6.480218
alf0_deg = 0.0

alf0 = alf0_deg*np.pi/180.0
omega = float(params['Omega'])
rho = float(params['density'])
Rad = float(params['radius'])
nb = float(params['nb'])

rhoMars = 0.022
vTip = Rad*omega

secSpan = np.array(data['secSpan'])
secCL = np.array(data['secCL'])
secArea = np.array(data['secArea'])
secAlpha = np.array(data['secAlpha'])

vinf = secSpan*omega

alphalist = (180.0/np.pi)*(secCL/CLa_lin + alf0)
# print(alphalist)

fig, ax = plt.subplots(2)
ax[0].plot(secSpan/Rad, alphalist, label='FELS')
ax[0].plot(secSpan/Rad, secAlpha, label='Ind. vel.')
ax[0].legend()
ax[0].set_ylabel('Alpha (deg)')
ax[0].grid(True)

c81File = "NACA5605_short.C81"
with open(c81File, 'r') as fh:
    c81Airfoil = c81.load(fh)

CL_nonlin = []
for alpha in alphalist:
    CL_nonlin.append(c81Airfoil.getCL(alpha, 0.5))

secLift = CL_nonlin*(0.5*rhoMars*secArea*vinf*vinf)
ThrustMars = nb*np.sum(secLift)
CTMars = ThrustMars / (rhoMars*np.pi*Rad*Rad*(vTip)**2.0)
print('Min/Max alpha (deg) = ' + str(np.min(alphalist)) +' / ' + str(np.max(alphalist)))
print('Thrust in Mars = ' + str(ThrustMars))
print('CT in Mars = ' + str(CTMars))

# Write distribution to file
outMat = np.column_stack((secSpan/Rad, alphalist, secLift))
np.savetxt('forceDist.dat', outMat, delimiter=',')

ax[1].plot(secSpan/Rad, secLift)
ax[1].set_ylabel('Lift')
ax[1].grid(True)
plt.xlabel('sec. span (r/R)')

plt.show()
