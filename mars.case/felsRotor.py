#!/usr/bin/python3

import numpy as np
import c81utils as c81
import parseResults as pr
import matplotlib.pyplot as plt

params = pr.getParams()
data = pr.getForceDist()

CLa_lin = 6.283185
CLo_lin = 0.710635

omega = float(params['Omega'])
rho = float(params['density'])
Rad = float(params['radius'])
nb = float(params['nb'])

rhoMars = 0.022
vTip = Rad*omega

secSpan = np.array(data['secSpan'])
secCL = np.array(data['secCL'])
secArea = np.array(data['secArea'])

vinf = secSpan*omega

alphalist = (180.0/np.pi)*(secCL - CLo_lin)/CLa_lin
print(alphalist)
quit()

c81File = "NACA5605.C81"
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

# plt.plot(secSpan/Rad, secLift)
# plt.show()
