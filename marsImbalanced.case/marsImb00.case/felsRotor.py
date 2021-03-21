#!/usr/bin/python3

import numpy as np
import c81utils as c81
import parseResults as pr
import matplotlib.pyplot as plt
import tabulate as tb

params = pr.getParams()
data = pr.getForceDist()
locals().update(params)
locals().update(data)

CLa_lin = 2.0*np.pi
alf0_deg = -6.480218

alf0 = alf0_deg*np.pi/180.0
vTip = radius*Omega
dx = secArea/secChord
vInf = secSpan*Omega

alphalist = (180.0/np.pi)*(secCL/CLa_lin + alf0)
machlist = secVel/velSound
# DEBUG
alphalist = secAlpha

fig, ax = plt.subplots(2)
ax[0].plot(secSpan/radius, alphalist, 'b*-', label='FELS')
ax[0].plot(secSpan/radius, secAlpha, 'r*-', label='Ind. vel.')
ax[0].legend()
ax[0].set_ylabel('Alpha (deg)')
ax[0].grid(True)

c81File = "NACA5605XFOIL.C81"
with open(c81File, 'r') as fh:
    c81Airfoil = c81.load(fh)

CL_nonlin = []
CD_nonlin = []
for i, alpha in enumerate(alphalist):
    CL_nonlin.append(c81Airfoil.getCL(alpha, machlist[i]))
    CD_nonlin.append(c81Airfoil.getCD(alpha, machlist[i]))

secLift = CL_nonlin*(0.5*density*secArea*vInf*vInf)
bladeThrust = np.sum(secLift)
ThrustMars = nb*np.sum(secLift)
CTMars = ThrustMars / (density*np.pi*radius*radius*(vTip)**2.0)
print('Min/Max alpha (deg) = ' + str(np.min(alphalist)) +' / ' + str(np.max(alphalist)))
print('Thrust in Mars = ' + str(ThrustMars))
print('CT in Mars = ' + str(CTMars))

# Write distribution to file
# outMat = np.column_stack((secSpan/radius, alphalist, secLift/dx))
# np.savetxt('loadVLM.dat', outMat, delimiter='  ')
outDict = {'rbyR': secSpan/radius, 'secArea': secArea, 'secAlpha': secAlpha, \
           'alphaLookup': alphalist, 'CL_lin': secCL, \
           'CL_nonlin': CL_nonlin, 'CD_nonlin': CD_nonlin}
outTable = tb.tabulate(outDict, headers='keys', tablefmt='tsv', \
                       showindex=False)
with open('loadVLM.dat', 'w') as fh:
    fh.write(outTable)

ax[1].plot(secSpan/radius, secLift/dx, label='FELS')
ax[1].set_ylabel('Lift per unit span')
ax[1].grid(True)
ax[1].legend()
plt.xlabel('sec. span (r/R)')

plt.show()
