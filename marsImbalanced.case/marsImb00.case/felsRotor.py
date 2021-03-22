#!/usr/bin/python3

import numpy as np
import c81utils as c81
import parseResults as pr
import matplotlib.pyplot as plt
import tabulate as tb


CLa_lin = 2.0*np.pi

def getNonlinear(alf0_deg, paramsFile=None, forceDistFile=None, c81File=None):
    params = pr.getParams(paramsFile)
    data = pr.getForceDist(forceDistFile)

    alf0 = alf0_deg*np.pi/180.0
    vTip = params['radius']*params['Omega']
    dx = data['secArea']/data['secChord']
    vInf = data['secSpan']*params['Omega']

    alphaLookup = (180.0/np.pi)*(data['secCL']/CLa_lin + alf0)
    machlist = data['secVel']/params['velSound']

    with open(c81File, 'r') as fh:
        c81Airfoil = c81.load(fh)

    CL_nonlin = []
    CD_nonlin = []
    for i, alpha in enumerate(alphaLookup):
        CL_nonlin.append(c81Airfoil.getCL(alpha, machlist[i]))
        CD_nonlin.append(c81Airfoil.getCD(alpha, machlist[i]))

    secLift = CL_nonlin*(0.5*params['density']*data['secArea']*vInf*vInf)
    secDrag = CD_nonlin*(0.5*params['density']*data['secArea']*vInf*vInf)
    Thrust = params['nb']*np.sum(data['secLift'])
    CT = Thrust / (params['density']*np.pi*params['radius']**2.0*(vTip)**2.0)


    sectDict = {'rbyR': data['secSpan']/params['radius'], \
                     'secArea': data['secArea'], \
                     'secAlpha': data['secAlpha'], \
                     'alphaLookup': alphaLookup, 'CL_lin': data['secCL'], \
                     'CL_nonlin': CL_nonlin, 'secLift': secLift, \
                     'CD_nonlin': CD_nonlin, 'secDrag': secDrag}
    return sectDict, Thrust, CT, dx


sectDict, Thrust, CT, dx = getNonlinear(c81File='NACA5605XFOIL.C81', \
                                         alf0_deg=-6.480218)
locals().update(sectDict)
print('Min/Max alpha (deg) = ' + \
      str(np.min(alphaLookup)) +' / ' + str(np.max(alphaLookup)))
print('Blade thrust = ' + str(np.sum(secLift)))
print('Thrust = ' + str(Thrust))
print('CT = ' + str(CT))

# Write distribution to file
outTable = tb.tabulate(sectDict, headers='keys', tablefmt='tsv', \
                       showindex=False)
with open('loadVLM.dat', 'w') as fh:
    fh.write(outTable)

# Plots
fig, ax = plt.subplots(2)
ax[0].plot(sectDict['rbyR'], alphaLookup, 'b*-', label='FELS')
ax[0].plot(sectDict['rbyR'], secAlpha, 'r*-', label='Ind. vel.')
ax[0].legend()
ax[0].set_ylabel('Alpha (deg)')
ax[0].grid(True)

ax[1].plot(sectDict['rbyR'], secLift/dx, label='FELS')
ax[1].set_ylabel('Lift per unit span')
ax[1].grid(True)
ax[1].legend()
plt.xlabel('sec. span (r/R)')

plt.show()
