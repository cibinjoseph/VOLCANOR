#!/usr/bin/python3
import numpy as np
import parseResults as pr

# navg = 72
# omega = 400.3  # rad/s
# dt = 1.0E-3
# nt = round(2*np.pi/(omega*dt))
# Tavg = [0,0,0,0,0]
#
# for i in range(5):
#     mat = np.loadtxt('Results/r0' + str(i+1) + 'ForceDim.dat', skiprows=1)
#     T = mat[:, 1]
#     Tavg = np.sum(T[-navg:])/navg
#
#     print('Current iteration = ' + str(mat[-1,0]))
#     print('Average Thrust = ' + str(Tavg))

nrevs = 5

def averageOver(values, nSteps, nTimes):
    valuesFlipped = np.flip(values)
    avg = []
    if nTimes*nSteps+1 >= len(values):
        print('Warning: Cannot compute ' + str(nTimes) + \
              ' times. Switching to 1')
        nTimes = 1
    for i in range(nTimes):
        iStart = i*nSteps
        iEnd = iStart + nSteps
        avg.append(np.average(valuesFlipped[iStart:iEnd]))
    return avg


# Get parameters to obtain navg
propParams = pr.getParams()
dt = float(propParams['dt'])

# Find inst. Tz of each propeller from ForceDim file
propLift = []
nsteps = []
for i in range(4):
    propLift.append(pr.getForceDim(file=pr.ResultsDir+ \
                                   'r0'+str(i+1)+'ForceDim.dat')['Lz'])
    omega = float(pr.getParams(file=pr.ResultsDir + \
                         'r0'+str(i+1)+'Params.dat')['Omega'])
    nsteps.append(round(2*np.pi/(omega*dt)))

# Print Tz for props
print('Average Fz over last 1 rev')
print('Prop 1 : ' + str(averageOver(propLift[0], nsteps[0], 1)))
print('Prop 2 : ' + str(averageOver(propLift[1], nsteps[1], 1)))
print('Prop 3 : ' + str(averageOver(propLift[2], nsteps[2], 1)))
print('Prop 4 : ' + str(averageOver(propLift[3], nsteps[3], 1)))

# Find inst. Lz of wing from ForceDim file
wingLift = pr.getForceDim(file=pr.ResultsDir+'r05ForceDim.dat')['Lz']

# Sum them to find net inst. L(x,y,z)
netLift = propLift[0] + propLift[1] + propLift[2] + propLift[3] + \
        wingLift

# Find average lift for last 5 revolutions
print('Average Fz over last ' + str(nrevs) + ' revs')
print(averageOver(netLift, nsteps[0], nrevs))

