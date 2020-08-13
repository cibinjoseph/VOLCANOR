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

# Get parameters to obtain navg
propParams = pr.getParams()

# Find inst. Tz of each propeller from ForceDim file
propLift = [[],[],[],[]]
for i in range(4):
    propLift[i] = pr.getForceDim(file=pr.ResultsDir+ \
                               'r0'+str(i+1)+'ForceDim.dat')['Lz']

# Find inst. Lz of wing from ForceDim file
wingLift = pr.getForceDim(file=pr.ResultsDir+'r05ForceDim.dat')['Lz']

# Sum them to find net inst. L(x,y,z)
netLift = propLift[0] + propLift[1] + propLift[2] + propLift[3] + \
        wingLift

# Find average lift for last 5 revolutions

# Print out values
