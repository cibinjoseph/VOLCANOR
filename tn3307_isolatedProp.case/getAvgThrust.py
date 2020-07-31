#!/usr/bin/python3
import numpy as np

mat = np.loadtxt('Results/r01ForceDim.dat', skiprows=1)
navg = 72
T = mat[:, 1]
Tavg = np.sum(T[-navg:])/navg

print('Current iteration = ' + str(mat[-1,0]))
print('Average Thrust = ' + str(Tavg))
