#!/usr/bin/python3
""" Extract pitch input and unsteady CT from case dir """

import numpy as np
import matplotlib.pyplot as plt
# import itertools
import sys

try:
    dirName = sys.argv[1] + '/'
except:
    dirName = './'

dt = 1.219509148257E-2
nt = 118
nb = 3
slopeCorrection = 5.73/(2*np.pi)

# Extract pitch input from log file
with open(dirName + 'volcanor.log') as fh:
    lines = fh.readlines()
pitch = []
count = -1
for line in lines:
    if ("pitch" in line):
        count += 1
        if count%nb == 0:
            pitch.append(float(line.split()[2]))

# Extract Ct data
ctmat = np.loadtxt(dirName + "Results/r01ForceNonDim.csv", skiprows=1)
ct = ctmat[:, 1]*slopeCorrection

t = np.arange(nt+1)*dt

fig, ax = plt.subplots(2, 1)

ax[0].plot(t, pitch)
ax[0].grid()

ax[1].plot(t, ct)
ax[1].grid()

fig.tight_layout()

# Write to tecplot file
print('TITLE = "Pitch ramp"')
print('VARIABLES = "t" "pitch" "CT"')
print('ZONE I=' + str(len(t)))
for ti, pi, cti in zip(t, pitch, ct):
    print(str(ti) + '\t' + str(pi) + '\t' + str(cti))


plt.show()
