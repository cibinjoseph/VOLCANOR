#!/usr/bin/python3
""" Creates C81 file from alpha CL columns """
""" First line is taken as airfoil name """

import c81utils
import numpy as np
import matplotlib.pyplot as plt
import sys


# Read input file
filename = sys.argv[-1]
with open(filename, 'r') as fh:
    lines = fh.readlines()

# Extract airfoil polars
alpha = []
CL = []
CD = []
CM = []
airfoilName = lines[0].strip()

for line in lines[1:]:
    cols = line.split()
    lenCol = len(cols)
    if len(cols) > 1:
        alpha.append(float(cols[0]))
        CL.append(float(cols[1]))
        if len(cols) == 2:
            CD = CL
            CM = CL
        elif len(cols) == 3:
            CD.append(float(cols[2]))
            CM = CL
        elif len(cols) == 4:
            CD.append(float(cols[2]))
            CM.append(float(cols[3]))

# Convert to C81 object
mach = [0.0, 0.5, 1.0]
CL = np.tile(CL, (3, 1)).T
CD = np.tile(CD, (3, 1)).T
CM = np.tile(CM, (3, 1)).T
c81airfoil = c81utils.C81(airfoilName, \
                         alpha, mach, CL, \
                         alpha, mach, CD, \
                         alpha, mach, CM)

# Extract filename without extension
filename = filename.split('.')[0]
with open(filename + '.C81', 'w') as fh:
    c81utils.dump(c81airfoil, fh)

# Print status message
print('Airfoil name: ' + airfoilName)
print('C81 data written to file: ' + filename + '.C81')
