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
airfoilName = lines[0].strip()
for line in lines[1:]:
    cols = line.split()
    alpha.append(float(cols[0]))
    CL.append(float(cols[1]))

# Convert to C81 object
mach = [0.0, 0.5, 1.0]
CL = np.tile(CL, (3, 1)).T
c81airfoil = c81utils.C81(airfoilName, \
                         alpha, mach, CL, \
                         alpha, mach, CL, \
                         alpha, mach, CL)

# Extract filename without extension
filename = filename.split('.')[0]
with open(filename + '.C81', 'w') as fh:
    c81utils.dump(c81airfoil, fh)
