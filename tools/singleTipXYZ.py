#!/usr/bin/python3
""" Converts wakeTipXYZ plots to single blade's """

import sys
import os

try:
    filenames = sys.argv[1:]
except:
    raise ValueError('Wrong filename')

for filename in filenames:
    num =2+int(filename[-9:-4])
    os.system("head -" + str(num) + " " + filename + " | sponge " + filename)
