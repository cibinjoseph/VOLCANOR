import numpy as np

arr = np.load('NACA0002CLCD.npy')
alphaList =arr[:, 0]

def getCLCD(alphaRad):
    CL = np.interp(alphaRad, alphaList, arr[:, 1])
    CD = np.interp(alphaRad, alphaList, arr[:, 2])
    return CL, CD
