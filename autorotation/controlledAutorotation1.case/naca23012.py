import numpy as np

def getCLCD(alphaList):
    # alphaList in radians
    alphaMax = 12.0*np.pi/180.0
    alphaC = -1.5*np.pi/180.0
    CLa = 5.73
    CLMax = 1.32
    CL0 = 0.15
    CD0 = 0.0079

    CL = []
    CD = []
    for alpha in alphaList:
        CD.append((1.03-((1.03-CD0)*np.cos(2*(alpha + alphaC)))))
        CL.append(CLa*np.sin(alpha)*np.cos(alpha)+CL0)

    CL = np.clip(CL, 2.0*CL0-CLMax, CLMax)

    return np.array(CL), np.array(CD)
