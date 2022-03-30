import numpy as np
from scipy import integrate
import c81utils as c81
import parseResults as pr
import naca0002
import os

params = pr.getParams()
data, forceDistFile = pr.getForceDist()
iternum = pr.getForceDim()['iter'][-1]

g = 8.87
m = 12.4
nb = params['nb']
Inertia = 0.163
c = params['chord']
dt = params['dt']
theta0 = params['theta0']
h0 = 80473  # Initial height (in m)

def getdw(CL, CD, r, w, Omega, vi, rho):
    integrand = nb*(CL*Omega*r+CD*(w-vi))*rho*0.5*np.sqrt((r*Omega)**2.0+(w-vi)**2.0)*c
    integralterm = integrate.simps(integrand, r)
    dw = g-integralterm/m
    return dw

def getdOmega(CL, CD, r, w, Omega, vi, rho):
    integrand = nb*(CL*(w-vi)-CD*Omega*r)*rho*0.5*np.sqrt((r*Omega)**2.0+(w-vi)**2.0)*c*r
    integralterm = integrate.simps(integrand, r)
    dOmega = integralterm/Inertia
    return dOmega

def getrho(h):
    # h in metres
    z = -abs(h)/1000-85
    a = [-1.7374e-03, -3.5660e-01, -1.0435e+01]
    polyz = (a[0]*(z**2)+a[1]*z+a[2])
    rho = 0.00390*np.exp(polyz)
    return rho

def getthetaRad(h):
    ztheta = np.load('ztheta.npy')
    # Piecewise interpolation using -h for increasing order
    thetaval = np.interp(-1*h, -1*ztheta[:, 0], ztheta[:, 1])
    return thetaval


def main():
# Read dynamics inputs
    with open('dynamics.dat', 'r') as fh:
        line = fh.readline().split()
        dynDataIn = np.array(line, dtype='float64')
        w, Omega = dynDataIn
        # Vertical axis convention is opposite
        w = -1.0*w

# Read current height from numpy array or initialize if iter = 0
    if iternum > 0:
        hHist = np.load('height.npy')
        h = hHist[-1]
    else:
        h = h0
        hHist = np.array([h0])
        try:
            os.remove('height.npy')
        except:
            pass

# Get alpha distribution
    vFree = data['secSpan']*params['Omega']
    alpha = data['secAlpha']
    r = data['secSpan']
# Vertical axis convention is opposite
    vi = -1.0*data['secViz']

# Find CL CD distribution
    CL, CD = naca0002.getCLCD(alpha*np.pi/180)

# Integrate to obtain next omega and w
    wNext = w + dt*getdw(CL, CD, r, w, Omega, vi, getrho(h))
    OmegaNext = Omega + dt*getdOmega(CL, CD, r, w, Omega, vi, getrho(h))
    hNext = h - dt*wNext
    thetaNext = getthetaRad(hNext)

# Bounds on variables
    # OmegaNext = min(OmegaNext, 200.0)

# Write out next omega and w
    print([-1.0*wNext, OmegaNext, hNext/1000, thetaNext*180.0/np.pi])
    dynDataOut = [-1.0*wNext, OmegaNext, thetaNext]
    np.savetxt('dynamics.dat', dynDataOut, delimiter=' ')

# Write hHist to numpy array
    hHist = np.append(hHist, hNext)
    np.save('height.npy', hHist)

if __name__ == "__main__":
    main()
