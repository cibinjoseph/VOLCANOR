import numpy as np
from scipy import integrate
import c81utils as c81
import parseResults as pr
import naca23012

params = pr.getParams()
data, forceDistFile = pr.getForceDist()

g = 9.81
m = 12.4
nb = params['nb']
rho = params['density']
Inertia = 0.163
c = params['chord']

def getdw(CL, CD, r, w, Omega):
    integrand = nb*(CL*Omega*r+CD*w)*rho*0.5*np.sqrt((r*Omega)**2.0+w**2.0)*c
    integralterm = integrate.simps(integrand, r)
    dw = g-integralterm/m
    return dw

def getdOmega(CL, CD, r, w, Omega):
    integrand = nb*(CL*w-CD*Omega*r)*rho*0.5*np.sqrt((r*Omega)**2.0+w**2.0)*c*r
    integralterm = integrate.simps(integrand, r)
    dOmega = integralterm/Inertia
    return dOmega

# Read dynamics inputs
with open('dynamics.dat', 'r') as fh:
    line = fh.readline().split()
    dynDataIn = np.array(line, dtype='float64')
    w, Omega, dt = dynDataIn
    w = abs(w)

# Get alpha distribution
vFree = data['secSpan']*params['Omega']
alpha = data['secAlpha']
r = data['secSpan']

# Find CL CD distribution
CL, CD = naca23012.getCLCD(alpha*np.pi/180)

# Integrate to obtain next omega and w
wNext = w + dt*getdw(CL, CD, r, w, Omega)
OmegaNext = Omega + dt*getdOmega(CL, CD, r, w, Omega)

# Write out next omega and w
print([-1.0*w, Omega, -1.0*wNext, OmegaNext])
dynDataOut = [-1.0*wNext, OmegaNext]
np.savetxt('dynamics.dat', dynDataOut, delimiter=' ')
