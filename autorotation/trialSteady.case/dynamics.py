import numpy as np
from scipy import integrate
import c81utils as c81
import parseResults as pr
import naca23012

params = pr.getParams()
data, forceDistFile = pr.getForceDist()

g = 8.87
m = 12.4
nb = params['nb']
rho = params['density']
Inertia = 0.163
c = params['chord']
dt = params['dt']
theta0 = params['theta0']

def getdw(CL, CD, r, w, Omega, vi):
    integrand = nb*(CL*Omega*r+CD*(w-vi))*rho*0.5*np.sqrt((r*Omega)**2.0+(w-vi)**2.0)*c
    integralterm = integrate.simps(integrand, r)
    dw = g-integralterm/m
    return dw

def getdOmega(CL, CD, r, w, Omega, vi):
    integrand = nb*(CL*(w-vi)-CD*Omega*r)*rho*0.5*np.sqrt((r*Omega)**2.0+(w-vi)**2.0)*c*r
    integralterm = integrate.simps(integrand, r)
    dOmega = integralterm/Inertia
    return dOmega

# Read dynamics inputs
with open('dynamics.dat', 'r') as fh:
    line = fh.readline().split()
    dynDataIn = np.array(line, dtype='float64')
    w, Omega = dynDataIn
    # Vertical axis convention is opposite
    w = -1.0*w

# Get alpha distribution
vFree = data['secSpan']*params['Omega']
alpha = data['secAlpha']
r = data['secSpan']
# Vertical axis convention is opposite
vi = -1.0*data['secViz']

# Find CL CD distribution
CL, CD = naca23012.getCLCD(alpha*np.pi/180)

# Integrate to obtain next omega and w
wNext = w + dt*getdw(CL, CD, r, w, Omega, vi)
OmegaNext = Omega + dt*getdOmega(CL, CD, r, w, Omega, vi)

# Write out next omega and w
print([-1.0*w, Omega, -1.0*wNext, OmegaNext])
dynDataOut = [-1.0*wNext, OmegaNext, theta0]
np.savetxt('dynamics.dat', dynDataOut, delimiter=' ')
