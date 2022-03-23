import numpy as np
from scipy import integrate
import c81utils as c81
import parseResults as pr

params = pr.getParams()
data, forceDistFile = pr.getForceDist()

g = 9.81
m = 12.4
nb = 4.0
rho = 0.448
Inertia = 0.163
c = params['chord']

def getCL(alpha):
    # alpha in radians
    return 2*np.pi*alpha+0.141

def getCD(alpha):
    # alpha in radians
    return -1.02*np.cos(2.0*alpha-0.05)+1.028

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
dx = data['secArea']/data['secChord']
vFree = data['secSpan']*params['Omega']
alpha = data['secAlpha']
r = data['secSpan']

# Find CL CD distribution
CL = getCL(alpha*np.pi/180)
CD = getCD(alpha*np.pi/180)

# Integrate to obtain next omega and w
wNext = w + dt*getdw(CL, CD, r, w, Omega)
OmegaNext = Omega + dt*getdOmega(CL, CD, r, w, Omega)

# Write out next omega and w
print([-1.0*w, Omega, -1.0*wNext, OmegaNext])
dynDataOut = [-1.0*wNext, OmegaNext]
np.savetxt('dynamics.dat', dynDataOut, delimiter=' ')
