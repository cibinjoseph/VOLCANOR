import numpy as np
import c81utils as c81

with open('Results/r01ForceDist00720.dat', 'r') as fh:
    lines = fh.readlines()[2:]

omega = 314.16
CLa_lin = 4.584
# CLa_lin = 5.8
CLo_lin = 0.8326
rho = 1.2
rhoMars = 0.0205
Rad = 0.518
vTip = Rad*omega

secSpan = []
secCL = []
secArea = []
for line in lines[0:17]:
    cols = line.split()
    secSpan.append(float(cols[1]))
    secCL.append(float(cols[2]))
    secArea.append(float(cols[4]))

secSpan = np.array(secSpan)
secCL = np.array(secCL)
secArea = np.array(secArea)

vinf = secSpan*omega

alphalist = (180.0/np.pi)*(secCL - CLo_lin)/CLa_lin

c81File = "NACA5605.C81"
with open(c81File, 'r') as fh:
    c81Airfoil = c81.load(fh)

CL_nonlin = []
for alpha in alphalist:
    CL_nonlin.append(c81Airfoil.getCL(alpha, 0.5))

secLift = CL_nonlin*(0.5*rhoMars*secArea*vinf*vinf)
ThrustMars = 2.0*np.sum(secLift)
CTMars = ThrustMars / (rhoMars*np.pi*Rad*Rad*(vTip)**2.0)
print('Thrust in Mars = ' + str(ThrustMars))
print('CT in Mars = ' + str(CTMars))
