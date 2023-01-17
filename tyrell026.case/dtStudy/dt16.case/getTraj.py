import numpy as np

# For the trajectory h = h0 + hs * sin(om*t)
# and unit airfoil chord
h0 = 0.5
hs = 0.1
k = 0.53
vinf = 10.0
nondimdt = 1.0/16.0
dt = nondimdt/vinf
ncyc = 5

om = 2.0*vinf*k
timePeriod = 2.0*np.pi/om
totTime = ncyc * timePeriod

nt = int(totTime/float(dt))
t = np.arange(1, nt+1, 1)*dt
w = om*hs*np.cos(om*t)

for i in range(nt):
    print("-10.0 0.0 " + str(w[i]) + " 0.0 0.0 0.0")

# print("Total time = ", end="")
# print(totTime)
print("dt, nt = ", end="")
print([dt, nt])
