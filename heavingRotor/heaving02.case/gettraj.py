import numpy as np
import matplotlib.pyplot as plt

ntInit = 300
nt = 720*2
dt = 2.77993E-3
R = 0.762

h0 = 1.0*R
h1 = 0.5*R
f = 1.0

t = np.linspace(0, nt*dt, nt)

y = h0 + h1*np.sin(2*np.pi*f*t)
dy = (2*np.pi*f)*h1*np.cos(2*np.pi*f*t)

yMat = np.zeros((ntInit+nt, 1))
dyMat = np.zeros((ntInit+nt, 6))

yMat[0:ntInit, 0] = h0/R
yMat[ntInit:, 0] = y/R
dyMat[ntInit:, 2] = dy

np.savetxt('trajectory01.in', dyMat)
np.savetxt('h01.in', yMat)

# plt.plot(t, y)
# plt.show()
