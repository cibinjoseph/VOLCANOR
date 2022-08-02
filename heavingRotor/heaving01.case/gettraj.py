import numpy as np
import matplotlib.pyplot as plt

ntInit = 300
nt = 720*2
dt = 2.77993E-3
f = 0.5

t = np.linspace(0, nt*dt, nt)

y = 1.0 + 0.5*np.sin(2*np.pi*f*t)
dy = 0.762*0.5*2*np.pi*f*np.cos(2*np.pi*f*t)

yMat = np.zeros((ntInit+nt, 6))
dyMat = np.zeros((ntInit+nt, 6))

yMat[ntInit:, 2] = y
dyMat[ntInit:, 2] = dy

np.savetxt('trajectory01.in', dyMat)
np.savetxt('h01.in', yMat)

# plt.plot(t, y)
# plt.show()
