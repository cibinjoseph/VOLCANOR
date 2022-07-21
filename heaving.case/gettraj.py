import numpy as np
import matplotlib.pyplot as plt

ntInit = 300
nt = 720
dt = 2.77993E-3
f = 0.5

t = np.linspace(0, nt*dt, nt)

y = 1.0 + 0.5*np.sin(2*np.pi*f*t)
yp = 0.762*0.5*2*np.pi*f*np.cos(2*np.pi*f*t)

yInit = np.zeros((ntInit, 6))
yMat = np.zeros((nt, 6))
yMat[:, 2] = yp

outmat = np.vstack((yInit, yMat))
print(outmat.shape)
np.savetxt('trajectory01.in', outmat)

# plt.plot(t, y)
# plt.show()
