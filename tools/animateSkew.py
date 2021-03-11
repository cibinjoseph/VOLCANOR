#!/usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import os

nframes = len(os.listdir('Results'))
fps = 20
ntForOneRev = 24

plt.style.use('seaborn-pastel')
fig = plt.figure()
ax = plt.axes(xlim=(0, (nframes+5)/ntForOneRev), ylim=(0, 1.0))
line, = ax.plot([], [], 'r-', lw=3)

ratio = 0.3
xleft, xright = ax.get_xlim()
ybottom, ytop = ax.get_ylim()
ax.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)

def init():
    line.set_data([], [])
    return line,
def animate(i):
    print(str(i+1) + '/' + str(nframes))
    n = i + 1
    filename = 'Results/r01b01skew' + str(n).zfill(5) + '.dat'
    mat = np.loadtxt(filename)
    if mat.size != 2:
        y = mat[:, 1]
    else:
        y = [0.0]
    x = np.arange(len(y))/ntForOneRev
    line.set_data(x, y)
    return line,

anim = FuncAnimation(fig, animate, init_func=init,
                               frames=nframes, interval=1000/fps, blit=True)


anim.save('skew.gif', writer='imagemagick')
