#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

dataFELS = np.loadtxt('vlm_climb.txt', delimiter=',')
secSpan = np.array(dataFELS[:, 0])
alpha = np.array(dataFELS[:, 1])
secForce = np.array(dataFELS[:, 2])
plt.plot(secSpan, secForce, 'ro-', label='VLM')

dataBEMT = np.loadtxt('bemt_climb.txt')
secSpanBEMT = np.array(dataBEMT[:, 0])
secAlphaBEMT = np.array(dataBEMT[:, 1])
secForceBEMT = np.array(dataBEMT[:, 2])
plt.plot(secSpanBEMT, secForceBEMT, 'bo-', label='BEMT')


plt.legend(['VLM', 'BEMT'])
plt.ylabel('Load per unit span')
plt.xlabel('r/R')
plt.grid()
plt.show()
