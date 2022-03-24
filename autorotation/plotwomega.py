import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('r01bladedynamics.csv', delim_whitespace=True)
iterations = data['iter'].tolist()
w = data['w'].tolist()
omega = data['omega'].tolist()
revs = np.array(iterations)*5.0/360.0

plt.subplot(1, 2, 1)
plt.plot(revs, w, 'r-')
plt.xlabel('revs')
plt.ylabel('w (m/s)')
plt.grid()

plt.subplot(1, 2, 2)
plt.plot(revs, omega, 'b-')
plt.xlabel('revs')
plt.ylabel('omega (rad/s)')
plt.grid()

plt.show()
