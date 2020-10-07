import numpy as np
import matplotlib.pyplot as plt
import parseResults as pr

data = pr.getForceDist()
radius = 5.0
dx = data['secArea']/data['secChord']
density = 0.002377
forcedist = density*data['secLift']/dx
r = data['secSpan']/radius

plt.plot(r, forcedist, 'bo-')

outmat = np.stack((r, forcedist), axis=-1)
np.savetxt('forcesVOLCANOR.dat', outmat, delimiter='  ')

# plt.show()
