import matplotlib.pyplot as plt
import numpy as np
from glob import glob

def getData(filename):
    """ Reads file and returns data """
    data = np.loadtxt(filename, skiprows=1)
    return data[:, 0], data[:, 1]

def getIterMax(filename):
    with open(filename) as fh:
        lines = fh.read().splitlines()
        last_line = lines[-1]
        iterMax = last_line.split()[-1]
        return iterMax

# Initial plot
it, CT = getData('Results/r01ForceNonDim.csv')
itMax = getIterMax('volcanor.log')

plt.ion()
plt.plot(it, CT)
plt.xlabel('iteration')
plt.ylabel('CT')

filelist = glob('Results/r??ForceNonDim.csv')

# Plot update
while True:
    for file in filelist:
        it, CT = getData(file)
        itLims = [it[0], it[-1]]
        CTmax = CT[-1]*1.05
        CTmin = CT[-1]*0.95
        plt.plot(it, CT, label=file[8:11])
        plt.legend()
        plt.title('Iteration: '+ str(int(it[-1])) + '/' + str(itMax))
        plt.plot(itLims, [CTmin, CTmin], 'r-', alpha = 0.3)
        plt.plot(itLims, [CTmax, CTmax], 'r-', alpha = 0.3)
    plt.pause(1.0)
    plt.clf()

plt.show(block=True)
