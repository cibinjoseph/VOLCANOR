#!/usr/bin/python3
# code to quickly plot wake, tip, lift and drag curves using VisIt

import argparse
import signal
import os
import sys
from subprocess import call
from time import sleep

scriptDir = os.path.dirname(os.path.realpath(__file__)) + '/'
currentDir = os.getcwd() + '/'

src_plotDir = scriptDir + 'src_plot/'
resultsDir = currentDir + 'Results/'

# Functions to be invoked for asynchronous keyboard signals ctrl+Z and ctrl+C


def ctrlZ_func(signum, frame):
    print('Reloading plots...')
    os.execl(sys.executable, 'python', __file__, *sys.argv[1:])  # Rerun code


def ctrlC_func(signum, frame):
    try:
        os.remove('visitlog.py')
    except OSError:
        pass
    print('Program exit...')  # Exit program
    sys.exit(0)


def wait4file(filename):
    if os.path.exists(filename) == False:
        print('Waiting for file creation...')
    while os.path.exists(filename) == False:
        sleep(1)


# Attach signal to the respective signal handlers (functions)
signal.signal(signal.SIGTSTP, ctrlZ_func)
signal.signal(signal.SIGINT, ctrlC_func)

# Define input arguments
parser = argparse.ArgumentParser(
    description=('Visualize plots using visit'),
    epilog='Author: Cibin Joseph')
parser.add_argument(
    '-w', '--wake', help='Plot wake structure', action='store_true')
parser.add_argument(
    '-f', '--force', help='Plot rotor force', action='store_true')
parser.add_argument('-s', '--span', help='Plot blade force',
                    action='store_true')
parser.add_argument(
    '-i', '--inflow', help='Plot blade inflow', action='store_true')
parser.add_argument('-t', '--tip', help='Plot wake tip', action='store_true')
parser.add_argument('-p', '--panel', help='Plot wing alone',
                    action='store_true')
parser.add_argument(
    '-g', '--gamma', help='Plot gamma sectional', action='store_true')
parser.add_argument(
    '-a', '--alpha', help='Plot alpha sectional', action='store_true')
parser.add_argument('-l', '--lift', help='Plot lift', action='store_true')
parser.add_argument('-d', '--drag', help='Plot drag', action='store_true')
parser.add_argument(
    '-c', '--custom', help='Use custom script', action='store_true')
parser.add_argument('-r', '--resultsDir', help='Specify Results directory')

# Parse args and obtain arguments
args = parser.parse_args()

# Convert to dict type for iterating through
argsDict = vars(args)

if args.resultsDir:
    resultsDir = currentDir + args.resultsDir
    argsDict['resultsDir'] = False

print('Reading results from: ' + resultsDir)

# Obtain filename for first argument that is True
pyFilename = 'plot_wake.py'  # default plot

# Check if custom script is used
if argsDict['custom'] == True:
    src_plotDir = currentDir
    argsDict.pop('custom')

for argName in argsDict:
    if argsDict[argName] == True:
        pyFilename = 'plot_' + argName + '.py'
        break

if pyFilename == 'plot_lift.py':
    wait4file(resultsDir + 'lift.curve')

elif pyFilename == 'plot_drag.py':
    wait4file(resultsDir + 'drag.curve')

if pyFilename == 'plot_force.py':
    call(['gnuplot', src_plotDir + 'plot_force.py'])

else:
    call(['visit', '-np', '4', '-s', '{}/{}'.format(src_plotDir, pyFilename)])

try:
    os.remove('visitlog.py')
except OSError:
    pass
