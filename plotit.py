#!/usr/bin/python
# code to quickly plot wake, tip, lift and drag curves using VisIt

import visit
import argparse
import signal, os, sys
from subprocess import call
from time import sleep

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

# Attach signal to the respective signal handlers (functions)
signal.signal(signal.SIGTSTP, ctrlZ_func)
signal.signal(signal.SIGINT, ctrlC_func)

# Define input arguments
parser = argparse.ArgumentParser(
        description = ('Visualize plots using visit'), 
        epilog = 'Author: Cibin Joseph')
parser.add_argument('-w', '--wake', help='Plot wake structure', action = 'store_true')
parser.add_argument('-f', '--force', help='Plot rotor force', action = 'store_true')
parser.add_argument('-s', '--span', help='Plot blade force', action = 'store_true')
parser.add_argument('-i', '--inflow', help='Plot blade inflow', action = 'store_true')
parser.add_argument('-t', '--tip', help='Plot wake tip', action = 'store_true')
parser.add_argument('-p', '--panel', help='Plot wing alone', action = 'store_true')
parser.add_argument('-g', '--gamma', help='Plot gamma sectional', action = 'store_true')
parser.add_argument('-a', '--alpha', help='Plot alpha sectional', action = 'store_true')
parser.add_argument('-l', '--lift', help='Plot lift', action = 'store_true')
parser.add_argument('-d', '--drag', help='Plot drag', action = 'store_true')

# Parse args and obtain arguments
args = parser.parse_args()

# Convert to dict type for iterating through
argsDict = vars(args)

# Obtain filename for first argument that is True
filename = 'plot_wake.py'  # default plot
for argName in argsDict:
    if argsDict[argName] == True:
        filename = 'plot_'+argName+'.py'
        break

src_plot_dir = 'src_plot'

if filename == 'plot_lift.py':
    if os.path.exists('Results/lift.curve') == False:
        print('Waiting for file creation...') 
    while os.path.exists('Results/lift.curve') == False:
        sleep(1)

elif filename == 'plot_drag.py':
    if os.path.exists('Results/drag.curve') == False:
        print('Waiting for file creation...') 
    while os.path.exists('Results/drag.curve') == False:
        print('Waiting for file creation...') 
        sleep(1)

if filename == 'plot_force.py':
    sys.path.insert(0, 'src_plot')  # Append src_plot/ to search path
    import plot_force

else:
    call(['visit', '-np', '4', '-s', '{}/{}'.format(src_plot_dir, filename)])

try:
    os.remove('visitlog.py')
except OSError:
    pass
