#!/usr/bin/python
# code to quickly plot wake, tip, lift and drag curves using VisIt

import visit
import argparse
import os
from subprocess import call
from time import sleep

parser = argparse.ArgumentParser(
        description=('Visualize plots using visit'),
        epilog='Author: Cibin Joseph')
parser.add_argument('-w', '--wake', help='Plot wake structure', action='store_true')
parser.add_argument('-t', '--tip', help='Plot wake tip', action='store_true')
parser.add_argument('-p', '--panel', help='Plot wing alone', action='store_true')
parser.add_argument('-g', '--gamma', help='Plot gamma sectional', action='store_true')
parser.add_argument('-l', '--lift', help='Plot lift', action='store_true')
parser.add_argument('-d', '--drag', help='Plot drag', action='store_true')

args = parser.parse_args()

src_plot_dir = 'src_plot'

if args.wake == True:
    filename = 'plot_wake.py'

elif args.tip == True:
    filename = 'plot_tip.py'

elif args.panel == True:
    filename = 'plot_panel.py'

elif args.gamma == True:
    filename = 'plot_gamma.py'

elif args.lift == True:
    filename = 'plot_lift.py'
    try:
        if os.path.exists('Results/lift.curve') == False:
            print('Waiting for file creation...') 
        while os.path.exists('Results/lift.curve') == False:
            sleep(1)
    except KeyboardInterrupt:
        print(' ...Program Exit')

elif args.drag == True:
    filename = 'plot_drag.py'
    try:
        if os.path.exists('Results/drag.curve') == False:
            print('Waiting for file creation...') 
        while os.path.exists('Results/drag.curve') == False:
            print('Waiting for file creation...') 
            sleep(1)
    except KeyboardInterrupt:
        print(' ...Program Exit')

else:
    print('Error: Wrong input arguments')
    raise ValueError


try:
    call(['visit','-s','{}/{}'.format(src_plot_dir,filename)])
    try:
        os.remove('visitlog.py')
    except OSError:
        pass
except KeyboardInterrupt:
    try:
        os.remove('visitlog.py')
    except OSError:
        pass
    print(' ...Program Exit')
