#!/usr/bin/python

import visit
import argparse
from subprocess import call
from os import remove

parser = argparse.ArgumentParser(
        description=('Visualize plots using visit'),
        epilog='Author: Cibin Joseph')
parser.add_argument('-w', '--wake', help='Plot wake structure', action='store_true')
parser.add_argument('-t', '--tip', help='Plot wake tip', action='store_true')
parser.add_argument('-l', '--lift', help='Plot lift', action='store_true')
parser.add_argument('-d', '--drag', help='Plot drag', action='store_true')

args = parser.parse_args()

src_plot_dir = 'src_plot'

if args.wake == True:
    filename = 'plot_wake.py'
elif args.tip == True:
    filename = 'plot_tip.py'
elif args.lift == True:
    filename = 'plot_lift.py'
elif args.drag == True:
    filename = 'plot_drag.py'
else:
    print('Error: Wrong input arguments')
    raise ValueError

try:
    call(['visit','-s','{}/{}'.format(src_plot_dir,filename)])
    try:
        remove('visitlog.py')
    except OSError:
        pass
except KeyboardInterrupt:
    try:
        remove('visitlog.py')
    except OSError:
        pass
    print(' ...Program Exit')
