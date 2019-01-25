#!/usr/bin/python

# time = timestep rpm ndeg [nrev]
# Calculates time(in s) taken for rotating 'ndeg' degrees, given rpm
import sys
import math

if len(sys.argv) != 3 and len(sys.argv) != 4:
    print('ERROR: Wrong number of arguments')
    print('USAGE:  timestep rpm ndeg [nrev]\n')
    raise RuntimeError

args=sys.argv
rpm =float(eval(args[1],{},{}))     # eval for evaluating minor expressions
ndeg=float(eval(args[2],{},{}))

omega=2.*math.pi*rpm/60.  # in rad/s
degs_time=(2.*ndeg*math.pi)/(omega*180.)

print('RPM              = {}'.format(rpm))
print('Omega (rad/s)    = {}'.format(omega))
print('Degrees          = {}'.format(ndeg))
print('Timestep (s)     = {}'.format(degs_time))

if (len(sys.argv) == 4):
    nrev=float(eval(args[3],{},{}))
    nt=nrev*360./ndeg
    print('Revs             = {}'.format(nrev))
    print('no. of timesteps = {}'.format(nt))
