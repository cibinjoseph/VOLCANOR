#!/usr/bin/python3
""" Computes FELS CL and CD from results file """

import numpy as np
import c81utils as c81
import parseResults as pr
import matplotlib.pyplot as plt
import tabulate as tb
import sys
import argparse


CLa_lin = 2.0*np.pi

parser = argparse.ArgumentParser(description='Force estimation from \
                                 linear secional lift')
parser.add_argument('-c', '--c81', action='store', help='C81 airfoil file')
parser.add_argument('-d', '--dir', default=pr.ResultsDir, metavar='Results/', \
                    action='store', help='Directory')
parser.add_argument('-r', '--rotor', default=pr.rotorNum, metavar='XX', \
                    action='store', help='Rotor num as string "XX"')
parser.add_argument('-b', '--blade', default=pr.bladeNum, metavar='XX', \
                    action='store', help='Blade num as string "XX"')
parser.add_argument('-i', '--iter', default=pr.iterNum, metavar='XXXXX', \
                    action='store', help='Iteration as string "XXXXX"')
parser.add_argument('-q', '--quiet', action='store_true', \
                    help='Suppress plots')

args = parser.parse_args()
print(args)

pr.ResultsDir = args.dir
pr.bladeNum = args.blade
pr.rotorNum = args.rotor
pr.iterNum = args.iter
c81File = args.c81

params = pr.getParams()
data = pr.getForceDist()

# Check if rotor or wing
isRotor = False
if abs(params['Omega']) > 0:
    isRotor = True

alf0 = params['alpha0']*np.pi/180.0
dx = data['secArea']/data['secChord']

if isRotor:
    vTip = params['radius']*params['Omega']
    vInf = data['secSpan']*params['Omega']
else:
    vFreestream = params['u']
    vInf = data['secVel']

alphaLookup = (180.0/np.pi)*(data['secCL']/CLa_lin + alf0)
machlist = data['secVel']/params['velSound']

if c81File == None:
    c81File = params['airfoilFile']

with open(c81File, 'r') as fh:
    c81Airfoil = c81.load(fh)

CL_nonlin = []
CD_nonlin = []
for i, alpha in enumerate(alphaLookup):
    CL_nonlin.append(c81Airfoil.getCL(alpha, machlist[i]))
    CD_nonlin.append(c81Airfoil.getCD(alpha, machlist[i]))

secLift_nonLin = CL_nonlin*(0.5*params['density']* \
                            data['secArea']*vInf*vInf)
secDrag_nonLin = CD_nonlin*(0.5*params['density']* \
                            data['secArea']*vInf*vInf)

if isRotor:
    secTorque = data['secSpan']*secDrag_nonLin
    Torque = np.sum(secTorque)
    denom = params['density']*np.pi*params['radius']**2.0*(vTip)**2.0
    CQ = Torque / (denom*params['radius'])
else:
    Drag = np.sum(secDrag_nonLin)
    denom = 0.5*params['density']*params['chord']*params['radius']*(vFreestream)**2.0
    CD = Drag / denom

Thrust = np.sum(secLift_nonLin)
CT = Thrust / denom

# Print inputs
print('C81 file = ' + c81File)
print()

# Print stats
print('Min/Max alpha (deg) = ' + \
      str(np.min(alphaLookup)) +' / ' + str(np.max(alphaLookup)))
print('Thrust = ' + str(Thrust))
print('CT = ' + str(CT))
if isRotor:
    print('Torque0 = ' + str(Torque))
    print('CQ0 = ' + str(CQ))
    print('Nb x Thrust = ' + str(params['nb']*Thrust))
    print('Nb x CT = ' + str(params['nb']*CT))
else:
    print('Drag0 = ' + str(Drag))
    print('CD0 = ' + str(CD))

# Write distribution to file
sectDict = {'rbyR': data['secSpan']/params['radius'], \
            'secArea': data['secArea'], \
            'secVel': vInf, 'dx': dx, \
            'secAlpha': data['secAlpha'], \
            'alphaLookup': alphaLookup, 'CL_lin': data['secCL'], \
            'CL_nonlin': CL_nonlin, 'secLift': secLift_nonLin, \
            'CD0_nonlin': CD_nonlin, 'secDrag0': secDrag_nonLin \
           }
outTable = tb.tabulate(sectDict, headers='keys', tablefmt='tsv', \
                       showindex=False)
with open('felsRotor.csv', 'w') as fh:
    fh.write(outTable)

# Plots
if args.quiet == False:
    fig, ax = plt.subplots(2)
    ax[0].plot(sectDict['rbyR'], alphaLookup, 'b*-', label='FELS')
    ax[0].plot(sectDict['rbyR'], data['secAlpha'], 'r*-', label='Ind. vel.')
    ax[0].legend()
    ax[0].set_ylabel('Alpha (deg)')
    ax[0].grid(True)

    ax[1].plot(sectDict['rbyR'], secLift_nonLin/dx, label='FELS')
    ax[1].set_ylabel('Lift per unit span')
    ax[1].grid(True)
    ax[1].legend()
    plt.xlabel('sec. span (r/R)')

    plt.show()
