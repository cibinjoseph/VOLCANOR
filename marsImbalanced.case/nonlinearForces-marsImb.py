#!/usr/bin/python3
""" Computes CL and CD from results file using lookup tables """

import numpy as np
import c81utils as c81
import parseResults as pr
import matplotlib.pyplot as plt
import tabulate as tb
import sys
import argparse
from warnings import warn
from scipy import integrate

# Assumptions
# 1. No sideslip in secVel
# 2. No unsteady motion

CLa_lin = 2.0*np.pi

parser = argparse.ArgumentParser(description='Force estimation from \
                                 linear secional lift')
parser.add_argument('-a', '--all', action='store_true', help='On all results')
parser.add_argument('-b', '--blade', default=pr.bladeNum, metavar='XX', \
                    action='store', help='Blade num as string "XX"')
parser.add_argument('-c', '--c81', action='store', help='C81 airfoil file')
parser.add_argument('-d', '--dir', default=pr.ResultsDir, metavar='Results/', \
                    action='store', help='Directory')
parser.add_argument('-f', '--filealpha', action='store_true', \
                    help='Use alpha from input file')
parser.add_argument('-i', '--iter', default=pr.iterNum, metavar='XXXXX', \
                    action='store', help='Iteration as string "XXXXX"')
parser.add_argument('-o', '--out', action='store', default=None, \
                    help='Output nonlinear results to file')
parser.add_argument('-q', '--quiet', action='store_true', \
                    help='Suppress plots')
parser.add_argument('-r', '--rotor', default=pr.rotorNum, metavar='XX', \
                    action='store', help='Rotor num as string "XX"')
parser.add_argument('-s', '--steady', action='store_true', \
                    help='Rotor num as string "XX"')
parser.add_argument('-z', '--zero', action='store', \
                    help='Override zero lift angle in degs')

args = parser.parse_args()

pr.ResultsDir = args.dir
pr.bladeNum = args.blade
pr.rotorNum = args.rotor
pr.iterNum = args.iter
c81File = args.c81

params = pr.getParams()
data, forceDistFile = pr.getForceDist()

# Print input files
print('ForceDist file: ' + forceDistFile)

# Check if rotor or wing
isRotor = False
if abs(params['Omega']) > 0:
    isRotor = True

if args.zero is not None:
    alf0 = float(args.zero)*np.pi/180.0
    if alf0 > 0.0:
        warn("Changing sign of alpha0 to negative", stacklevel=2)
        alf0 = -1.0*alf0
else:
    alf0 = params['alpha0']*np.pi/180.0

dx = data['secArea']/data['secChord']

if isRotor:
    vRef = params['radius']*params['Omega']
    vRes = data['secVel']
    vFree = data['secSpan']*params['Omega']
else:
    vRef = np.linalg.norm([params['u'], params['v'], params['w']])
    vRes = abs(data['secVel'])
    vFree = abs(np.ones(vRes.shape)*vRef)

# for vRes1, vFree1 in zip(vRes, vFree):
#     if vRes1-vFree1 > 0:
#         phi.append(np.arctan2(np.sqrt(np.abs(vRes1**2-vFree1**2)), vFree1))
#     else:
#         warn('vRes < vFree: ' + str(vRes1) + ' < ' + str(vFree1), stacklevel=2)
#         phi.append(0)

# induced angle, phi
# Computing induced angle directly does not give correct results
# for rotary wings for some reason
# phi = data['secPhi']*np.pi/180.0
# phi here includes the climb/descent velocity besides the infow vel component
try:
    phi = (data['secTheta'] - data['secAlpha'])*np.pi/180.0
except:
    phi = data['secAlpha']*0.0

if args.filealpha:
    alphaLookup = data['secAlpha']
else:
    if args.steady:
        # Plus appears to be necessary at times
        CLvals = data['secCL']+data['secCLu']
    else:
        CLvals = data['secCL']
    alphaLookup = (180.0/np.pi)*(CLvals/CLa_lin + alf0)

# Override velSound with earth atmosphere since inside chamber
velSound = params['velSound']
machlist = vRes/velSound

if c81File == None:
    try:
        c81File = params['airfoilFile']
    except KeyError:
        pass

print('C81 file: ' + str(c81File))
print()

try:
    with open(c81File, 'r') as fh:
        c81Airfoil = c81.load(fh)

    # Lift and drag are w.r.t resultant velocity direction
    CL_nonlin = []
    CD0_nonlin = []
    for i, alpha in enumerate(alphaLookup):
        CL_nonlin.append(c81Airfoil.getCL(alpha, machlist[i]))
        CD0_nonlin.append(c81Airfoil.getCD(alpha, machlist[i]))
except (FileNotFoundError, TypeError):
    warn('c81File not found. Using linear curve.', stacklevel=2)
    CL_nonlin = data['secCL']
    CD0_nonlin = np.zeros(CL_nonlin.shape)

secLift_nonLin = CL_nonlin*(0.5*params['density']* \
                            data['secArea']*vRes*vRes)
secDrag0_nonLin = CD0_nonlin*(0.5*params['density']* \
                            data['secArea']*vRes*vRes)

# Assuming X along and Z perpendicular to freestream
secFz_nonLin = secLift_nonLin*np.cos(phi) - secDrag0_nonLin*np.sin(phi)
secFx_nonLin = secLift_nonLin*np.sin(phi) + secDrag0_nonLin*np.cos(phi)

if isRotor:
    secTorque = data['secSpan']*secFx_nonLin
    Torque = np.sum(secTorque)
    denom = params['density']*np.pi*params['radius']**2.0*(vRef)**2.0
    CQ = Torque / (denom*params['radius'])
    vi = np.tan(phi)*data['secSpan']*params['Omega']-params['w']
    viMean = 2*integrate.simps(vi*data['secSpan'], data['secSpan']) \
            /((1.0-(params['root_cut'])**2)*(params['radius'])**2.0)
else:
    Drag0 = np.sum(secDrag0_nonLin)
    vi = np.tan(phi)*vRef
    viMean = np.mean(vi)
    Fx = np.sum(secFx_nonLin)
    denom = 0.5*params['density']*params['chord']*params['radius']*(vRef)**2.0
    CFx = Fx / denom
    CD0 = Drag0 / denom

Fz = np.sum(secFz_nonLin)
CFz = Fz / denom

# Print stats
print('Min/Max alpha (deg) = ' + \
      str(np.min(alphaLookup)) +' / ' + str(np.max(alphaLookup)))
print('Fz = ' + str(Fz))
print('CFz = ' + str(CFz))
print('Mean downwash = ' + str(viMean))
if isRotor:
    print('Torque = ' + str(Torque))
    print('CQ = ' + str(CQ))
    print('Nb x Fz = ' + str(params['nb']*Fz))
    print('Nb x CFz = ' + str(params['nb']*CFz))
else:
    print('Fx = ' + str(Fx))
    print('CFx = ' + str(CFx))

# Write distribution to file
sectDict = {'rbyR': data['secSpan']/params['radius'], \
            'area': data['secArea'], \
            'velRes': vRes, 'dx': dx, 'velFree': vFree, \
            'alpha': data['secAlpha'], 'phi': phi, 'vi': vi, \
            'alphaLookup': alphaLookup, 'CL_lin': data['secCL'], \
            'CL_nonlin': CL_nonlin, 'Lift': secLift_nonLin, \
            'CD0_nonlin': CD0_nonlin, 'Drag0': secDrag0_nonLin, \
            'Fz_nonLin': secFz_nonLin, 'Fx_nonLin': secFx_nonLin
           }
outTable = tb.tabulate(sectDict, headers='keys', tablefmt='tsv', \
                       showindex=False)
if args.out:
    with open(args.out, 'w') as fh:
        fh.write(outTable)

# Plots
if args.quiet == False:
    fig, ax = plt.subplots(2)
    ax[0].plot(sectDict['rbyR'], alphaLookup, 'b*-', label='CL interpolated')
    ax[0].plot(sectDict['rbyR'], data['secAlpha'], 'r*-', label='Ind. vel.')
    ax[0].legend()
    ax[0].set_ylabel('Alpha (deg)')
    ax[0].grid(True)

    ax[1].plot(sectDict['rbyR'], secFz_nonLin/dx, label='CL interpolated')
    ax[1].set_ylabel('Lift per unit span')
    ax[1].grid(True)
    ax[1].legend()
    plt.xlabel('sec. span (r/R)')

    plt.show()
