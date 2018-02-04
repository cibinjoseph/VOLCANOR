#!/usr/bin/python

import visit
import argparse

parser = argparse.ArgumentParser(
        description=("Visualize plots using visit"),
        epilog="Author: Cibin Joseph")
parser.add_argument("-w","--wake",help="Plot wake structure",action="store_true")
parser.add_argument("-t","--tip" ,help="Plot wake tip",action="store_true")
parser.add_argument("-l","--lift",help="Plot lift",action="store_true")
parser.add_argument("-d","--drag",help="Plot drag",action="store_true")

args = parser.parse_args()

if args.wake == True:
    filename="~/WorkInProgress/VLM_Wing2.0/Results/wake*.tec database"
elif args.tip == True:
    filename="./Results/tip*.tec database"
elif args.lift == True:
    filename="./Results/lift.curve"
elif args.drag == True:
    filename="./Results/drag.curve"
else:
    print("Error: Wrong input arguments")
    raise ValueError
print(filename)

visit.AddArgument("-cli -gui -forceinteractivecli")
visit.Launch()
visit.OpenDatabase(filename)#"./Results/lift.curve")

visit.AddPlot("Curve", "# Lift", 1, 1)
visit.DrawPlots()
raw_input("Press any key to exit")
#ResetView()
## Begin spontaneous state
#View3DAtts = View3DAttributes()
#View3DAtts.viewNormal = (0.00409547, -0.932535, 0.361057)
#View3DAtts.focus = (0.0209413, 0.0427835, -0.909041)
#View3DAtts.viewUp = (0.015927, 0.361075, 0.932401)
#View3DAtts.viewAngle = 30
#View3DAtts.parallelScale = 8.21639
#View3DAtts.nearPlane = -16.4328
#View3DAtts.farPlane = 16.4328
#View3DAtts.imagePan = (0, 0)
#View3DAtts.imageZoom = 1
#View3DAtts.perspective = 1
#View3DAtts.eyeAngle = 2
#View3DAtts.centerOfRotationSet = 0
#View3DAtts.centerOfRotation = (0.0209413, 0.0427835, -0.909041)
#View3DAtts.axis3DScaleFlag = 0
#View3DAtts.axis3DScales = (1, 1, 1)
#View3DAtts.shear = (0, 0, 1)
#View3DAtts.windowValid = 1
#SetView3D(View3DAtts)
## End spontaneous state
#
