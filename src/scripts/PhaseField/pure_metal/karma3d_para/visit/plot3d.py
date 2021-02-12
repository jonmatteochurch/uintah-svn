#!/bin/sh
VISIT_PATH="/home/jonmatteochurch/Developer/visit/trunk/opt"
VISIT_VERS="3.1.2/linux-x86_64"
PATH="$PATH:$VISIT_PATH/bin"
LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$VISIT_PATH/$VISIT_VERS/lib"
"exec" "python2" "$0" "$@"

import sys
import os

print VISIT_PATH+os.sep+VISIT_VERS+"/lib/site-packages"
sys.path.append(VISIT_PATH+os.sep+VISIT_VERS+"/lib/site-packages")
sys.path.append(VISIT_PATH+os.sep+VISIT_VERS+"/lib/site-packages/visit")

from visit import *

assert len(sys.argv) > 1
db = sys.argv[1]
out = os.path.basename(db);
out = out.replace(".uda","")
if len(sys.argv) > 2:
	out = sys.argv[2]
out = out.replace(".","_")
if out[-1].isdigit():
	out = out + "_"

Launch()
OpenDatabase(db + "/index.xml", 0)
n=TimeSliderGetNStates()
SetTimeSliderState(n-1)

AddPlot("Pseudocolor", "u/0")
AddOperator("ThreeSlice")

AddPlot("Contour", "psi/0")
ContourAtts = ContourAttributes()
ContourAtts.colorType = ContourAtts.ColorBySingleColor  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
ContourAtts.legendFlag = 0
ContourAtts.singleColor = (153, 153, 153, 255)
ContourAtts.contourValue = (0)
ContourAtts.contourMethod = ContourAtts.Value  # Level, Value, Percent
SetPlotOptions(ContourAtts)

View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (2, 1, 1)
View3DAtts.focus = (75, 75, 75)
View3DAtts.viewUp = (0, 0, 1)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 150
View3DAtts.nearPlane = -300
View3DAtts.farPlane = 300
View3DAtts.imagePan = (0.05, 0)
View3DAtts.centerOfRotation = (75, 75, 75)
SetView3D(View3DAtts)

LightAtts = LightAttributes(0)
LightAtts.brightness = .6
SetLight(0, LightAtts)

LightAtts.type = LightAtts.Ambient
SetLight(1, LightAtts)

AnnotationAtts = GetAnnotationAttributes()
AnnotationAtts.databaseInfoExpansionMode = AnnotationAtts.SmartDirectory
AnnotationAtts.userInfoFlag = 0
AnnotationAtts.axes3D.visible = 0
SetAnnotationAttributes(AnnotationAtts)
DrawPlots()

r=GetRenderingAttributes()
r.specularFlag = 1
r.specularCoeff = 0.3
r.specularPower = 5
SetRenderingAttributes(r)

SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.format = SaveWindowAtts.BMP
SaveWindowAtts.fileName = out
#SaveWindowAtts.width, SaveWindowAtts.height = 3000,3000
SaveWindowAtts.screenCapture = 1
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()
