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

DrawPlots()
SetActivePlots(0)
Query("SpatialExtents")
Extents = GetQueryOutputValue()
Width = (Extents[1]-Extents[0], Extents[3]-Extents[2], Extents[5]-Extents[4])
Center = ((Extents[1]+Extents[0])/2, (Extents[3]+Extents[2])/2, (Extents[5]+Extents[4])/2)

View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (2, 1, 1)
View3DAtts.focus = Center
View3DAtts.viewUp = (0, 0, 1)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = Width[0]
View3DAtts.nearPlane = -2*Width[0]
View3DAtts.farPlane = 2*Width[0]
View3DAtts.imagePan = (0.05, 0)
View3DAtts.centerOfRotation = Center
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
SaveWindowAtts.format = SaveWindowAtts.VTK
SaveWindowAtts.fileName = out
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()

SaveWindowAtts.format = SaveWindowAtts.BMP
SaveWindowAtts.fileName = out
SaveWindowAtts.screenCapture = 1
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()

SetActivePlots(1)
HideActivePlots()

AnnotationAtts.legendInfoFlag=0
AnnotationAtts.timeInfoFlag=0
AnnotationAtts.databaseInfoFlag=0
AnnotationAtts.axes3D.bboxFlag=0
AnnotationAtts.axes3D.triadFlag=0
SetAnnotationAttributes(AnnotationAtts)

View3DAtts.viewNormal = (0, 0, 1)
View3DAtts.viewUp = (0, 1, 0)
View3DAtts.parallelScale = Center[0]
View3DAtts.imagePan = (0, 0)
View3DAtts.perspective = 0
SetView3D(View3DAtts)

SaveWindowAtts.fileName = out + "_xy_"
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()

View3DAtts.viewNormal = (0, -1, 0)
View3DAtts.viewUp = (0, 0, 1)
SetView3D(View3DAtts)

SaveWindowAtts.fileName = out + "_xz_"
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()

View3DAtts.viewNormal = (1, 0, 0)
View3DAtts.viewUp = (0, 0, 1)
SetView3D(View3DAtts)

SaveWindowAtts.fileName = out + "_xz_"
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()
