# Cut the groove for Henry's HSC wordworking project.

# The whole cut is a rounded rectangle, but the length of it exceeds the CNC
# machine size, so we'll have to do it in two passes & carefully line them up.
# But obviously both halves of the cut look the same, so we'll just program one
# carve and use it twice.

import cnc
import numpy as np

##  PARAMETERS  ###############################################################
# Measuring rectangle at centre of groove!
# Zero will be the centre of the 'bottom' groove (portrait orientation).

# INPUTS:
longStraight = 435  #  a little over half the length of the long straight
shortStraight = 309.8
cornerRadius = 85
grooveWidth = 10
bitWidth = 3.175  #  not really sure about this
numParallelCuts = 4

# CALCULATED:
offsets = np.linspace(-1,1,numParallelCuts) * (grooveWidth - bitWidth) / 2
halfRectWidth = shortStraight / 2 + cornerRadius

##  CUT PATH  #################################################################
# longStraight = 2*longStraight - shortStraight  #  fix earlier blunder

cuts = [
  np.array(
    [ [-halfRectWidth-offset, longStraight+cornerRadius] ] +
    [ [-shortStraight/2 + (cornerRadius+offset) * np.cos(phi),
      cornerRadius + (cornerRadius+offset) * np.sin(phi)]
      for phi in np.linspace(np.pi, 1.5*np.pi, 91)] + 
    [ [shortStraight/2 + (cornerRadius+offset*0.8) * np.cos(phi),
      cornerRadius + (cornerRadius+offset) * np.sin(phi)]
      for phi in np.linspace(1.5*np.pi, 2*np.pi, 91)] + 
    [ [halfRectWidth+offset*0.8, longStraight+cornerRadius] ] )
  for offset in offsets]

# (x,y) = np.vstack((cuts[0], cuts[1][::-1])).T
(x,y) = np.vstack([cuts[k][::((-1)**k)] for (k,cut) in enumerate(cuts)]).T

# cnc.cutPath(x, y, 31, 1, 500, -25, "Henry.gcode")

cnc.cutPath(x, y, 35, 4, 750, -25, "Henry.gcode", 1.5)
cnc.cutPath(-x, y, 35, 4, 750, -25, "HenryWiderAtR.gcode", 1.5)

# START FROM 30mm ABOVE SURFACE!!!
