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
bitWidth = 6  #  not really sure about this

# CALCULATED:
offsets = np.array([-1,1]) * (grooveWidth - bitWidth) / 2
halfRectWidth = shortStraight / 2 + cornerRadius

##  CUT PATH  #################################################################
cuts = [
  np.array(
    [ [-halfRectWidth-offset, shortStraight+cornerRadius] ] +
    [ [-shortStraight/2 + (cornerRadius+offset) * np.cos(phi),
      cornerRadius + (cornerRadius+offset) * np.sin(phi)]
      for phi in np.linspace(np.pi, 1.5*np.pi, 91)] + 
    [ [shortStraight/2 + (cornerRadius+offset) * np.cos(phi),
      cornerRadius + (cornerRadius+offset) * np.sin(phi)]
      for phi in np.linspace(1.5*np.pi, 2*np.pi, 91)] + 
    [ [halfRectWidth+offset, shortStraight+cornerRadius] ] )
  for offset in offsets]

(x,y) = np.vstack((cuts[0], cuts[1][::-1])).T

cnc.cutPath(x, y, 35, 4, 300, -25, "Henry.gcode", 1.3)
# START FROM 30mm ABOVE SURFACE!!!
