# create g-code to level a rectangular region

import cnc
import sys
import numpy as np

if len(sys.argv) != 9:
  print("Usage: python level.py lox hix loy hiy z safez offset feedrate")
else:
  xran = np.array([float(a) for a in sys.argv[1:3]])
  yran = np.array([float(a) for a in sys.argv[3:5]])
  zht = float(sys.argv[5])
  sht = float(sys.argv[6])
  offset = float(sys.argv[7])
  feedrate = float(sys.argv[8])
  cnc.thicknesser(xran, yran, zht, sht, offset, feedrate, "level.gcode")
