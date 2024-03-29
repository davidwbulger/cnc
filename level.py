# create g-code to level a rectangular region

import cnc
import sys
import numpy as np

if len(sys.argv) not in [9, 10, 12]:
  print("Usage: python level.py lox hix loy hiy z safez offset feedrate "
    "[quick [numPasses extraDepthPerPass]]")
else:
  xran = np.array([float(a) for a in sys.argv[1:3]])
  yran = np.array([float(a) for a in sys.argv[3:5]])
  zht = float(sys.argv[5])
  sht = float(sys.argv[6])
  offset = float(sys.argv[7])
  feedrate = float(sys.argv[8])
  quick = False if len(sys.argv)<10 else bool(sys.argv[9])
  numPasses = 1 if len(sys.argv)<12 else float(sys.argv[10])
  extraDepthPerPass = 1 if len(sys.argv)<12 else float(sys.argv[11])
  if numPasses>1:
    zht = zht - extraDepthPerPass * np.arange(numPasses)
  cnc.thicknesser(xran, yran, zht, sht, offset, feedrate, "level.gcode", quick)
