# create g-code to cut an x-y path to a given depth

import cnc
import sys
import numpy as np

if len(sys.argv) < 9 or len(sys.argv)%2 < 1:
  print("Usage: python cutPath.py finalDepth numPasses feedRate safeZ "
    "x0 y0 ... xK yK")
else:
  params = tuple(float(a) for a in sys.argv[1:])
  (finalDepth, numPasses, feedRate, safeZ) = params[:4]
  numPasses = int(numPasses)
  xy = np.reshape(params[4:],(-1,2))

  cnc.cutPath(xy[:,0], xy[:,1], finalDepth, numPasses, feedRate, safeZ,
    "cutpath.gcode")
