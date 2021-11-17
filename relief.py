# Load an stl file, transform it into relief and save gcode.

import cnc
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import sys

if len(sys.argv) != 5:
  print("Usage: python relief.py modName height depth bkgRad")
else:
  modName = sys.argv[1]
  (height, depth, bkgRad) = tuple(float(a) for a in sys.argv[2:])
  pt = cnc.polyTri(f"./{modName}.stl")
  print(pt)

  # Assume raw model height & depth are both 1. Scale as required:
  pt = pt.afxform(np.array(
    [[height,0,0,0],[0,height,0,0],[0,0,depth,0],[0,0,0,1]]))

  # Also add a floor at (0,0,-depth) of radius bkgRad:
  for n in range(30):
    pt.addConvexPolygon(np.array([
      [0,0,-depth],
      [bkgRad*np.cos(n/30),bkgRad*np.sin(n/30),-depth],
      [bkgRad*np.cos((n+1)/30),bkgRad*np.sin((n+1)/30),-depth]]))

  print("Final shape:")
  print(pt)

  pg = pt.toPG(0.3)

  # WRITE THE GCODE:
  tp = pg.SingleToolNoOpt(1)
  tp.PathToGCode(1000, f"{modName}.gcode")
