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
  sector = 2*np.pi/30
  for n in range(30):
    pt.addConvexPolygon(np.array([
      [0,0,-depth],
      [bkgRad*np.cos(n*sector),bkgRad*np.sin(n*sector),-depth],
      [bkgRad*np.cos((n+1)*sector),bkgRad*np.sin((n+1)*sector),-depth]]))

  print("Final shape:")
  print(pt)

  (offset, ballrad, speed) = (0.35, 1.0, 2000)
  # (offset, ballrad, speed) = (0.15, 0.75, 1200)

  pg = pt.toPG(offset)

  # WRITE THE GCODE:
  tp = pg.SingleToolNoOpt(ballrad)
  tp.PathToGCode(speed, f"{modName}.gcode")
