# Load an stl file, transform it into relief and save gcode to cut multiple
# passes.

import cnc
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import sys

if len(sys.argv) != 9:
  print("Usage: python relief.py modName height depth bkgRad stckHt bitRad" +
    " offset numPasses")
  print("where")
  print(" modName    is the name, sans extension, of stl model;")
  print(" height     is the required horizontal extent of the carving;")
  print(" depth      is the carving's vertical depth (flat bkgd to proudest);")
  print(" bkgRad     is the radius of the flat circular background;")
  print(" stckHt     is the height of the orig stock wrt the bit's home pos;")
  print(" bitRad     is the radius of the ball-nosed butting bit;")
  print(" offset     is the lateral offset between toolpaths;")
  print(" numPasses  is the number of passes to do.")
  print("All measurements are in mm.")
else:
  modName = sys.argv[1]
  (height, depth, bkgRad, stckHt, bitRad, offset) = tuple(
    float(a) for a in sys.argv[2:8])
  numPasses = int(sys.argv[8])
  if stckHt >= 0:
    raise ValueError("The stckHt parameter should be negative.")
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

  # # Also move down to stock height:
  # pt = pt.afxform(np.array([[1,0,0,0],[0,1,0,0],[0,0,1,stckHt],[0,0,0,1]]))

  print("Final shape:")
  print(pt)

  pg = pt.toPG(offset)

  # WRITE THE GCODE:

  # Not this way, because we want multiple layers:
  # tp = pg.SingleToolNoOpt(bitRad)

  # Not this way, because something's awry with MultiToolGreedy:
  # cude = depth/(numPasses-0.3)  #  make final pass a little shallower
  # tp = pg.MultiToolGreedy(stckHt,[{'bitrad':bitRad,'cude':cude,'ds':1}],
  #   ymono=True)[0]

  tpOnePass = pg.SingleToolNoOpt(bitRad)
  cude = np.concatenate(
    (depth - np.arange(1,numPasses)*depth/(numPasses-0.3),[0]))
  tpAll = cnc.catToolPaths(list(tpOnePass.afxform(np.array(
    [[1,0,0,0],[0,1,0,0],[0,0,1,cd],[0,0,0,1]])) for cd in cude))
  # Also move down to stock height:
  tp = tpAll.afxform(np.array([[1,0,0,0],[0,1,0,0],[0,0,1,stckHt],[0,0,0,1]]))

  tp.PathToGCode(1000, f"{modName}_{bitRad}mm.gcode")
