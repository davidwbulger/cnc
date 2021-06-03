# Load an stl file exported from Blender and manipulate it to produce code for CNC carving.

import cnc
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

##################################################################################################################
# PARAMETERS:

ozr = np.array([-5.4821,-2.1783])  #  old z range, in the stl file ("metres")
mzr = np.array([-2.4,-2.1783])     #  mid z range, for bas-relief compression ("metres")
nzr = np.array([-4, 0.5])         #  new z range, for the gcode (mm)

mxr = np.array([-0.32,0.5])
nxr = np.array([-50,50])

myr = np.array([-0.8,0.8])
nyr = np.array([5,5+(np.diff(nxr)*np.diff(myr)/np.diff(mxr))[0]])  #  to keep the xy aspect unchanged

multiplier = 1.21

##################################################################################################################
# LOAD AND MANIPULATE:

pt = cnc.polyTri("./BasRelief.stl")
print(pt)

oldz = pt.facets[:,:,2,None]
midz = mzr[0]+(mzr[1]-mzr[0])*(oldz-ozr[0])/(ozr[1]-ozr[0])
pt.facets *= midz/oldz

scjac = np.array([(nxr[1]-nxr[0])/(mxr[1]-mxr[0]),  #  scale Jacobian
  (nyr[1]-nyr[0])/(myr[1]-myr[0]), (nzr[1]-nzr[0])/(mzr[1]-mzr[0])])
scint = np.array([nxr[0],nyr[0],nzr[0]]) - np.array([mxr[0],myr[0],mzr[0]]) * scjac  #  scale intercept
pt.facets = multiplier * (scint[None,None,:] + scjac[None,None,:] * pt.facets)

pt = pt.afxform(np.array([[0,-1,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]]))

print("Final shape:")
print(pt)

print(-multiplier*nyr[::-1])
print(multiplier*nxr)
pg = pt.toPG(0.4,xrange=-multiplier*nyr[::-1],yrange=multiplier*nxr)

mx = min([np.min(xzp[0,:]) for xzp in pg.xz])
Mx = max([np.max(xzp[0,:]) for xzp in pg.xz])
my = np.min(pg.y)
My = np.max(pg.y)
mz = min([np.min(xzp[1,:]) for xzp in pg.xz])
Mz = max([np.max(xzp[1,:]) for xzp in pg.xz])
print(f"Cut box ({mx},{Mx})x({my},{My})x({mz,Mz}).")

##################################################################################################################
# PLOT:

if False:
  fig = plt.figure()
  ax = fig.add_subplot(projection="3d")
  pg.plot(ax)
  cnc.hackaspect(ax)
  mng = plt.get_current_fig_manager()
  mng.window.state("zoomed")
  plt.show()

##################################################################################################################
# WRITE THE GCODE:
if True:
  yinc=False
  # tp = pg.SingleToolNoOpt(1)
  # tp.PathToGCode(1800, "steamboat.gcode" if multiplier==1 else f"steamboat{multiplier}.gcode")
  tp = pg.SingleToolNoOpt(1.3, yinc=yinc)
  tp.PathToGCode(1800, "steam" + ("boat" if yinc else "back") + f"{multiplier}.gcode")
