# Load an stl file exported from Blender and manipulate it to produce code for CNC carving.

import cnc
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

##################################################################################################################
# PARAMETERS:

ozr = np.array([-5.4821,-2.1783])  #  old z range, in the stl file ("metres")
mzr = np.array([-2.4,-2.1783])     #  mid z range, for bas-relief compression ("metres")
nzr = np.array([-12, 0.8])         #  new z range, for the gcode (mm)

mxr = np.array([-0.32,0.5])
nxr = np.array([-50,50])

myr = np.array([-0.8,0.8])
nyr = np.array([15,15+(np.diff(nxr)*np.diff(myr)/np.diff(mxr))[0]])  #  to keep the xy aspect unchanged

multiplier = 1 # 0.25

##################################################################################################################
# LOAD AND MANIPULATE:

pt = cnc.polyTri("./BasRelief.stl")
print(pt)

# for i in range(pt.facets.size[0]):
#   for j in range(3):
#     newz = nzr[0]+(nzr[1]-nzr[0])*(pt.facets[i,j,2]-ozr[0])/(ozr[1]-ozr[0])
#     pt.facets[i,j,:] *= newz/pt.facets[i,j,2]

oldz = pt.facets[:,:,2,None]
midz = mzr[0]+(mzr[1]-mzr[0])*(oldz-ozr[0])/(ozr[1]-ozr[0])
pt.facets *= midz/oldz

scjac = np.array([(nxr[1]-nxr[0])/(mxr[1]-mxr[0]),  #  scale Jacobian
  (nyr[1]-nyr[0])/(myr[1]-myr[0]), (nzr[1]-nzr[0])/(mzr[1]-mzr[0])])
scint = np.array([nxr[0],nyr[0],nzr[0]]) - np.array([mxr[0],myr[0],mzr[0]]) * scjac  #  scale intercept
pt.facets = multiplier * (scint[None,None,:] + scjac[None,None,:] * pt.facets)

print("Final shape:")
print(pt)

pg = pt.toPG(0.5,xrange=multiplier*nxr,yrange=multiplier*nyr)

##################################################################################################################
# PLOT:

if True:
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
  tp = pg.SingleToolNoOpt(1)
  tp.PathToGCode(1000, "steamboat.gcode" if multiplier==1 else f"steamboat{multiplier}.gcode")
