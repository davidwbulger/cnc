# Load an stl file designed for 3D printing (downloaded from the internet) and manipulate it to produce code for
# CNC carving.

import cnc
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

##################################################################################################################
# PARAMETERS:

##################################################################################################################
# LOAD AND MANIPULATE:

pt = cnc.polyTri("./Full_Bus.stl")  #  Full_Bus.stl")
print(pt)

# A = np.identity(4)
# A[:3,:3] = scipy.stats.special_ortho_group.rvs(3)

r = np.sqrt(0.5)
s = 0.1 ; t = 0.47  #  bas-relief shrinkage; overall shrinkage
theta = 0.1
A1 = np.array([[np.cos(theta),0,np.sin(theta),0], [0,1,0,0], [-np.sin(theta),0,np.cos(theta),0], [0,0,0,1]])  #  tilt down a tiny bit
A2 = np.array([[r,r,0,0],[0,0,1,0],[r,-r,0,0],[0,0,0,1]])  #  rotate to 45 degrees ("half-profile"?)
A3 = np.array([[t,0,0,0], [0,t,0,0], [0,0,s*t,0], [0,0,0,1]])  #  squash z direction for bas-relief
A = np.matmul(A3, np.matmul(A2,A1))

pt = pt.afxform(A)

bb = pt.bbox()
rad = 0.1 * np.min(np.diff(bb[:2,:],axis=1))
cartoucheVertices = np.row_stack([[bb[0,j]+rad*np.cos(theta),bb[1,k]+rad*np.sin(theta),bb[2,0]] for (j,k,th0) in
  [(1,1,0),(0,1,np.pi/2),(0,0,np.pi),(1,0,1.5*np.pi)] for theta in np.linspace(th0,th0+np.pi/2,31)])
pt.addConvexPolygon(cartoucheVertices)

pg = pt.toPG(2)

##################################################################################################################
# PLOT:

fig = plt.figure()
ax = fig.add_subplot(projection="3d")
pg.plot(ax)
cnc.hackaspect(ax)
mng = plt.get_current_fig_manager()
mng.window.state("zoomed")
# ax.set_xlim(-128, -121)
# ax.set_ylim(-452, -444)
# ax.set_zlim(43, 50)
plt.show()

##################################################################################################################
# WRITE THE GCODE:

