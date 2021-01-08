# create g-code for clamps for dode

import cnc
import numpy as np

# This creates a gcode file to cut a clamping jig layer (presumably out of mdf).

# PARAMETERS:
cude = 4  #  safe cutting depth for a single pass
ballrad = 3.0
offset = 3.0
feedrate = 600
cutdepth = 55
rf = 67

# FIRSTLY, ROUGHLY CUT THE FACES:
phi = 0.5+np.sqrt(1.25)
pentag = 0.5*np.array([[phi,1-phi,-2,1-phi,phi],
  [np.sqrt(3-phi),np.sqrt(phi+2),0,-np.sqrt(phi+2),-np.sqrt(3-phi)],
  [0,0,0,0,0]]).T  #  facial planar exradius = 1
verts = np.vstack((rf*pentag-np.array([[0,0,cutdepth]]), (rf+cutdepth*(phi-1))*pentag))
verts += np.random.normal(scale=1e-6,size=verts.shape)

pt = cnc.polyTri(facets=np.empty((0,3,3)))
pt.addConvexPolygon(verts[:5,:])
pt.addConvexPolygon(verts[[0,1,6,5],:])
pt.addConvexPolygon(verts[[1,2,7,6],:])
pt.addConvexPolygon(verts[[2,3,8,7],:])
pt.addConvexPolygon(verts[[3,4,9,8],:])
pt.addConvexPolygon(verts[[4,0,5,9],:])

pg = pt.toPG(offset)
tp0 = pg.MultiToolGreedy(0,[{'bitrad':ballrad,'cude':cude,'ds':1}], yinc=True)[0]

# NOW CUT TO THE SLOPING EDGES (star pattern to avoid collisions):
nodes = np.zeros((3,1))  #  origin
nodes = np.hstack((nodes, (verts.T)[:,[2,7,4,9,1,6,3,8,0,5]]))
taxis = np.array([range(10),5*[0,1]])

# NOW ADD A FINER PASS:
shrinkage = 2/phi*ballrad*np.sin(0.5*np.arctan(2)) # reduce tooltip exradius to account for bit radius
for d in np.arange(0,cutdepth-5,2.0):
  radius = rf + (d-cutdepth)/phi - shrinkage
  nodes = np.hstack((nodes, (radius*pentag-np.array([[0,0,d]])).T))

for radius in np.arange(rf,0,-2.0):
  nodes = np.hstack((nodes, (radius*pentag-np.array([[0,0,cutdepth]])).T))

nodes = np.hstack((nodes, np.array([[0,0],[0,0],[-cutdepth,0]])))
taxis = np.hstack((taxis,np.array([[nodes.shape[1]-2],[0]])))
tp1 = cnc.ToolPath(nodes, taxis)
tp = cnc.catToolPaths([tp0,tp1])
tp.PathToGCode(600,"clampJig.gcode")

##################################################################################################################
#####     APPENDIX A:  VERTICES' COORDINATES     #################################################################
##################################################################################################################
"""
This script creates gcode files for cutting jigs for a dodecahedral box. The size parameter we're using, rf,
is defined as the planar exradius of the exterior faces, i.e., the distance from the centre of a face to a vertex
of that face. To obtain coordinates of the exterior vertices, multiply the following vectors by rf/2, and note
that phi denotes the "golden ratio" (1+sqrt(5))/2:

BOTTOM FACE:
          phi                 -phi+1              -2                  -phi+1              phi
          sqrt(3-phi)         sqrt(phi+2)         0                   -sqrt(phi+2)        -sqrt(3-phi)
          -phi-1              -phi-1              -phi-1              -phi-1              -phi-1

'TROPIC OF CAPRICORN':
          phi+1               -1                  -2*phi              -1                  phi+1
          sqrt(phi+2)         sqrt(4*phi+3)       0                   -sqrt(4*phi+3)      -sqrt(phi+2)
          -phi+1              -phi+1              -phi+1              -phi+1              -phi+1

'TROPIC OF CANCER':
2*phi               1                   -phi-1              -phi-1              1
0                   sqrt(4*phi+3)       sqrt(phi+2)         -sqrt(phi+2)        -sqrt(4*phi+3)
phi-1               phi-1               phi-1               phi-1               phi-1

TOP FACE:
2                   phi-1               -phi                -phi                phi-1
0                   sqrt(phi+2)         sqrt(3-phi)         -sqrt(3-phi)        -sqrt(phi+2)
phi+1               phi+1               phi+1               phi+1               phi+1
"""
