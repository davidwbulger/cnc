# create g-code for small dode

import cnc
import numpy as np
import matplotlib.pyplot as plt

# PARAMETERS
rf = 50  #  radius of face at outside
th = 9  #  thickness of board

ballrad = 1 - 0.05  #  [hack to leave glue gap]  #  radius of drillbit
offset = 2  #RESTORE  ballrad/4  #  gap between adjacent toolpaths
numTeeth = 3  #RESTORE  12  #  number of teeth per joined edge (on each piece)

edgePropJoined = 0.85  #  proportion of edge used for finger joint
mitreMargin = 1  #  vertically, for simplicity
numJointEdges = 3  #  this many edges will be cut for SMORF joints; the others will be plain mitre cuts. Varies by face.

# Calculated values:
edgeHalfLen = rf*np.sin(np.pi/5)
toothWidth = edgeHalfLen*edgePropJoined/numTeeth
cutsPerTooth = int(toothWidth/offset + 0.95)  #  ensure offset is at least nearly as small as it should
offset = toothWidth/cutsPerTooth  #  update offset for integer number of cuts per tooth
halfNumCuts = int((edgeHalfLen-offset/2)/offset + 0.95)
y = np.linspace((0.5-halfNumCuts)*offset, (halfNumCuts-0.5)*offset, 2*halfNumCuts)

ird = rf*np.sqrt((7+3*np.sqrt(5))/8)  # solid inradius to exterior face
irf = rf*np.cos(np.pi/5)  #  planar inradius to exterior edge
v = th*(1-0.5*np.sqrt(ird*ird+irf*irf)/ird)
v = 0.5*(th-mitreMargin)*(1-(irf/ird)**2)
xzmast = np.vstack((irf*(1-np.array([th+2*ballrad,th,mitreMargin,th,mitreMargin,0,-2*ballrad])/ird),
  np.array([th+2*ballrad,th,th-v,mitreMargin+v,mitreMargin,0,-2*ballrad])))  #  all nodes used on any path
xzflat = xzmast[:,[0,6]]        #  the straight paths (near the corners)
xzup = xzmast[:,[0,1,2,4,6]]    #  the tenons (or "teeth" or "fingers")
xzdown = xzmast[:,[0,1,3,4,6]]  #  the mortises
numCornerCuts = halfNumCuts - cutsPerTooth*numTeeth  #  number of cuts at each end, beyond the teeth
xz = [xzflat]*numCornerCuts + ([xzup]*cutsPerTooth + [xzdown]*cutsPerTooth)*numTeeth + [xzflat]*numCornerCuts

# Figure out how to "flip" the joint edge:
# [Note: this geometry is specific to this oblique joint, so it's not suitable for cnc.py.]
slopedir = xzmast[:,0]-xzmast[:,4]
slopedir = slopedir / np.linalg.norm(slopedir)
M = 2*np.outer(slopedir,slopedir)-np.identity(2)
moff = xzmast[:,4,None] - np.matmul(M,xzmast[:,4,None])

# Set up the loop to ensure the joint edge is round and meets:
rounded = cnc.PathGrid(y,xz)
print(str(rounded))
ctol=0.03
abserr=ctol+1
while abserr > ctol:
  flipt = cnc.PathGrid(rounded.y, [moff+np.matmul(M,tp) for tp in reversed(rounded.xz)])
  roundabove = rounded.roundJoint(ballrad,ctol)
  roundbelow = flipt.roundJoint(ballrad,ctol)
  abserr = (roundabove-roundbelow).maxabs()
  rounded = 0.5 * (roundabove + roundbelow)
  print(str(rounded))
cutpath = rounded.castToMold(ballrad, ctol, np.NINF)
print(str(cutpath))

tooth = cutpath.pacePathGrid(th+ballrad, th+ballrad+10, 3)
print(str(tooth))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d') 
rounded.plot(ax, 'blue')
cutpath.plot(ax, 'green')
tooth.plot(ax, 'red')
plt.show()
