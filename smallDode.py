# create g-code for small dode

import cnc
import numpy as np
import matplotlib.pyplot as plt

# This will create five separate gcode files:
#   Blank.gcode: cuts a blank, using a wider drill bit


##################################################################################################################
#### ##    PARAMETERS:    ## #####
##################################################################################################################

blankCutDepth = 6 # 14  #  board seems to be 10mm, and bit radius is 3mm; this is 1mm extra
blankNumPasses = 3 # 3
blankbitrad = 1 # 3
excess = 0 # 2  #  excess, in mm, around the margin of the pentagon blanks cut in rough phase
rf = 57 # 67  #  planar exradius of face at outside [see Appendix A]
th = 9  #  thickness of board
originHeight = 20
ballrad = 1.0
gcodeToolSeq = [{'bitrad':ballrad,'cude':1.5,'ds':1}]

boolPlot = False

## TOOTH WIDTH/SHAPE STUFF:

# After a failed first try (I believe due to numerical precision errors) I'm now fixing the tooth width as the
# width of the drill bit---essentially. The true bit width will be the mortisse width, and the tooth width will be
# slightly smaller, with the difference determined by glueGapHack.
glueGapHack = 0.001 # 0.03
realBallRad = 1.17  #  it's REALLY 1.0, but this seems to be the half-width of the cut, between teeth at least

# realBallRad = 1.00
# glueGapHack = 0.06 # 0.03
# mortiseWidthHack = 1.19 # 0.17
rf = 24

#### ##    CALCULATED VALUES:    ## #####

ird = rf*np.sqrt((7+3*np.sqrt(5))/8)  # solid inradius to exterior face
irf = rf*np.cos(np.pi/5)  #  planar inradius to exterior edge
phi = (1+np.sqrt(5))/2

##################################################################################################################
#### ##    CREATE THE GCODE TO CUT OUT THE BLANKS:    ## #####
##################################################################################################################

# This would be used 12 times to cut pentagonal blanks that are a bit too big.

verts = 0.5*(rf+blankbitrad+excess)*np.array([[phi, np.sqrt(3-phi)], [1-phi,np.sqrt(phi+2)], [-2,0],
  [1-phi,-np.sqrt(phi+2)], [phi,-np.sqrt(3-phi)]])
bnodes = np.zeros((3,7+6*blankNumPasses))
bnodes[:,:4] = np.array([[0,0,0],[0,0,5],[*verts[0,:],5],[*verts[0,:],0]]).T
bnodes[:,-3:] = np.array([[*verts[0,:],5], [0,0,5], [0,0,0]]).T 
bnodes[:2,4] = verts[0,:]
bnodes[:2, 5:10] = verts.T
bnodes[2,5:10] = -blankCutDepth/blankNumPasses
for k in range(10,(4+6*blankNumPasses)):
  bnodes[:,k] = bnodes[:,(k-6)] - np.array([0,0,blankCutDepth/blankNumPasses])

btaxis = np.array([[0, 3, 4+6*blankNumPasses], [0,1,0]])

btp = cnc.ToolPath(bnodes, btaxis)
btp.PathToGCode(500, "Blank.gcode")

#### ##    CREATE THE GCODE TO MARK THE SPOILBOARD WITH A PENTAGON FOR POSITIONING:    ## #####

# After cutting the blanks, secure a big piece of spoiler board to the base board and then carve this shallow
# "guide" into it to help position the first blank. The first joint cut will obliterate the guide, but it will
# leave an equally effective visual indication of where to put the remaining blanks.

# Remember:
#   don't change the XY position from cutting the guide to the end of cutting ALL blanks' edges;
#   manually select grain orientation for 2JE, 3JE & 4JE blanks (joint edges are anticlockwise from right edge).

verts = 0.5*(rf+1.3+excess)*np.array([[phi, np.sqrt(3-phi)], [1-phi,np.sqrt(phi+2)], [-2,0],
  [1-phi,-np.sqrt(phi+2)], [phi,-np.sqrt(3-phi)]])
bnodes = np.zeros((3,7+6*1))
bnodes[:,:4] = np.array([[0,0,0],[0,0,5],[*verts[-1,:],5],[*verts[-1,:],0]]).T
bnodes[:,-3:] = np.array([[*verts[-1,:],5], [0,0,5], [0,0,0]]).T 
bnodes[:2,4] = verts[-1,:]
bnodes[:2, 5:10] = verts.T
bnodes[2,5:10] = -0.6

btaxis = np.array([[0, 3, 10], [0,1,0]])

btp = cnc.ToolPath(bnodes, btaxis)
btp.PathToGCode(300, "Guide.gcode")

##################################################################################################################
#### ##    CREATE THE TOOLPATH FOR THE EDGES WITH TEETH:    ## #####
##################################################################################################################

# ## Tooth width/shape stuff:
# # We'll tell the cnc module that this is the bit radius, to persuade it to leave room for glue:
# ballrad = realBallRad-glueGapHack
# xzmast = np.vstack((irf*(1-np.array([th+2*ballrad,th,mitreMargin,th,mitreMargin,0,-2*ballrad])/ird),
#   np.array([2*ballrad,0,-v,-th+mitreMargin+v,-th+mitreMargin,-th,-th-2*ballrad])))  # all nodes used on any path
# gcodeToolSeq = [{'bitrad':ballrad,'cude':1.5,'ds':1}]
# toothWidth= mortiseWidthHack*(2*realBallRad-1.9*glueGapHack) # 2*(rBR-gGH) would be "perfect"; this leaves a tol
# edgeHalfLen = rf*np.sin(np.pi/5)  #  exact at exterior (bottom) face
# numTeeth = int(np.floor(edgeHalfLen*edgePropJoined/toothWidth))
# 
# # Update offset for integer number of cuts per tooth (ACTUAL offset, not hacked for mortise width):
# offset = 2*toothWidth/cutsPerTooth
# xzflat = xzmast[:,[0,6]]        #  the straight paths (near the corners)
# xzup = xzmast[:,[0,1,2,4,6]]    #  the tenons (or "teeth" or "fingers")
# xzdown = xzmast[:,[0,1,3,4,6]]  #  the mortises
# 
# numCornerCuts = int(np.ceil(edgeHalfLen/offset-(cutsPerTooth//2)*numTeeth))
# xz= [xzflat]*numCornerCuts+([xzup]*(cutsPerTooth//2)+[xzdown]*(cutsPerTooth//2))*numTeeth+[xzflat]*numCornerCuts
# maxAbsY = offset*0.5*(len(xz)-1) / mortiseWidthHack
# y = np.linspace(-maxAbsY, maxAbsY, len(xz))
# 
# # Figure out how to "flip" the joint edge:
# # [Note: this geometry is specific to this oblique joint, so it's not suitable for cnc.py.]
# slopedir = xzmast[:,0]-xzmast[:,4]
# slopedir = slopedir / np.linalg.norm(slopedir)
# M = 2*np.outer(slopedir,slopedir)-np.identity(2)
# moff = xzmast[:,4,None] - np.matmul(M,xzmast[:,4,None])
# 
# # Set up the loop to ensure the joint edge is round and meets:
# rounded = cnc.PathGrid(y,xz)
# ctol=0.03
# abserr=ctol+1
# while abserr > ctol:
#   flipt = cnc.PathGrid(rounded.y, [moff+np.matmul(M,tp) for tp in reversed(rounded.xz)])
#   roundabove = rounded.roundJoint(ballrad,ctol)
#   roundbelow = flipt.roundJoint(ballrad,ctol)
#   # roundabove = rounded.castToMold(ballrad,ctol).moldToCast(ballrad,ctol)
#   # roundbelow = flipt.castToMold(ballrad,ctol).moldToCast(ballrad,ctol)
#   abserr = (roundabove-roundbelow).maxabs()
#   rounded = 0.5 * (roundabove + roundbelow)
#   # rounded = roundabove.min(roundbelow)  #  elementwise min of two PathGrids
# breakpoint()
# 
# teethtooth = (rounded+th-originHeight).MultiToolGreedy(th-originHeight, gcodeToolSeq, yinc=True)[0]
# teethtooth.nodes[1,:] *= mortiseWidthHack

# Parameters for tenon & mortise shape:
edgePropJoined = 0.5 # 0.9  #  max proportion of edge used for finger joint
mitreMargin = 1  #  vertically, for simplicity
cutsPerTooth = 6
if cutsPerTooth%4 != 2: raise ValueError("cutsPerTooth must be twice an odd number.")
mortiseWidth = 2.35  #  This should be 2ce the bit radius, but in practice the bit carves out a little more.
toothWidth = 2.15
offset = (mortiseWidth + toothWidth) / cutsPerTooth
# Make the tooth "shorter" (in the direction normal to the "joint plane" described by xzflat) by this amount:
brachidont = 0.3

if (1-2/cutsPerTooth) * (mortiseWidth+toothWidth) / 2 >= 2*ballrad:
  raise ValueError("These values might give multiple cuts within the mortise.")

# Calculated values:
edgeHalfLen = rf*np.sin(np.pi/5)  #  exact at exterior (bottom) face
numTeeth = int(np.floor(2*edgeHalfLen*edgePropJoined/(mortiseWidth+toothWidth)))
v = 0.5*(th-mitreMargin)*(1+(irf/ird)**2)  #  length of vertical (i.e., face-normal) edge of mortise
# Determine the centre of the arc at the tooth's apex:
brachshift = (brachidont/np.linalg.norm([ird,irf])) * np.array([ird,irf])
ctc = (np.array([irf*(1-mitreMargin/ird), -th+mitreMargin+v])  #  mirror image of nadir of mortise
  - ballrad * np.array([1, irf/ird])  #  centre of circle tangent to edges of mortise's mirror image
  - brachshift)  #  shorten tooth for fit tolerance
toothBaseShrink = brachshift[0] * np.array([-1,ird/irf])  #  shift of tenon's base vertex relative to mortise's
rho = np.linspace(2*np.arctan(irf/ird), 0, 8)

xzflat = np.vstack((irf*(1-np.array([th+2*ballrad,-2*ballrad])/ird), np.array([2*ballrad,-th-2*ballrad])))
xzdown = np.vstack((irf*(1-np.array([th+2*ballrad,th,th,mitreMargin,-2*ballrad])/ird),
  np.array([2*ballrad,0,-v,-th+mitreMargin,-th-2*ballrad])))
xzup = np.column_stack((xzdown[:,0], xzdown[:,1]-toothBaseShrink,
  ctc[:,None] + ballrad*np.vstack((np.cos(rho), np.sin(rho))),
  xzdown[:,3]+toothBaseShrink, xzdown[:,4]))
xzup[0,-3] = min(xzup[0,-3:])  #  since otherwise numerical error could spoil non-decreasingness

numCornerCuts = int(np.ceil(edgeHalfLen/offset-(cutsPerTooth//2)*numTeeth))
xz= [xzflat]*numCornerCuts+([xzup]*(cutsPerTooth//2)+[xzdown]*(cutsPerTooth//2))*numTeeth+[xzflat]*numCornerCuts
xz = ([xzflat]*numCornerCuts
  + ([xzflat]*((cutsPerTooth-2)//4) +[xzup] +[xzflat]*((cutsPerTooth-2)//4) +[xzdown]*(cutsPerTooth//2)) *numTeeth
  + [xzflat]*numCornerCuts)
maxAbsY = offset*0.5*(len(xz)-1)
y = np.linspace(-maxAbsY, maxAbsY, len(xz))
teethtooth = (cnc.PathGrid(y,xz)+(th-originHeight)).MultiToolGreedy(th-originHeight, gcodeToolSeq, yinc=True)[0]

#### ##    CREATE THE TOOLPATH FOR THE EDGES WITHOUT TEETH:    ## #####

# Actual vertices are [xzmast[0,{0,6}], {+,-}edgeHalfLen, xzmast[1,{0,6}].
# We want the cuts in the long direction, though, so we'll construct this at 90deg & then rotate it back.
# Thus transformed vertices are [{+,-}edgeHalfLen, xzmast[0,{0,6}], xzmast[1,{0,6}].
s = np.linspace(0, 1, 1*int(np.ceil((xzflat[0,1]-xzflat[0,0])/offset)))

flat = cnc.PathGrid(xzflat[0,0]+(xzflat[0,1]-xzflat[0,0])*s,
  [np.array([[-edgeHalfLen,edgeHalfLen],2*[xzflat[1,0]+(xzflat[1,1]-xzflat[1,0])*sval]]) for sval in s])
flattooth = (flat+th-originHeight).MultiToolGreedy(th-originHeight, gcodeToolSeq, yinc=True)[0]
flattooth = flattooth.afxform(np.array([[0,1,0,0],[-1,0,0,0],[0,0,1,0],[0,0,0,1]]))


##################################################################################################################
#### ##    NOW COMBINE FIVE EDGES TO CREATE THE PIECE:    ## #####
##################################################################################################################

theta = 0.4*np.pi
c = np.cos(theta)
s = np.sin(theta)
R = np.array([[c,-s,0,0],[s,c,0,0],[0,0,1,0],[0,0,0,1]])

for numJointEdges in [2,3,4,5]:
  alltooth = cnc.catToolPaths(list(teethtooth.afxform(np.linalg.matrix_power(R,e)) for e in range(numJointEdges))
    + list(flattooth.afxform(np.linalg.matrix_power(R,e)) for e in range(numJointEdges,5)))
  if boolPlot:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d') 
    rounded.plot(ax, 'blue')
    # cutpath.plot(ax, 'green')  #  no longer works since I overwrite cutpath in the flat edge section.
    alltooth.plot(ax, 'red', linewidth=1)
    cnc.hackaspect(ax)
    plt.show()
  alltooth.PathToGCode(300, f"SDFace{numJointEdges}.gcode")
  if numJointEdges==2:
    oops = cnc.ToolPath(alltooth.nodes * np.array([[1],[-1],[1]]), alltooth.taxis)
    oops.PathToGCode(300, "Oops.gcode")

# Actually rotate these into position to look at the geometry of the equatorial clamping jig.

##################################################################################################################
#####     APPENDIX A:  VERTICES' COORDINATES     #################################################################
##################################################################################################################
"""
This script creates gcode files for cutting the faces of a dodecahedral box. The size parameter we're using, rf,
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
