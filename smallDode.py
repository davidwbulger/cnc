# create g-code for small dode

import cnc
import numpy as np
import matplotlib.pyplot as plt

# This will create six separate gcode files:
#   Blank.gcode: cuts a blank, using a wider drill bit
#   Guide.gcode: cuts a shallow line for positioning the blank for the next four programs
#   SDFace2.gcode to SDFace5.gcode: cuts the edges onto a blank
# Note that the four programs SDFace#.gcode assume a starting position that's 20mm above the centre of the bottom
# of the pentagonal panel (so, 11mm above the top of the panel, if it's exactly 9mm thick).

##################################################################################################################
#### ##    PARAMETERS:    ## #####
##################################################################################################################

# Blank stuff (not really used, in the end):
blankCutDepth = 6 # 14  #  board seems to be 10mm, and bit radius is 3mm; this is 1mm extra
blankNumPasses = 3 # 3
blankbitrad = 1 # 3
excess = 0 # 2  #  excess, in mm, around the margin of the pentagon blanks cut in rough phase

# Main programs (SDFace#.gcode) stuff:
rf = 67  #  planar exradius of face at outside [see Appendix A]
th = 9  #  thickness of board
originHeight = 20
ballrad = 1.0
gcodeToolSeq = [{'bitrad':ballrad,'cude':1.5,'ds':1}]

# Parameters for tenon & mortise shape:
edgePropJoined = 0.9  #  max proportion of edge used for finger joint
mitreMargin = 1  #  vertically, for simplicity
cutsPerTooth = 6
if cutsPerTooth%4 != 2: raise ValueError("cutsPerTooth must be twice an odd number.")
mortiseWidth = 2.35  #  This should be 2ce the bit radius, but in practice the bit carves out a little more.
toothWidth = 2.15
# Make the tooth "shorter" (in the direction normal to the "joint plane" described by xzflat) by this amount:
brachidont = 0.4
if (1-2/cutsPerTooth) * (mortiseWidth+toothWidth) / 2 >= 2*ballrad:
  raise ValueError("These values might give multiple cuts within the mortise.")
toothRoundingFactor = 1.9  #  I don't know why 1.0 doesn't work, but this seems better.

boolPlot = False

## TOOTH WIDTH/SHAPE STUFF:

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

# Calculated values:
offset = (mortiseWidth + toothWidth) / cutsPerTooth
edgeHalfLen = rf*np.sin(np.pi/5)  #  exact at exterior (bottom) face
numTeeth = int(np.floor(2*edgeHalfLen*edgePropJoined/(mortiseWidth+toothWidth)))
v = 0.5*(th-mitreMargin)*(1+(irf/ird)**2)  #  length of vertical (i.e., face-normal) edge of mortise
rho = np.linspace(2*np.arctan(irf/ird), 0, 8)

xzflat = np.vstack((irf*(1-np.array([th+2*ballrad,-2*ballrad])/ird), np.array([2*ballrad,-th-2*ballrad])))
xzdown = np.vstack((irf*(1-np.array([th+2*ballrad,th,th,mitreMargin,-2*ballrad])/ird),
  np.array([2*ballrad,0,-v,-th+mitreMargin,-th-2*ballrad])))
xzup = [None] * (cutsPerTooth//2)
for k in range(cutsPerTooth//2):
  # The following variable round's the tooth's edge by reducing its height away from its central ridge:
  extraShortening = mortiseWidth/2 - np.sqrt((mortiseWidth/2)**2-(offset*(k-(cutsPerTooth//4)))**2)
  # Determine the centre of the arc at the tooth's apex:
  brachshift = ((brachidont+extraShortening)/np.linalg.norm([ird,irf])) * np.array([ird,irf])
  ctc = (np.array([irf*(1-mitreMargin/ird), -th+mitreMargin+v])  #  mirror image of nadir of mortise
    - toothRoundingFactor*ballrad * np.array([1, irf/ird])  #  centre of circ tangt to mortise reflection's ridge
    - brachshift)  #  shorten tooth for fit tolerance
  toothBaseShrink = brachshift[0] * np.array([-1,ird/irf])  #  shift of tenon's base vertex relative to mortise's
  xzup[k] = np.column_stack((xzdown[:,0], xzdown[:,1]-toothBaseShrink,
    ctc[:,None] + toothRoundingFactor*ballrad*np.vstack((np.cos(rho), np.sin(rho))),
    xzdown[:,3]+toothBaseShrink, xzdown[:,4]))
  xzup[k][0,-3] = min(xzup[k][0,-3:])  #  since otherwise numerical error could spoil non-decreasingness

numCornerCuts = int(np.ceil(edgeHalfLen/offset-(cutsPerTooth//2)*numTeeth))
xz = ([xzflat] * numCornerCuts + (xzup + [xzdown] * (cutsPerTooth//2)) * numTeeth + [xzflat] * numCornerCuts)
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
flattooth = (flat+th-originHeight+0.8*gcodeToolSeq[0]['cude']
  ).MultiToolGreedy(th-originHeight, gcodeToolSeq, yinc=True)[0]
# flattooth = (flat+th-originHeight).MultiToolGreedy(th-originHeight, gcodeToolSeq, yinc=True)[0]

# This edit is intended to cut the final pass of the flat edges with a finer offset, to get them smoother. Using
# the finer offset all along is unnecessarily slow.
s = np.linspace(0, 1, 1*int(np.ceil((xzflat[0,1]-xzflat[0,0])*(3/offset))))
flat = cnc.PathGrid(xzflat[0,0]+(xzflat[0,1]-xzflat[0,0])*s,
  [np.array([[-edgeHalfLen,edgeHalfLen],2*[xzflat[1,0]+(xzflat[1,1]-xzflat[1,0])*sval]]) for sval in s])

finerToolSeq = [{k:(40 if k=='cude' else v) for (k,v) in gcodeToolSeq[0].items()}]
finerflattooth = (flat+th-originHeight
  ).MultiToolGreedy(th-originHeight, finerToolSeq, yinc=True)[0]
flattooth = cnc.catToolPaths([flattooth, finerflattooth])

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
