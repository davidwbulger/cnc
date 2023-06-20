# Project with Nikau to cut blocks for printing. Based on BibleCovers.py.
# Started 2023 June 15.

"""
Ideas:
  Remember to reverse the image, in two ways:
    flip it left to right
    cut where ink SHOULDN'T appear
  Function to cut region:
    find (outer?) boundary
    cut along it (adjusting cut depth as necessary)
    compute remainder (with cut path overlap)
    partition remainder into regions
    order regions by proximity (now?)
    recursively call same function on each region
    concerns:
      what about non-simply-connected regions?
      linear, cyclic & hexamorphic paths
  Find previous images (mousepad & bible covers) & determine how different they
    really are. Maybe not much adaptation is needed. Actually, I've looked at
    mousepad.jpg, & it IS very much a line drawing, better suited to
    skeletonisation. We'll definitely need to adapt the method a bit here.
  Remind myself: is there an easy way to visualise the determined path?

17 June:
Simply cut a raster path? It would leave zigzag edges for edges perpendicular
to the grooves. Could do two passes, in x & y directions, but diagonal edges
would still zigzag.
How about this:
  Compute boundary
  Cut paths along all boundaries
  Determine what other cutting is necessary
  Do remaining cut by raster
A bit more elaborate, but will give a much better result.

For cutting the paths along the boundaries, one challenge is to determine the
depth and lateral offset at each point along the path, since there may for
instance be sharp corners or tiny regions where we should cut at reduced depth.
A simple and fairly quick method would be to discretise the boundary and, at
each point, use bisection to find the largest tangent circle (up to bit
diameter) not intersecting the inked region.
"""

import cnc
import numpy as np
from scipy import ndimage
from scipy.signal import decimate
from skimage.filters import gaussian
from skimage import io
from skimage import measure
from skimage.morphology import skeletonize
from matplotlib.path import Path
import matplotlib.pyplot as plt

# PARAMETERS:
imfile = "cicada.jpg"
outputWidth = 72  #  width of desired output, in millimetres
outfile = "printBlock.gcode"  #  widest point is about 72 pixels
downSampleRate = 8
blurRad = 1.2  #  for antialiasing the boundary cut paths, if dsRate==1

vBitAngle = 90  #  bit angle in degrees
maxDepth = 3.2  #  3.2 mm
sanddepth = 0.2  #  cut this much deeper, so a groove remains after sanding 

# CALCULATED FROM PARAMETERS:
# Note: vBitAngle is the dihedral angle between the two sides of a groove cut
# by the bit, i.e., TWICE the angle between a groove wall (or bit edge) and the
# vertical.
# bit-angle multiplier; ratio of depth of cut to halfwidth of cut:
bam = 1/np.tan(vBitAngle*np.pi/360)
maxRad = maxDepth / bam
image = io.imread(imfile, as_gray=True).astype(float)[:,::-1]  #  flip LR
if downSampleRate > 1:
  image = decimate(decimate(image, downSampleRate, axis=0),
    downSampleRate, axis=1)
else:
  image = gaussian(image,sigma=blurRad)
dpmm = image.shape[1]/outputWidth  #  dots per millimetre

# FIRST CALCULATE THE BOUNDARY CUTS, FOR SMOOTHER EDGES:
contours = measure.find_contours(image, 0.5)
# Note, by default, the contours are positively oriented around low-valued,
# i.e., inked regions.
# Precalculate some stuff:
# Normals pointing into the cut (non-ink) region (note that in contours, the
# loop is explicit [the last vertex is repeated at the start] whereas this
# doesn't happen in normals):
normals = [
  vertNormal/np.linalg.norm(vertNormal,axis=1)[:,None] for vertNormal in (
    edgeNormal+np.roll(edgeNormal,-1,axis=0) for edgeNormal in (
      edgeParallel/np.linalg.norm(edgeParallel,axis=1)[:,None] @
      np.array([[0,-1],[1,0]]) for edgeParallel in (
        vertex[1:]-vertex[:-1] for vertex in contours)))]
if (checkContours := False):
  k = np.argmax([l>10 and l<20 for l in (n.shape[0] for n in normals)])
  print(len(contours))
  print(np.array([[f([c.shape[z] for c in contours])
    for f in [np.min, np.mean, np.max]] for z in range(2)]))
  # Display the image and plot all contours found
  (fig, ax) = plt.subplots()
  ax.imshow(image, cmap=plt.cm.gray)
  for contour in contours:
    ax.plot(contour[:, 1], contour[:, 0], linewidth=1, color='r')
  for (v,n) in zip(contours[k][1:], normals[k]):
    ax.plot([v[1],v[1]+n[1]], [v[0],v[0]+n[0]], linewidth=1, color='b')
  ax.axis('image')
  ax.set_xticks([])
  ax.set_yticks([])
  plt.show()
  exit()

"""
Alright, so I'm thinking if we precalculate the quadratic terms, it should
speed up checking which points are in which circles. In particular we're gonna
be interested in whether (u,v) is inside the circle of radius t that's centred
on (x,y)+(p,q), where (x,y) is a vertex on the contour & (p,q) is the normal.
It's inside if
t^2 > [(x+tp)-u]^2+[(y+tq)-v]^2
 = x^2+2txp+t^2p^2-2xu-2tpu+u^2 + y^2+2tyq+t^2q^2-2yv-2tqv+v^2 
but note that (p,q) is normalised, so t^2=t^2(p^2+q^2), so (u,v) is inside if
0 > x^2+2txp-2xu-2tpu+u^2 + y^2+2tyq-2yv-2tqv+v^2, i.e., if
[u,u^2,v,v^2] @ ([2x,-1,2y,-1] + t[2p,0,2q,0]) > x^2+y^2 + 2t(xp+yq).
Can simply try with max radius first, and thereby identify any at risk pixels;
then calculate reduced radius for each pixel by solving this linear eqn, then
use the lowest.
Further manipulation:
[u,u^2,v,v^2] @ ([2x,-1,2y,-1] + t[2p,0,2q,0]) > x^2+y^2 + 2t(xp+yq) <==>
x^2+y^2 - U@[2x,-1,2y,-1] < 2t[(u-x)p+(v-y)q]  <==>
(x-u)^2+(y-v)^2 < 2t[(u-x).n]  <==>
t > (u-x).(u-x) / 2(u-x).n
Are we tacitly assuming a 90 degree bit?
Not so far; this "t" (max (u-x).(u-x)/2(u-x).n) is the biggest radius circle
tangent to the contour that doesn't overlap the inked region.
"""
X = np.array([np.arange(0.5,image.shape[0])/dpmm]*image.shape[1])
Y = np.array([np.arange(image.shape[1]-0.5,0,-1)/dpmm]*image.shape[0]).T
W = np.array([X[image<0.5],Y[image<0.5]]).T  #  all inked pixels as rows
boundaryCuts = [bocut(W,c,n) for (c,n) in zip(contours,normals)]

def bocut(W,cmat,nmat):
  # Determine boundary cut given contour cmat & normals nmat.
  t = maxRad*np.ones((cmat.shape[0],1))
  for (k,(c,n)) in enumerate(zip(cmat,nmat)):
    Wx = W-c
    obx = np.sum(Wx*Wx,axis=1) < 2*maxRad*Wx@n
    if np.any(obx):
      Wo = Wx[obx]
      t[k] = np.min(np.sum(Wo*Wo,axis=1) / (2*Wo@n))
  # Path of bit vertex:
  return np.concatenate((cmat+t*nmat,-np.tan(vBitAngle/2)*t))

exit()

# EXPERIMENT:
copa = Path.make_compound_path(*[Path(vertices=c) for c in contours])
ink = gaussian(image,sigma=blurRad)
# xyi = np.array([np.tile(np.arange(image.shape[1])/dpmm,image.shape[0]),
#   np.repeat(np.arange(image.shape[0])/dpmm,image.shape[1]),
#   ink.flatten()])  #  a 3-row matrix with each vertex's (x,y,blurred ink).
# interior = np.full(xyi.shape[1],False)
# for (k,pix) in enumerate(xyi.T):


# A 2-column matrix with the (x,y) coords of each image point in the same order
# as ink.flatten() will give:
xyflat = np.array([np.tile(np.arange(image.shape[1])/dpmm,image.shape[0]),
  np.repeat(np.arange(image.shape[0])/dpmm,image.shape[1])]).T
interior = copa.contains_points(xyflat)
print(np.min(ink[interior]))
print(np.max(ink[interior]))
print(np.min(ink[~interior]))
print(np.max(ink[~interior]))
exit()

# HALFWIDTHS OF BOUNDARY CUTS:
contours = [np.concatenate((c,np.full((c.shape[0],1),maxRad)),axis=1)
  for c in contours]


# Now process it in some way to partition into "ink" & "no ink." A small amount
# of smoothing followed by a threshold would simplify the topology. On the
# other hand, it would "extremise" locally light grey or dark grey regions.
# What about this:
#   two scales: small & medium (still fairly small)
#   calculate local mean value at medium scale
#   do until converge:
#     smooth at small scale
#     threshold at 0.5
#     calculate local mean value at medium scale
#     subtract difference
#   Quite possibly there's already something like this or better in python's
#     image tools.
# Note: coords are array-style: [vert from top, horz from left]
(J,K) = image.shape

# perform skeletonization
print("Skeletonising...")
skeleton = skeletonize(image)

plt.ion()
(fig, ax) = plt.subplots(figsize=(9, 9))
ax.imshow(skeleton)
ax.set_xticks([]), ax.set_yticks([])
ax.axis([0, image.shape[1], image.shape[0], 0])
fig.canvas.flush_events()

# get strokewidth:
sw = ndimage.distance_transform_edt(image) * skeleton  #  what is '*' here?

scale = outputWidth/K  #  millimetres per pixel
(X,Y) = np.meshgrid(scale*np.arange(K), scale*np.arange(J)[::-1])
Z = -sw * (scale * bam)

#  shifts cuts downward so that the origin/safe-height can be a little higher:
Z = (Z-sanddepth) * (Z<0)
numToDo = np.sum(Z<0)  # number of pixels needing cutting still to be scheduled
XYZ = np.stack((X,Y,Z),2)
toDo = (Z<0)

def neighYet(j,k):
  # The number of neighbouring points still needing scheduling.
  return np.sum(toDo[max(0,j-1):min(J,j+2),max(0,k-1):min(K,k+2)]) - toDo[j,k]

def spiralAround(j,k):
  yield (j,k)
  for searchRad in range(1,max(J,K)):
    for jo in [j-searchRad,j+searchRad]:
      if jo>=0 and jo<J:
        for ko in range(max(0,k-searchRad+1),min(K,k+searchRad)):
          yield (jo,ko)
    for ko in [k-searchRad,k+searchRad]:
      if ko>=0 and ko<K:
        for jo in range(max(0,j-searchRad),min(J,j+searchRad+1)):
          yield (jo,ko)

curpos = (0,0)
paths = []
while numToDo>0:
  curpath = [] # np.zeros((2,0))
  # Try to find a loose end:
  for (jo,ko) in spiralAround(*curpos):
    if toDo[jo,ko] and neighYet(jo,ko)==1:
      curpath.append(np.array([jo,ko]))
      toDo[jo,ko]=False
      break
  if not len(curpath):  #  couldn't find a loose end; just grab any point
    for (jo,ko) in spiralAround(*curpos):
      if toDo[jo,ko]:
        curpath.append(np.array([jo,ko]))
        toDo[jo,ko]=False
        break
  # Now we can assume that we've found a starting point.
  while len(curpath)<numToDo and neighYet(jo,ko):
    for (jo,ko) in spiralAround(jo,ko):  # terminates at an immediate neighbour
      if toDo[jo,ko]:
        curpath.append(np.array([jo,ko]))
        toDo[jo,ko]=False
        break
  paths.append(curpath)
  numToDo -= len(curpath)
  curpos = (jo,ko)
  arcjk = np.array(curpath)
  ax.plot(arcjk[:,1], arcjk[:,0], '-r', lw=3)
  fig.canvas.flush_events()  #  plt.draw()
  print(numToDo)

print("All the paths have been found. Writing the gcode...")
nodes = np.row_stack(([0,0,0],[0,0,sfht]) +
  tuple(np.row_stack(([X[tuple(path[0])],Y[tuple(path[0])],sfht],
  XYZ[tuple(np.array(path).T)],
  [X[tuple(path[-1])],Y[tuple(path[-1])],sfht])) for path in paths) +
  ([0,0,sfht],)).T
taxis = np.row_stack((np.insert(np.cumsum([k for path in paths
  for k in (2,len(path))])-0,0,0),
  np.arange(2*len(paths)+1)%2))
cnc.ToolPath(nodes,taxis).PathToGCode(900,outfile)

# BG:301808 ; FG:CCBB88
