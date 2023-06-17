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
"""

import cnc
import numpy as np
from scipy import ndimage
from skimage import io
from skimage.morphology import skeletonize
import matplotlib.pyplot as plt

# PARAMETERS:
imfile = "../../../stamp.png"
outputWidth = 72  #  width of desired output, in millimetres
outfile = "stamp.gcode"  #  widest point is about 72 pixels

vBitAngle = 90  #  bit angle in degrees
vBitWidth = 3.175
sfht = 3  #  3 mm ain't much but ought to suffice
sanddepth = 0.2  #  cut this much deeper, so a groove remains after sanding 

image = io.imread(imfile, as_gray=True).astype(float)
# Now process it in some way to partition into "ink" & "no ink." A small amount
# of smoothing followed by a threshold would simplify the topology. On the
# other hand, it would "extremise" locally light grey or dark grey regions.
# What about this:
#   two scales: small & medium (still fairly small)
#   calculate local mean value at medium scale
#   do until converge:
#     smooth at small scale
$     threshold at 0.5
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

# Note: vBitAngle is the dihedral angle between the two sides of a groove cut
# by the bit, i.e., TWICE the angle between a groove wall (or bit edge) and the
# vertical.
# bit-angle multiplier; ratio of depth of cut to width of cut:
bam = 0.5/np.tan(vBitAngle*np.pi/360)

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
