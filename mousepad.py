# Program to cut Nikau's mousepad.

import cnc
import numpy as np
from scipy import ndimage
from skimage import io
from skimage.morphology import skeletonize
import matplotlib.pyplot as plt

image = io.imread("mousepad.jpg", as_gray=True).astype(float) < 0.5
# Note: coords are array-style: [vert from top, horz from left]
(J,K) = image.shape

# perform skeletonization
print("Skeletonising...")
skeleton = skeletonize(image)

# (labeled, numComps) = ndimage.label(skeleton, structure=np.ones((3,3)))

plt.ion()
(fig, ax) = plt.subplots(figsize=(9, 9))
ax.imshow(skeleton)
ax.set_xticks([]), ax.set_yticks([])
ax.axis([0, image.shape[1], image.shape[0], 0])
# ax.plot([10,1000],[500,1500],'-b', lw=3)
fig.canvas.flush_events()  #  plt.show()

# get strokewidth:
sw = ndimage.distance_transform_edt(image) * skeleton

# bam = 2.5  #  placeholder guess
vBitAngle = 20  #  in degrees ; probably need a less acute bit
bam = 0.5/np.tan(vBitAngle*np.pi/360)  #  bit-angle multiplier; ratio of depth of cut to width of cut

outputWidth = 200  #  width of desired output, in millimetres
orht = 5  #  start with bit 5mm above surface
scale = outputWidth/K  #  millimetres per pixel
(X,Y) = np.meshgrid(scale*np.arange(K), scale*np.arange(J)[::-1])
Z = -sw * (scale * bam)
Z = (Z-orht) * (Z!=0)  #  shifts cuts downward so that the origin/safe-height can be a little higher
numToDo = np.sum(Z<0)  #  number of pixels needing cutting that remain to be scheduled
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
    for (jo,ko) in spiralAround(jo,ko):  #  will terminate at an immediate neighbour
      if toDo[jo,ko]:
        curpath.append(np.array([jo,ko]))
        toDo[jo,ko]=False
        break
  paths.append(curpath)
  numToDo -= len(curpath)
  curpos = (jo,ko)
  arcjk = np.array(curpath)
  # arcxy = XYZ[tuple(arcjk.T)]
  ax.plot(arcjk[:,1], arcjk[:,0], '-r', lw=3)
  fig.canvas.flush_events()  #  plt.draw()
  print(numToDo)

print("All the paths have been found. Writing the gcode...")
# prenodes = ([0,0,0],
#   (np.row_stack(([X[tuple(path[0])],Y[tuple(path[0])],0], XYZ[tuple(np.array(path).T)],
#   [X[tuple(path[-1])],Y[tuple(path[-1])],0])) for path in paths),
#   [0,0,0])
# breakpoint()
nodes = np.row_stack(([0,0,0],) +
  tuple(np.row_stack(([X[tuple(path[0])],Y[tuple(path[0])],0], XYZ[tuple(np.array(path).T)],
  [X[tuple(path[-1])],Y[tuple(path[-1])],0])) for path in paths) +
  ([0,0,0],)).T
taxis = np.row_stack((np.insert(np.cumsum([k for path in paths for k in (2,len(path))])-1,0,0),
  np.arange(2*len(paths)+1)%2))
cnc.ToolPath(nodes,taxis).PathToGCode(1200,"Mousepad.gcode")