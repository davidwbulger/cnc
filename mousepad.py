# Program to cut Nikau's mousepad.

import cnc
import numpy as np
from scipy import ndimage
from skimage import io
from skimage.morphology import skeletonize
import matplotlib.pyplot as plt
from skimage.util import invert

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

bam = 2.5  #  bit-angle multiplier; ratio of depth of cut to width of cut
outputWidth = 200  #  width of desired output, in millimetres
scale = outputWidth/K  #  millimetres per pixel
(X,Y) = np.meshgrid(scale*np.arange(K), scale*np.arange(J)[::-1])
Z = -sw * (scale * bam)
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

"""
path = np.zeros((3,0))
for ii in range(numToDo):
  # Find a nearby place to start a cutpath:
  searchRadius = 1  #  in the supremum norm
  while not Z[curpos]:
    searchRadius *= 2
    searchBlock = np.array([[max(0,curpos[0]-searchRadius), max(0,curpos[1]-searchRadius)],
      [min(J-1,curpos[0]+searchRadius), min(K-1,curpos[1]+searchRadius)]])
    zMinRelLoc = (
      np.unravel_index(Z[searchBlock[0,0]:(1+searchBlock[1,0]), searchBlock[0,1]:(1+searchBlock[1,1])].argmin(),
      diff(searchBlock,axis=0)[0]+1))
    if curpos[0]-searchRadius >= 0 and np.any(Z[curpos[0]-searchRadius, 
  
if False:
  # display results
  fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 4),
                           sharex=True, sharey=True)
  
  ax = axes.ravel()
  
  ax[0].imshow(image, cmap=plt.cm.gray)
  ax[0].axis('off')
  ax[0].set_title('original', fontsize=20)
  
  ax[1].imshow(skeleton, cmap=plt.cm.gray)
  ax[1].axis('off')
  ax[1].set_title('skeleton', fontsize=20)
  
  fig.tight_layout()
  plt.show()

"""
