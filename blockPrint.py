# Project with Nikau to cut blocks for printing. Based on BibleCovers.py.
# Started 2023 June 15.

# Find unused imports:
import cnc
import numpy as np
from scipy import ndimage
from scipy.signal import decimate
from skimage.filters import gaussian
from skimage import io
from skimage import measure
# from skimage.morphology import skeletonize
from matplotlib.path import Path
import matplotlib.pyplot as plt

# PARAMETERS:
# imfile = "cicada.jpg" ; downSampleRate = 8
# imfile = "W.jpg" ; downSampleRate = 1
imfile = "steamboat.jpg" ; downSampleRate = 8
outputWidth = 72  #  width of desired output, in millimetres
outfile = "printBlock.gcode"  #  widest point is about 72 pixels
blurRad = 1.2  #  for antialiasing the boundary cut paths, if dsRate==1

vBitAngle = 90  #  bit angle in degrees
maxDepth = 3.2  #  3.2 mm
sanddepth = 0.2  #  cut this much deeper, so a groove remains after sanding 

# CALCULATED FROM PARAMETERS:
# Note: vBitAngle is the dihedral angle between the two sides of a groove cut
# by the bit, i.e., TWICE the angle between a groove wall (or bit edge) and the
# vertical.
# bit-angle multiplier; ratio of halfwidth of cut to depth of cut:
bam = np.tan(vBitAngle*np.pi/360)
maxRad = maxDepth * bam
image = io.imread(imfile, as_gray=True).astype(float)[:,::-1]  #  flip LR
if imfile=="steamboat.jpg":  #  hack to standardise range
  image -= np.min(image)
  image /= np.max(image)
if downSampleRate > 1:
  image = decimate(decimate(image, downSampleRate, axis=0),
    downSampleRate, axis=1)
  # Somehow, the above expands the range of colour values.
  image = (image-0.5)*0.5/np.max(np.abs(image-0.5)) + 0.5
else:
  image = gaussian(image,sigma=blurRad)
dpmm = image.shape[1]/outputWidth  #  dots per millimetre
M = max(image.shape)/dpmm

# FUNCTIONS:
def bocut(W,cmat,nmat):
  # Determine boundary cut given contour cmat & normals nmat.
  t = maxRad*np.ones((cmat.shape[0],1))
  for (k,(c,n)) in enumerate(zip(cmat,nmat)):
    Wx = W-c
    obx = np.sum(Wx*Wx,axis=1) < 2*maxRad*Wx@n
    if np.any(obx):
      Wo = Wx[obx]
      t[k] = np.min(np.sum(Wo*Wo,axis=1) / (2*Wo@n))
      if c[1]<-55 and t[k] < 0.4*maxRad:
        breakpoint()  #  some ts are too shallow; work out why.
        # I bet find_contours has an off-by-one effect like len(diff(x)).
  # Path of bit vertex:
  try:
    retval = np.concatenate((cmat[1:]+t[1:]*nmat,-bam*t[1:]),axis=1)
  except:
    breakpoint()
  return retval

# FIRST CALCULATE THE BOUNDARY CUTS, FOR SMOOTHER EDGES:
# contours = [c[:,::-1]/dpmm for c in measure.find_contours(image, 0.5)]
# contours = [c for c in measure.find_contours(image, 0.5)]
contours=[c@[[0,-1/dpmm],[1/dpmm,0]] for c in measure.find_contours(image,0.5)]
image = (image>=0.5)
# path = Path(vertices=contours[0])
# iv = np.array([[path.contains_point((j,k)), image[j,k]]
#   for j in range(image.shape[0]) for k in range(image.shape[1])])
# print([f(iv[iv[:,0]==bv][:,1]) for f in (np.min,np.max) for bv in (0,1)])
# exit()
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
  k = 0  #  np.argmax([l>10 and l<20 for l in (n.shape[0] for n in normals)])
  print(len(contours))
  print(np.array([[f([c.shape[z] for c in contours])
    for f in [np.min, np.mean, np.max]] for z in range(2)]))
  # Display the image and plot all contours found
  (fig, ax) = plt.subplots()
  ax.imshow(image, cmap=plt.cm.gray)
  for contour in contours:
    ax.plot(dpmm*contour[:, 0], -dpmm*contour[:, 1], linewidth=1, color='r')
  for (v,n) in zip(dpmm*contours[k][1:], 0.03*M*normals[k]):
    ax.plot([v[0],v[0]+n[0]], [-v[1],-v[1]-n[1]], linewidth=1, color='b')
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
X = np.array([np.arange(0,image.shape[1])/dpmm]*image.shape[0])
Y = np.array([-np.arange(0,image.shape[0])/dpmm]*image.shape[1]).T
W = np.array([X[~image],Y[~image]]).T  #  all inked pixels as rows
boundaryCuts = [bocut(W,c,n) for (c,n) in zip(contours,normals)]

# COMPUTE RASTER CUTS:
dpo = int(dpmm*maxRad)  #  dots per offset
targetDepth =np.minimum(maxRad,ndimage.distance_transform_edt(image)/dpmm)/bam
# Note that
#   np.split(s,1+(np.diff(s==0)).nonzero()[0])[(1 if s[0]==0 else 0)::2]
# gives a list of contiguous arrays of indices of nonzero values. Each will
# yield a raster cut.
# Generator for raster cuts' k indices & targetDepth columns:
ktG = ((k,targetDepth[:,k]) for k in range(0,image.shape[1],dpo))
jkG = ((jvec,k) for (k,t) in ktG for jvec in np.split(np.arange(len(t)),
  1+(np.diff(t==0)).nonzero()[0])[(1 if t[0]==0 else 0)::2])
rasterCuts = [np.array([X[jvec,k],Y[jvec,k],-targetDepth[jvec,k]]).T
  for (jvec,k) in jkG]

# Now to check, plot the image & the boundary cuts in 3D:
if (viewBoundaryCuts := False):
  ax = plt.figure().add_subplot(projection='3d')
  # ax.imshow(image, cmap=plt.cm.gray)
  print(image.shape)
  faco = np.stack([image]*3+[0.6+0*image],axis=-1)
  print((X.shape,Y.shape,faco.shape))
  surf = ax.plot_surface(X, Y, 0*X, facecolors=faco,
    # cmap=plt.cm.gray,
    linewidth=0, antialiased=False)
  for boc in boundaryCuts:
    ax.plot(*np.concatenate((boc,boc[0:1])).T, color='r')
  for rac in rasterCuts:
    ax.plot(*rac.T, color='b')
  ax.set_xlim([0,M])
  ax.set_ylim([-M,0])
  ax.set_zlim([-M/2,M/2])
  plt.xlabel('x')
  plt.ylabel('y')
  plt.show()
  exit()

exit()

# NOW COMBINE AND ORDER ALL CUTS:
allCuts = [(True,bc) for bc in boundaryCuts]+[(False,rc) for rc in rasterCuts]
# Create table of all points at which a cut can be started. For loop cuts,
# this is anywhere along the cut; otherwise, just the two ends.
xyjk = np.concatenate([
  np.concatenate([ac[:,:2],j*np.ones((len(ac),1)),np.arange(len(ac))[:,None]],
  axis=1) if loop else
  np.array([[ac[0,0],ac[0,1],j,0],[ac[-1,0],ac[-1,1],j,len(ac)-1]])
  for (loop,ac) in allCuts])
fovea = np.array([0,0])  #  initial "current point"
orderedCuts = []
# yts = np.arange(len(allCuts),dtype=int)  #  yet to schedule
yts = np.full(len(xyjk),True)
# while np.any(yts):  #  do above first so I can look at this -i style
#   rx = np.argmin([np.sum((row[:2]-fovea)**2 HIPPO np.inf



###############################################################################
exit()
###############################################################################



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
Z = -sw * (scale / bam)

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
