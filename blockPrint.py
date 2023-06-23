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
imfile = "W.jpg" ; downSampleRate = 1
# imfile = "steamboat.jpg" ; downSampleRate = 8
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
      # if c[1]<-55 and t[k] < 0.4*maxRad:
      #   breakpoint()  #  some ts are too shallow; work out why.
      #   # I bet find_contours has an off-by-one effect like len(diff(x)).
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
jkG = ([(jvec,k) for jvec in np.split(np.arange(len(t)),
  1+(np.diff(t==0)).nonzero()[0])[(1 if t[0]==0 else 0)::2]] for (k,t) in ktG)
rasterCuts = [[np.array([X[jvec,k],Y[jvec,k],-targetDepth[jvec,k]]).T
  for (jvec,k) in col] for col in jkG]

# Just a debugging utility:
l = lambda cl: np.concatenate([rc[[0,-1]] for rc in cl])

# NOW COMBINE AND ORDER ALL CUTS:
# Heuristic, but a good one:
#   Firstly plan the raster cuts as though we're only doing them:
#     each constant-x group is done in order, y either increasing or decresing
#     greedily go to whichever end is closer at end of each x
#   For each loop cut, find the shortest detour to select insertion point.
scheduledCuts = rasterCuts[0].copy()
for sublist in rasterCuts[1:]:
  # Firstly don't bother to schedule any raster cuts that are too short, since
  # the boundary cuts will handle such areas:
  sublist = [cut for cut in sublist if np.abs(cut[0,1]-cut[-1,1])>3*maxRad]
  # Now schedule the remaining ones in a mostly back&forth order:
  if np.abs(scheduledCuts[-1][-1,1]-sublist[-1][-1,1]) < np.abs(
    scheduledCuts[-1][-1,1]-sublist[0][0,1]):
    scheduledCuts += [rc[::-1] for rc in reversed(sublist)]
  else:
    scheduledCuts += sublist

# Find optimal insertion points for each loop cut (again, greedily):
for bc in boundaryCuts:
  # DL[j,k] is the extra motion added by inserting this loop after the
  # currently jth scheduled cut, starting & ending at its kth vertex.
  DL = np.array([np.linalg.norm(bc[:,:2]-scheduledCuts[j][-1,:2],axis=1) +
    np.linalg.norm(bc[:,:2]-scheduledCuts[j+1][0,:2],axis=1) 
    for j in range(len(scheduledCuts)-1)])
  opt = np.unravel_index(np.argmin(DL, axis=None), DL.shape)
  scheduledCuts.insert(opt[0]+1,np.concatenate((bc[opt[1]:],bc[:(opt[1]+1)])))

# Shift cuts downward a tiny bit so that you can sand after cutting. Also use
# raised origin to minimise damage if the early termination problem recurs.
extra = 30
sfht = 2 - extra  #  2mm above surface
scheduledCuts = [cut-[0,0,extra+sanddepth] for cut in scheduledCuts]

print("Writing the gcode...")
nodes = np.row_stack([
  np.row_stack([[*path[0,:2],sfht], path, [*path[-1,:2],sfht]])
  for path in scheduledCuts]).T
taxis = np.row_stack((np.insert(np.cumsum([k for path in scheduledCuts
  for k in (2,len(path))])-1,0,0),
  np.arange(2*len(scheduledCuts)+1)%2))
# cnc.ToolPath(nodes,taxis[:,:-1]).PathToGCode(900,outfile)
taxis = np.row_stack((np.cumsum([k for path in scheduledCuts
  for k in (2,len(path))])-2,
  np.arange(1,2*len(scheduledCuts)+1)%2))
cnc.ToolPath(nodes,taxis).PathToGCode(900,outfile)

# Now to check, plot the image & the boundary cuts in 3D:
if (viewCuts := True):
  ax = plt.figure().add_subplot(projection='3d')
  print(image.shape)
  faco = np.stack([image]*3+[0.6+0*image],axis=-1)
  print((X.shape,Y.shape,faco.shape))
  surf = ax.plot_surface(X, Y, 0*X, facecolors=faco,
    linewidth=0, antialiased=False)
  for boc in scheduledCuts:
    ax.plot(*boc.T, color='r')
  ax.set_xlim([0,M])
  ax.set_ylim([-M,0])
  ax.set_zlim([-M/2,M/2])
  plt.xlabel('x')
  plt.ylabel('y')
  plt.show()
  exit()

# BG:301808 ; FG:CCBB88
