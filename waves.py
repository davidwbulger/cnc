# Create g-code for "Raindrops" carve.

import cnc
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

########################################################################################################################
# PARAMETERS:

initSeed = 10
action_animate = False
action_gcode = True
gcodeFrame = 111
gcodeToolSeq = [{'bitrad':1.5,'cude':3,'ds':3},{'bitrad':0.5,'cude':0.7,'ds':1}]
transpose = True  #  Rotate ellipse after construction so that cuts are the short way
prain = 0.05  #  proportion of frames with a raindrop
gridres = 0.5 # 0.75  #  offset should be a multiple of this
damping = 0.03*gridres
wavespeed = 0.23  #  actually "stability parameter s" = c^2(delta t)^2/(delta x)^2 from Grigoryan
numFrames = 200
dropsize = 50*gridres**1.5
dropradius = 0.6*gridres**0.25

R = 120  #  semimajor radius
r = 26  #  semiminor radius
maxdepth = 10
mindepth = 0.5

# THIS ISN'T QUITE SO SIMPLE. SHOULD ROUGH CUT WITH LARGER BIT BEFORE FINAL CUT(S?) WITH SMALLER. (OTHERWISE IT'S 20
# CUTS!)
dpc = 0.5  #  depth per cut
frate = 150
offset = 4  #  0.7

########################################################################################################################
# FUNCTIONS NEEDED:

# The following three functions are needed for ellipse boundary calculations, and adapted from David Eberly's
# https://www.geometrictools.com/Documentation/DistancePointEllipseEllipsoid.pdf.

def robustLength(x,y):
  if (x := abs(x)) < (y := abs(y)):
    (x,y) = (y,x)
  return x*np.sqrt(1+(y/x)**2)

def getRoot(sr, z0, z1, g):
  n0 = sr*z0
  s0 = z1-1
  s1 = 0 if g<0 else robustLength(n0,z1)-1
  s = 0
  for i in range(60):
    s = (s0+s1)/2
    ratio0 = n0/(s+sr)
    ratio1 = z1/(s+1)
    g = ratio0**2 + ratio1**2 - 1
    if g>0:
      s0 = s
    else:
      s1 = s
  return s

def reflectGhost(x,y,R,r):
  # For simplicity, we assume that the ellipse is axial and centred.
  xneg = (x<0) ; x = abs(x)
  yneg = (y<0) ; y = abs(y)
  if y==0:
    im = [2*R-x, 0]
  elif x==0:
    im = [0, 2*r-y]
  else:
    z0 = x/R
    z1 = y/r
    g = z0**2+z1**2-1
    if g==0:
      im = [x,y]
    else:
      sr = (R/r)**2
      sbar = getRoot(sr, z0, z1, g)
      im = [x*(sr-sbar)/(sr+sbar), y*(1-sbar)/(1+sbar)]
  if xneg: im[0] = -im[0]
  if yneg: im[1] = -im[1]
  return im

# Add a raindrop:
def oneDrop(zgrid):
  # We don't want every raindrop identical, so let's randomly perturb the heights at the points in a 2x2 grid.
  # Randomly choose a grid in the ellipse:
  hitInterior = False
  while not hitInterior:
    jx = np.random.randint(0,2*N)
    iy = np.random.randint(0,2*n)
    hitInterior = interior[iy:iy+2, jx:jx+2].all()
  return zgrid + dropsize * (np.exp(-((xgrid-xgrid[iy,jx])**2+(ygrid-ygrid[iy,jx])**2)/(2*dropradius**2))
    - 0.25 * np.exp(-((xgrid-xgrid[iy,jx])**2+(ygrid-ygrid[iy,jx])**2)/(8*dropradius**2)))

# Animation functions:

# The callback used to advance time one step & draw:
def stepRipples(framenum, animlines):
  # Now draw it:
  zg = zdeck[:,:,framenum]
  ax.set_title(f"Seed {initSeed}; Frame {framenum}")
  for (ix,(xl,yl,zlix)) in enumerate(zip(downx,downy,downzix)):
    animlines[ix].set_data(xl,yl)
    animlines[ix].set_3d_properties(zg.ravel()[zlix])
  return animlines

# Two functions allowing crude user control of the animation (space to start/stop; cursor L/R to set time direction):
def update_time():
  t = 0
  t_max = numFrames
  while True:
    t = (t + anim.direction) % numFrames
    yield t

def on_press(event):
  if event.key.isspace():
    if anim.running:
      anim.event_source.stop()
    else:
      anim.event_source.start()
    anim.running ^= True
  elif event.key == 'left':
    anim.direction = -1
  elif event.key == 'right':
    anim.direction = +1

########################################################################################################################
# ADJUSTMENTS:

# Ensure the top, bottom, left & right of the ellipse are neat:
R = gridres * (np.floor(R/gridres) - 0.01)
r = gridres * (np.floor(r/gridres) - 0.01)

########################################################################################################################
# SET UP THE ELLIPSE'S ARRAYS:

# Just once, at initialisation, we need to find the image of each ghost point when reflected into the ellipse through
# the nearest point on the ellipse's boundary.
N = int(np.ceil(R/gridres))  #  2N+1 columns (including ghost margins)
n = int(np.ceil(r/gridres))  #  2n+1 rows (including ghost margins)
xmargin = gridres * np.arange(-N,N+1)
ymargin = gridres * np.arange(-n,n+1)
(xgrid, ygrid) = np.meshgrid(xmargin, ymargin)
interior = (xgrid/R)**2 + (ygrid/r)**2 < 1
inlist = np.flatnonzero(interior)  #  np.argwhere(interior.ravel())
xlist = xgrid.ravel()[inlist]
ylist = ygrid.ravel()[inlist]

# Calculate "ghost" nodes, the positions of their reflected images ("spectres"), and the linear combinations of
# interior nodes that we'll use to circumpolate them:
ghost = np.logical_and(np.logical_not(interior), np.stack((
  np.roll(interior,(1,0),(0,1)), np.roll(interior,(1,1),(0,1)), np.roll(interior,(0,1),(0,1)), 
  np.roll(interior,(-1,1),(0,1)), np.roll(interior,(-1,0),(0,1)), np.roll(interior,(-1,-1),(0,1)), 
  np.roll(interior,(0,-1),(0,1)), np.roll(interior,(1,-1),(0,1))
  ),axis=2).any(axis=2))
ghostlist = np.flatnonzero(ghost)  #  np.argwhere(ghost.ravel())
spectres = np.array([reflectGhost(x,y,R,r) for (x,y) in zip(xgrid.ravel()[ghostlist], ygrid.ravel()[ghostlist])])

# excrapolator_ind[i,j] will be the linear index to the jth vertex of the triangle approximating the ith spectre, and
# excrapolator_wt[i,j] will be its corresponding barycentric coordinate.
excrapolator_ind = np.zeros((len(ghostlist), 3),dtype=np.intc)
excrapolator_wt = np.zeros((len(ghostlist), 3))
corners = np.array([[0,0], [0,1], [1,0], [1,1]])
for (sx,spec) in enumerate(spectres):
  distance = np.Inf
  for jx in range(max(0,int(np.floor(spec[0]/gridres))+N-1), min(2*N,int(np.ceil(spec[0]/gridres))+N+1)):
    for iy in range(max(0,int(np.floor(spec[1]/gridres))+n-1), min(2*n,int(np.ceil(spec[1]/gridres))+n+1)):
      for omit in range(4):
        xyind = np.ravel_multi_index(tuple(([iy,jx]+corners[list(range(omit))+list(range(omit+1,4)),:]).T), xgrid.shape)
        xy = np.vstack((xgrid.ravel()[xyind], ygrid.ravel()[xyind])).T
        centroid = np.mean(xy, axis=0)
        candis = np.linalg.norm(centroid-spec)
        if candis < distance:
          distance = candis
          excrapolator_ind[sx,:] = xyind
  xyind = excrapolator_ind[sx,:]
  xy = np.vstack((xgrid.ravel()[xyind], ygrid.ravel()[xyind])).T
  excrapolator_wt[sx,:] = np.linalg.solve(np.vstack((xy.T, [1,1,1])), np.vstack((spec[0],spec[1],1))).T

zgrid = 0 * xgrid
lag1z = zgrid.copy()
lag2z = zgrid.copy()
zgrid = oneDrop(zgrid)  #  hit the pool with one random raindrop

# Calculate the deck of zgrids:
zdeck = np.zeros(zgrid.shape + (numFrames,))
np.random.seed(initSeed)
for fnum in range(numFrames):
  # Calculate reflected heights for ghost points:
  zgrid *= interior
  for (eachghost, ind, wt) in zip(ghostlist, excrapolator_ind, excrapolator_wt):
    zgrid.ravel()[eachghost] = sum(zgrid.ravel()[ind] * wt)
  # Now do a step of the PDE:
  if True:
    lag2z = lag1z.copy()
    lag1z = zgrid.copy()
    zgrid = ((2-6*wavespeed)*lag1z - (1-damping)*lag2z +
      wavespeed*(np.roll(lag1z,1,0) + np.roll(lag1z,-1,0) + np.roll(lag1z,1,1) + np.roll(lag1z,-1,1) +
      0.5*(np.roll(lag1z,(1,1),(0,1)) + np.roll(lag1z,(-1,1),(0,1)) + np.roll(lag1z,(1,-1),(0,1)) +
      np.roll(lag1z,(-1,-1),(0,1))))) / (1+damping)
  # Now maybe a raindrop:
  if np.random.rand() < prain: zgrid = oneDrop(zgrid)
  zdeck[:,:,fnum] = zgrid.copy()

########################################################################################################################
# SET UP THE ANIMATION:

if action_animate:
  # Create downsampled lines for plotting:
  downyvals = ymargin[range(1,2*n,1)]
  downx = [xlist[ylist==yval] for yval in downyvals]
  downy = [ylist[ylist==yval] for yval in downyvals]
  downzix = [inlist[ylist==yval] for yval in downyvals]
  
  # Now actually create the figure & pass the above callback to FuncAnimation:
  fig = plt.figure()
  fig.canvas.mpl_connect('key_press_event', on_press)
  ax = fig.add_subplot(projection="3d")
  animlines = [ax.plot(xl,yl,zgrid.ravel()[zlix],'black',linewidth=0.5)[0] for (xl,yl,zlix) in zip(downx,downy,downzix)]
  ax.set_xlim3d([-0.6*R, 0.6*R])
  ax.set_ylim3d([-0.6*R, 0.6*R])
  ax.set_zlim3d([-0.6*R, 0.6*R])
  fig.tight_layout()
  ax.set_clip_on(False)
  anim = animation.FuncAnimation(fig, stepRipples, frames=update_time, fargs=(animlines,), interval=42, repeat=True)
  anim.running = True
  anim.direction = -1
  mng = plt.get_current_fig_manager()
  mng.window.state("zoomed")
  plt.show()

########################################################################################################################
# WRITE THE GCODE:

if action_gcode:
  # Construct PathGrid:
  union = np.logical_or(interior, ghost)
  zproc = zdeck[:,:,gcodeFrame]
  zmax = np.max(zproc[interior])
  zmin = np.min(zproc[interior])
  zproc = ghost * 2 + interior * (-maxdepth + (maxdepth-mindepth) / (zmax-zmin) * (zproc-zmin))  #  2>0 prevents cutting
  if transpose:
    (xgrid, ygrid) = (ygrid.T, xgrid.T)
    zproc = zproc.T
    union = union.T
  whichRows = union.any(axis=1)
  y = ygrid[whichRows,0]  #  ymargin[whichRows]
  x = [xgrid[ro,union[ro,:]] for ro in np.flatnonzero(whichRows)]
  z = [zproc[ro,union[ro,:]] for ro in np.flatnonzero(whichRows)]
  xz = [np.vstack((xp,zp)) for (xp,zp) in zip(x,z)]
  pg = cnc.PathGrid(y,xz)
  regions = [(lambda x,y:x<1),(lambda x,y:x>-1)]
  tplist = pg.MultiToolPG(0, 5, gcodeToolSeq, regions)
  tplist[0].PathToGCode(500, "coarseEllipse.gcode")
  tplist[1].PathToGCode(300, "fineEllipse.gcode")

