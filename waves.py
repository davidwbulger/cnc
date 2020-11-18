# Create g-code for "Raindrops" carve.

import cnc
import numpy as np
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

########################################################################################################################
# PARAMETERS:

np.random.seed(1)
prain = 0.01  #  proportion of frames with a raindrop
damping = 0.002
wavespeed = 0.25  #  actually "stability parameter s" = c^2(delta t)^2/(delta x)^2 from Grigoryan
numFrames = 200
dropsize = 60
dropradius = 0.36

R = 120  #  semimajor radius
r = 26  #  semiminor radius
gridres = 0.75  #  offset should be a multiple of this
maxdepth = 10
mindepth = 0.5
# snaptime

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

# The callback used to advance time one step & draw:
def stepRipples(framenum, animlines):
  # Now draw it:
  zg = zdeck[:,:,framenum]
  ax.set_title(f"Frame {framenum} ; {np.mean(zg)}")
  for (ix,(xl,yl,zlix)) in enumerate(zip(downx,downy,downzix)):
    animlines[ix].set_data(xl,yl)
    animlines[ix].set_3d_properties(zg.ravel()[zlix])
  return animlines

# Add a raindrop:
def oneDrop(zgrid):
  # We don't want every raindrop identical, so let's randomly perturb the heights at the points in a 2x2 grid.
  # Randomly choose a grid in the ellipse:
  hitInterior = False
  while not hitInterior:
    xix = np.random.randint(0,2*N)
    yix = np.random.randint(0,2*n)
    hitInterior = interior[yix:yix+2, xix:xix+2].all()
  #for i in [yix,yix+1]:
  #  for j in [xix,xix+1]:
  #    zgrid[i,j] += dropsize # * (np.random.rand() - 0.5)
  #zgrid[yix,xix] += dropsize
  return zgrid + dropsize * (np.exp(-((xgrid-xgrid[yix,xix])**2+(ygrid-ygrid[yix,xix])**2)/(2*dropradius**2))
    - 0.25 * np.exp(-((xgrid-xgrid[yix,xix])**2+(ygrid-ygrid[yix,xix])**2)/(8*dropradius**2)))

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
ghost = np.logical_and(np.logical_not(interior), np.logical_or(np.roll(interior,-1,0),
  np.logical_or(np.roll(interior,1,0), np.logical_or(np.roll(interior,-1,1), np.roll(interior,1,1)))))
ghostlist = np.flatnonzero(ghost)  #  np.argwhere(ghost.ravel())
spectres = np.array([reflectGhost(x,y,R,r) for (x,y) in zip(xgrid.ravel()[ghostlist], ygrid.ravel()[ghostlist])])
##  print(spectres)
##  print(reflectGhost(xgrid[0,0], ygrid[0,0], R, r))
##  print(xgrid.ravel()[ghostlist])
##  print(ghostlist)
zgrid = 0 * xgrid
lag1z = zgrid.copy()
lag2z = zgrid.copy()
zgrid = oneDrop(zgrid)  #  hit the pool with one random raindrop

# Calculate the deck of zgrids:
zdeck = np.zeros(zgrid.shape + (numFrames,))
for fnum in range(numFrames):
  # Calculate reflected heights for ghost points:
  if False:
    excrapolator = interp2d(xlist, ylist, zgrid.ravel()[inlist])
    for (eachghost, im) in zip(ghostlist, spectres):
      zgrid.ravel()[eachghost] = excrapolator(im[0],im[1])
  # Now do a step of the PDE:
  if True:
    lag2z = lag1z.copy()
    lag1z = zgrid.copy()
    zgrid = ((2-4*wavespeed)*lag1z - (1-damping)*lag2z +
      wavespeed*(np.roll(lag1z,1,0) + np.roll(lag1z,-1,0) + np.roll(lag1z,1,1) + np.roll(lag1z,-1,1))) / (1+damping)
  # Now maybe a raindrop:
  if np.random.rand() < prain: zgrid = oneDrop(zgrid)
  zdeck[:,:,fnum] = zgrid.copy()

# Also create downsampled lines for plotting:
downyvals = ymargin[range(1,2*n,3)]
downx = [xlist[ylist==yval] for yval in downyvals]
downy = [ylist[ylist==yval] for yval in downyvals]
downzix = [inlist[ylist==yval] for yval in downyvals]

########################################################################################################################
# SET UP THE ANIMATION:

# Now actually create the figure & pass the above callback to FuncAnimation:
fig = plt.figure()
ax = fig.add_subplot(projection="3d")
animlines = [ax.plot(xl,yl,zgrid.ravel()[zlix],'black',linewidth=0.5)[0] for (xl,yl,zlix) in zip(downx,downy,downzix)]
ax.set_xlim3d([-R, R])
ax.set_ylim3d([-R, R])
ax.set_zlim3d([-R, R])
framenum = 0
perfunctoryGlobalVariable = animation.FuncAnimation(fig, stepRipples, numFrames, fargs=(animlines,), interval=42)  #  unlikely to achieve this speed!
plt.show()
########################################################################################################################

##  # CREATE SINGLE GRID:
##  y = np.linspace(offset-semiminor, semiminor-offset, num=nrows)
##  x = [(semimajor/semiminor)*np.arange(-np.sqrt(semiminor*semiminor-yp*yp), np.sqrt(semiminor*semiminor-yp*yp), min(1,offset)) for yp in y]
##  # print([len(xp) for xp in x])
##  dist = [np.sqrt((xp-semimajor)**2+yp**2)*nripples*np.pi/semimajor for (xp,yp) in zip(x,y)]  #  in radians of ripples
##  # z = [(np.cos(distp)-1)*depth/2 for distp in dist]
##  # z = [(np.cos(distp)-1)*depth/2 + depth-dpc for distp in dist]
##  z = [np.sqrt(1.001 - (yp/semiminor)**2 - (xp/semimajor)**2) + 0.2*np.cos(0.32*xp-0.2*yp) for (xp,yp) in zip(x,y)]
##  minz = min([min(zp) for zp in z])
##  maxz = max([max(zp) for zp in z])
##  slope = (maxdepth-mindepth)/(maxz-minz)
##  intercept = -(slope*minz+maxdepth)
##  z = [intercept + slope*zp for zp in z]
##  xz = [np.vstack((xp,zp)) for (xp,zp) in zip(x,z)]
##  taPG = cnc.PathGrid(y,xz)
##  
##  # CONVERT TO TOOLPATH:
##  # Note that, currently, PacePathGrid is the only way to convert a PathGrid to a ToolPath. This code 'manually' handles the layers,
##  # so we'll start by pretending that the safe cut per depth is huge to get the target ToolPath:
##  taTP = taPG.pacePathGrid(0, 2, 20)
##  
##  fig = plt.figure()
##  ax = fig.add_subplot(111, projection='3d') 
##  taTP.plot(ax, 'red')
##  cnc.hackaspect(ax)
##  plt.show()
##  
##  # LAYERS:
##  # This would automate the whole layers thing:
##  # layers = cnc.catToolPaths([taTP.afxform(np.array([[1,0,0,0],[0,1,0,0],[0,0,1,zadj],[0,0,0,1]])) for zadj in np.flip(np.arange(0,depth,dpc))])
##  # layers.PathToGCode(frate, "ripples.gcode")
##  
##  # But instead I'll just do the one cut & manually jog & repeat.
##  # taTP.PathToGCode(frate, "waves.gcode")
