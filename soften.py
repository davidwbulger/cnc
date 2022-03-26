from matplotlib import pyplot as plt
import numpy as np
from scipy import interpolate as terp, ndimage as nd
from skimage import io, measure

"""
Repeat until no change:
  replace black with union of circles entirely contained within black
  replace white with union of circles entirely contained within white

Ideally, cut boundary and then a grid to fill in interiors.
"""

ppm = 4.1  #  pixels per millimetre
bitrad = 3.175/2  #  in millimetres
rip = bitrad*ppm  #  radius in pixels
inra = int(rip)

def makeGCode(cuts, fname):
  pass

##  PHASE 1: ROUND THE BOUNDARIES SO ALL CUTS CAN BE REACHED  #################

lingrid = np.arange(-inra,inra+1)
(px,py) = np.meshgrid(lingrid,lingrid)
kernel = (px**2+py**2<rip**2).astype(np.uint8)

grsc=io.imread("wedjat.png",as_gray=True).astype(float)
image=(grsc>0.25).astype(np.uint8)
# Note: coords are array-style: [vert from top, horz from left]
(J,K) = grsc.shape

plt.ion()
(fig,ax) = plt.subplots()
ax.imshow(image,cmap=plt.cm.gray)
ax.set_xticks([]), ax.set_yticks([])
ax.axis([0, image.shape[1], image.shape[0], 0])
fig.canvas.flush_events()

flipper = 0
while True:
  newImage = image
  if flipper:  #  spread white & then black by flipping at start & middle
    newImage = 1-newImage
  # Spread one colour (and flip):
  newImage = (nd.convolve(newImage,kernel,mode='nearest')==0).astype(np.uint8) 
  # Spread the other colour (but don't flip):
  newImage = (nd.convolve(newImage,kernel,mode='nearest')>0).astype(np.uint8) 
  if not flipper:  #  spread black & then white by flipping at middle & end
    newImage = 1-newImage
  ax.imshow(newImage,cmap=plt.cm.gray)
  fig.canvas.flush_events()
  if np.array_equal(image, newImage):
    break
  image = newImage
  flipper = 1-flipper

# The real-time image is only for visualising infinite loop errors. If we've
# reached equilibrium, we can turn that off.
plt.ioff()

##  PHASE 2: CALCULATE CUTTING PATHS  #########################################

unew = np.arange(0,1.001,0.001)

# "Negative," i.e., the part that needs to be cut out of the background piece:
# Convolve background with kernel, then find boundary; may need to repeat to
# cut away interiors.
good = newImage - (grsc*(1-grsc)>0.125)
cuts = []
while not np.all(good):
  good = (nd.convolve(good,kernel,mode='nearest')>0).astype(np.uint8)
  cuts += measure.find_contours(good, 0.5)

ax.imshow(grsc,cmap=plt.cm.gray)
for (k,cu) in enumerate(cuts):
  cu = cu[:,::-1]
  if len(cu)>6:
    (tck,u) = terp.splprep(cu.T, s=0.1*len(cu))
    cu = terp.splev(unew,tck)
  else:
    cu = [cu[:,0], cu[:,1]]
  cuts[k] = cu
  ax.plot(cu[0], cu[1], color='r')
plt.show()
makeGCode(cuts, "Negative.gcode")

# Iris, the cut around the iris piece:
good = grsc*(1-grsc)>0.125
good = (nd.convolve(good, kernel, mode='nearest')>0).astype(np.uint8)
cuts = measure.find_contours(good, 0.5)
(fig,ax) = plt.subplots()
ax.imshow(grsc,cmap=plt.cm.gray)
ax.set_xticks([]), ax.set_yticks([])
ax.axis([0, image.shape[1], image.shape[0], 0])
for (k,cu) in enumerate(cuts):
  cu = cu[:,::-1]
  if len(cu)>6:
    (tck,u) = terp.splprep(cu.T, s=0.1*len(cu))
    cu = terp.splev(unew,tck)
  else:
    cu = [cu[:,0], cu[:,1]]
  cuts[k] = cu
  ax.plot(cu[0], cu[1], color='r')
plt.show()
makeGCode(cuts, "Iris.gcode")

# Brow & Wedjat, the cuts around the wedjat pieces:
good = grsc<0.25
good = (nd.convolve(good, kernel, mode='nearest')>0).astype(np.uint8)
cuts = measure.find_contours(good, 0.5)
(fig,ax) = plt.subplots()
ax.imshow(grsc,cmap=plt.cm.gray)
ax.set_xticks([]), ax.set_yticks([])
ax.axis([0, image.shape[1], image.shape[0], 0])
for (k,cu) in enumerate(cuts):
  cu = cu[:,::-1]
  if len(cu)>6:
    (tck,u) = terp.splprep(cu.T, s=0.1*len(cu))
    cu = terp.splev(unew,tck)
  else:
    cu = [cu[:,0], cu[:,1]]
  cuts[k] = cu
  ax.plot(cu[0], cu[1], color='r')
plt.show()

# Identify & separate the disjoint brow piece:
brix = np.argmax([np.mean(cu[1]) for cu in cuts])
makeGCode([cuts[brix]], "Brow.gcode")
cuts[brix:brix+1] = []
makeGCode(cuts, "Wedjat.gcode")
