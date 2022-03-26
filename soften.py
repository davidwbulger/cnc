import cnc
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

def makeGCode(cuts, depth, fname):
  # Note the change from matrix (& image) orientation to Cartesian:
  path3s = [np.concatenate(([cu[0]/ppm], [-cu[1]/ppm],
    np.full((1,len(cu[0])),-paz))).T
    for paz in depth*np.arange(1,4)/3 for cu in cuts]
  # Offset so tight box's lower right corner is origin:
  boxmin = np.min(np.concatenate(path3s), axis=0)
  # but hack Iris.gcode so its origin is pupil centreL
  if fname=="Iris.gcode":
    lix = np.argmin(path3s[0][:,0])
    rix = np.argmax(path3s[0][:,0])
    boxmin=np.array([(path3s[0][lix,0]+path3s[0][rix,0])/2,path3s[0][lix,1],0])
    irrad=(path3s[0][rix,0]-path3s[0][lix,0])/2
    print(f"Iris radius is {irrad}mm.")  #  57mm.
  path3s = [p-boxmin*[1,1,0] for p in path3s]
  cnc.cutSequence(path3s, 1000, 6, fname)

##  PHASE 1: ROUND THE BOUNDARIES SO ALL CUTS CAN BE REACHED  #################

lingrid = np.arange(-inra,inra+1)
(px,py) = np.meshgrid(lingrid,lingrid)
kernel = (px**2+py**2<rip**2).astype(np.uint8)

grsc=io.imread("RecordCase/wedjat.png",as_gray=True).astype(float)
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
makeGCode(cuts, 5, "Negative.gcode")

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
makeGCode(cuts, 6, "Iris.gcode")  #  extra depth for positives to ensure free

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
brix = np.argmin([np.mean(cu[1]) for cu in cuts])  #  recall, matrix indexing
makeGCode([cuts[brix]], 6, "Brow.gcode")
cuts[brix:brix+1] = []
makeGCode(cuts, 6, "Wedjat.gcode")
