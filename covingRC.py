# One of three CNC tasks for the record cabinet. Here we're doing the fancy
# coving around the edge of the top piece.

import cnc as cnc
import numpy as np
import matplotlib.pyplot as plt

phi = 0.5+np.sqrt(1.25)
###############################################################################
# Parameters:
# Board is currently 19x457x614
# Changing to 17x457x380
# Side panel width is 380
bitRad = 2
originHeight = 60
w = 30         #  "tooth" base width
a = 0.5*w      #  "tooth" length
h = 17+bitRad  #  actual width of top plus bit radius
b = h*4/3      #  extra tooth extension at bottom vs top of board
c = (2+np.sqrt(8))/5  #  coefs in quadratic determining tooth shape
d = (3+np.sqrt(8))/5  #  coefs in quadratic determining tooth shape
tas = int(np.ceil(380/w-np.sqrt(0.125)))  #  number of teeth along side
taf = int(np.floor((457-2*(a+b))/w-np.sqrt(2))) 
print((tas,taf,a,b,(taf+np.sqrt(2))*w+2*(a+b)))
safeHt = 0

###############################################################################
class Tooth:

  def __init__(self,p0,p1):
    # Tooth is assumed to point to the right as we move from p0 to p1 (each 2D)
    self.p0 = p0
    self.p1 = p1
    self.R = np.array([(p1-p0)@np.array([[0,-1],[1,0]]),p1-p0])/w

  def pathAtTheta(self, theta, neighbourHeight):
    # First create the path in the tooth's own coordinates (that is, the coords
    # with origin at (p0,0), and rotated around the z-direction so that
    # p1-p0=(0,w,0)):
    t = np.linspace(-b/(2*a),1+b/(2*a),99)
    x = 2*a*np.minimum(t,1-t) + b*np.sin(theta)
    y = w * (t - c*x/(a+b) + d*(x/(a+b))**2)
    xy = np.concatenate(([x],[y])).T

    # Now transform it to world coordinates:
    xy = self.p0 + xy@self.R

    # Now remove any vertices that would cut into neighbouring teeth:
    z = h*np.cos(theta)
    # At corners we're getting some strange behaviour. I thought this would be
    # simpler ... but essentially we want the long run of True in the middle of
    # the following list:
    inclu = np.array([neighbourHeight(v)<z for v in xy])
    prix = int(len(inclu)/2)
    ultx = prix
    while prix>0 and inclu[prix-1]:
      prix -= 1
    inclu[:prix] = False
    while ultx<len(inclu)-1 and inclu[ultx+1]:
      ultx += 1
    inclu[ultx+1:] = False
    xy = np.array([[*v,z-h-originHeight] for (v,wh) in zip(xy,inclu) if wh])
    return xy

  def heightAtXY(self,xy):
    # First convert xy to the tooth's own coordinates:
    xy = (xy-self.p0) @ self.R.T

    # Now solve for t & theta:
    t = xy[1]/w + c*xy[0]/(a+b) - d*(xy[0]/(a+b))**2
    sintheta = (xy[0] - 2*a*np.min([t,1-t]))/b
    costheta = (1 if sintheta<0 else (
      -np.inf if sintheta>=1 else np.sqrt(1-sintheta**2)))
    return h*costheta

###############################################################################
# zip neighbours utility:
def zipNeigh(L):
  return [(a,b+c) for (a,b,c) in zip(
    L, [[]]+[[d]for d in L[:-1]], [[d] for d in L[1:]]+[[]])]
# Converts, e.g., ['Albania','Brazil', 'China', 'Denmark'] into
# [('Albania', ['Brazil']),
#  ('Brazil', ['Albania', 'China']),
#  ('China', ['Brazil', 'Denmark']),
#  ('Denmark', ['China'])]

###############################################################################
# Create the list of teeth:
r2 = np.sqrt(2)
steps = np.array([[0,0]] + [[0,-w]]*tas + [[w/r2,-w/r2]] + [[w,0]]*taf +
  [[w/r2,w/r2]] + [[0,w]]*tas)
joins = np.cumsum(steps, axis=0)
teeth = [Tooth(p0,p1) for (p0,p1) in zip(joins[:-1], joins[1:])]

# Now create the paths & the gcode:
toothPaths = [[tooth.pathAtTheta(theta,
  lambda x:max(n.heightAtXY(x) for n in neighbours))
  for (tooth,neighbours) in zipNeigh(teeth)]
  for theta in np.linspace(np.pi/2,0,15,endpoint=False)]
thetaPaths = [np.concatenate(sl) for sl in toothPaths]
cnc.cutSequence(thetaPaths, 1000, safeHt, "Truffules.gcode")

# Also a primer to initiate the cut:
thetaPaths = [np.concatenate(toothPaths[0])]
cnc.cutSequence(thetaPaths, 1000, safeHt, "TrufPrimer.gcode")

print(np.min(thetaPaths[0],axis=0))
print(np.max(thetaPaths[0],axis=0))
