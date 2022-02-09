# Produce necklace. Two gcode programs are needed:
#  cut oval contour for front of piece
#  cut path with v-bit
#  cut oval contour for back of piece (1st gcode program again)

import cnc as cnc
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
# Parameters:
# Radii (x,y,z); centres (x,y,z); margin multipliers of ovoids:
ovoids = [(40,31,4.8,0,0,0,1.03),(4,4,4,0,39,0,2.5)]
mep = -0.6  #  norm parameter used to merge ovoids' surfaces

vWaste = 1  #  extra depth to be cut away; orig depth should be oD+2*vW
imWidth = 62
imXlate = np.array([1,-2.5])
rad = 1.6  #  radius of ball nose
graveDepth = 2  #  but try 10deg bit this time, and run it repeatedly
safeHt = 3
cutVant = 10  #  height above (0,0,z)st pot(0,0,z)=0 whence engraving starts
delta = 0.01  #  step size for estimating gradients

###############################################################################
# imported paths from Inkscape:
impPaths = [
np.array([
[105.18656, 127.56475],
[108.17635, 125.02795],
[109.35415, 121.94755],
[111.89095, 121.76635],
[112.79695000000001, 120.76975],
[114.60895000000001, 118.23295],
[119.50135, 115.51495],
[123.21595, 113.70295],
[124.74595000000001, 111.53535],
[123.48775, 108.26695],
[120.67915, 106.36435999999999],
[118.05175000000001, 105.27716],
[112.16275000000002, 102.83095999999999],
[106.72676000000001, 97.84795799999999],
[104.37116000000002, 88.33496299999999],
[102.83096000000002, 82.62716499999999],
[103.10276000000002, 77.00996799999999],
[104.91476000000002, 74.83556899999999],
[104.81396000000001, 73.83896899999999],
[102.36436, 71.75517099999999],
[99.29755800000001, 69.21837],
[94.13335900000001, 69.76196999999999],
[89.87516300000001, 72.02696999999999],
[86.06996300000002, 76.19456699999999],
[84.07676400000001, 82.89896599999999],
[82.89896600000002, 88.78796399999999],
[81.72116500000001, 92.32135999999998],
[77.19116800000002, 91.95896099999999],
[67.22517100000002, 89.784562],
[59.25237600000002, 86.341764],
[55.99077800000002, 85.616964],
[53.997576000000024, 86.62376499999999],
[51.82317700000002, 87.61016299999999],
[51.188979000000025, 95.40175899999998],
[51.00777900000003, 100.92835999999998],
[52.36677700000003, 107.36094999999999],
[55.62837500000003, 109.89774999999999],
[58.61817500000003, 110.62254999999999],
[61.24557500000003, 110.53195],
[62.60457200000003, 112.70635],
[63.51057200000003, 113.88415],
[61.78917400000003, 117.41755],
[58.70877600000003, 124.48435],
[55.99077800000003, 129.01435],
[53.72577800000003, 134.08794],
[53.54457800000003, 139.07094],
[57.440377000000026, 141.78894],
[62.87637400000003, 143.51033999999999],
[65.50377300000002, 142.69493999999997],
[69.21837000000002, 139.97693999999998],
[73.11416900000002, 137.53073999999998],
[77.46296600000002, 131.91353999999998],
[81.17756600000003, 129.01434999999998],
[84.25796400000003, 128.38015],
[81.99296600000002, 133.09133999999997],
[84.71096500000003, 141.51713999999998],
[89.05976200000003, 145.50354],
[93.31796100000004, 148.04034],
[98.11976000000004, 148.22153999999998],
[100.83776000000005, 145.86593999999997],
[101.20016000000004, 141.15473999999998],
[102.01556000000004, 134.99393999999998],
[102.46856000000004, 130.19214999999997],
[103.64636000000004, 126.83994999999997],
[105.18656, 127.56475],
]), np.array([
[90.781163, 91.234161],
[90.41876300000001, 95.220559],
[91.86836200000002, 99.841157],
[93.95215900000002, 102.46856],
[97.66675800000003, 104.64296],
[100.83776000000003, 106.90796],
[103.28396000000004, 110.53195000000001],
[104.91476000000003, 111.52855000000001],
[108.62935000000003, 111.89095],
]), np.array([
[67.134572, 103.46516],
[73.56717, 103.46516],
[78.459567, 102.64976],
[84.167365, 101.65316],
[88.516163, 99.931758],
[91.868362, 99.84115800000001],
]), np.array([
[84.710965, 114.42775],
[86.975962, 110.98495],
[89.51276299999999, 108.26695],
[92.593161, 105.91135],
[95.764161, 103.64636],
]), np.array([
[92.32136, 132.91014],
[93.13676099999999, 127.38355000000001],
[95.85476, 121.85695000000001],
[98.116358, 112],
[101, 108],
]), np.array([
[126.38695, 83.261366],
[127.38355, 85.798164],
[125.84335, 91.868362],
[123.75955, 97.123159],
[118.68595, 105.73016],
]), np.array([
[102.83096, 82.627165],
[107.81396000000001, 80.633966],
[116.23975000000002, 79.184367],
[120.58855000000001, 79.00316699999999],
[123.39715000000001, 79.99976799999999],
[124.48435, 80.724567],
]), np.array([
[123.39715, 79.999768],
[127.74595, 77.009968],
[133.90673999999999, 74.654369],
[138, 72.020168],
[142, 72],
[144.41633999999996, 75.194565],
[145.23173999999997, 78.912564],
[144.05393999999998, 81.724565],
[140.79234, 83.08696499999999],
[135.26574, 82.358763],
[131.73234, 82.44596299999999],
[127.38355, 85.79816199999999],])]

###############################################################################
# PROCESS THE PATHS:
# Flip vertically:
procPaths = [path*[1,-1] for path in impPaths]

# Rotate all the paths:
theta = -0.4
R = np.array([[np.cos(theta),np.sin(theta)],[-np.sin(theta),np.cos(theta)]])
procPaths = [np.matmul(path,R) for path in procPaths]

# Calculate bounding box:
mins = np.min(np.array([np.min(p,0) for p in procPaths]),0)
maxs = np.max(np.array([np.max(p,0) for p in procPaths]),0)

# Centre the image:
procPaths = [path-(mins+maxs)/2 for path in procPaths]

# Scale it:
procPaths = [imXlate+path*imWidth/(maxs[0]-mins[0]) for path in procPaths]

###############################################################################
# Now how to control ball end to cut the oval?
# Unlike in the first attempt, we'll allow a more general shape, defined as a
# level set, and calculate the gradients numerically.

# "Potential": negative inside shape, positive outside.
# It is assumed that each (x,y) has at most one +ve z for which pot(x,y,z)=0.
def pot(x):
  potlist = np.array([(((x[0]-xc)/xr)**2+((x[1]-yc)/yr)**2)**2 +
    ((x[2]-zc)/zr)**4 for (xr,yr,zr,xc,yc,zc,mm) in ovoids])
  if np.min(potlist)<1e-6:
    return -1
  else:
    return np.sum(potlist**mep)**(1/mep)-1

# Iterated interpolation to find z s.t. pot(x,y,z)=0.
# We make the simplifying assumption that pot is increasing with z^2.
def xytoxyz(xy):
  # Here [LHI][PZ] stand for low,high,interpolated pot,z values. Typically we
  # expect that the interval [LZ,HZ] will bound the sought z value, but we
  # won't rely on that.
  LP = pot(np.array((*xy,0)))
  if LP > 0:
    return np.array((*xy,0))
  LZ = 0
  HZ = 10
  HP = pot(np.array((*xy,HZ)))
  IP = 1
  while np.abs(IP) > 1e-6:
    IZ = (LZ+HZ)/2
    IP = pot(np.array((*xy,IZ)))
    if IP>0:
      (HP,HZ) = (IP,IZ)
    else:
      (LP,LZ) = (IP,IZ)
  return np.array((*xy,IZ))

def unitNormal(xyz):
  # Not hugely efficient but symmetrical at least.
  grad = np.array([pot(xyz+[delta if k==j else 0 for k in range(3)]) -
    pot(xyz-[delta if k==j else 0 for k in range(3)]) for j in range(3)])
  if np.abs(np.linalg.norm(grad)) < 1e-12:
    breakpoint()
  return grad/np.linalg.norm(grad)

def bitPos(x,y,bitrad,boolPos=True):
  # Inputs x & y coordinates for a point on the oval's surface.
  # Outputs the corresponding (x,y,z) for the nose of the bit.
  # Set boolPos to False for carving the recess.
  xyz = xytoxyz((x,y))
  if bitrad==0:  #  engraving step!
    return xyz + np.array([0,0,-cutVant-graveDepth])
  else:
    grad = unitNormal(xyz)
    if not boolPos:
      # Carving the recess:
      xyz *= [-1,1,-1]
      grad *= [1,-1,1]
      if xyz[2] >= 0:
        # happens if no contact
        xyz[2] = 1+bitrad  #  bitrad will be subtracted again, leaving clearance of 1mm
    # return xyz + bitrad*grad - np.array([0, 0, bitrad + (
    #   vWaste+np.max([o[2] for o in ovoids]) if boolPos else 0)])
    retval = xyz + bitrad*grad
    retval[2] -= bitrad  #  if retval[2]<0: retval[2] -= bitrad
    if boolPos:
      retval[2] -= (vWaste + np.max([o[2] for o in ovoids]))
    return retval

def tanPath(xr,yr,xc,yc,mm):
  # Returns the path in 2D of the tangent point to cut.
  stepangle = 0.01  #  radians; about 0deg34'23"
  offset = 0.18
  N = int(mm*2*np.pi*xr/(stepangle*offset))  #  num nodes in path
  r = mm*(1-np.linspace(1,0,N)**1.8)
  th = stepangle*np.arange(N)
  return np.vstack((xc+xr*r*np.cos(th), yc+yr*r*np.sin(th))).T

def nosePaths():
  # Returns the path in 3D of the nose of the bit while cutting the upper
  # surface.
  altov = [ovoids[0][:4]+(ovoids[0][4]+0.8,)+ovoids[0][5:],
    ovoids[1][:4]+(ovoids[1][4]-3.2,)+ovoids[1][5:]]
  retval =[np.array([bitPos(x,y,rad) for (x,y) in tanPath(xr,yr,xc,yc,mm)]) for
    (xr,yr,zr,xc,yc,zc,mm) in altov]
  # This bit is a total hack; no idea why it's needed.
  retval[1][:,2] += np.linspace(1,1.4,retval[1].shape[0])**2
  return retval

def recessPaths():
  # Returns the path in 3D of the nose of the bit while cutting the recess.
  altov = [ovoids[0], ovoids[1][:-1]+(3,)]
  return [np.array([bitPos(x,y,rad,False) for(x,y) in tanPath(xr,yr,xc,yc,mm)])
    for (xr,yr,zr,xc,yc,zc,mm) in altov]

def pointPaths():
  # Returns the path in 3D of the nose of the v-bit while engraving the image.
  return [np.array([bitPos(x,y,0) for (x,y) in path]) for path in procPaths]

def writePathListToGCode(pathList,progName):
  nodes = np.zeros((3,1))
  taxis = np.zeros((2,1),dtype=int)
  for path in pathList:
    if np.linalg.norm(nodes[:2,-1]-path[0,:2])>0.3:
      if taxis[-1,-1]==1:
        taxis = np.concatenate((taxis, np.array([[nodes.shape[1]],[0]])),
          axis=1)
      nodes = np.concatenate((nodes, np.array([[nodes[0,-1], path[0,0]],
        [nodes[1,-1],path[0,1]], [safeHt,safeHt]])), axis=1)
    if taxis[-1,-1]==0:
      taxis = np.concatenate((taxis, np.array([[nodes.shape[1]],[1]])), axis=1)
    nodes = np.concatenate((nodes, path.T), axis=1)
  taxis = np.concatenate((taxis, np.array([[nodes.shape[1]],[0]])), axis=1)
  nodes = np.concatenate((nodes, np.array([[nodes[0,-1],0,0],
    [nodes[1,-1],0,0], [safeHt,safeHt,0]])), axis=1)
  tp = cnc.ToolPath(nodes, taxis)
  tp.PathToGCode(1200, progName)

boolGCode = True
if boolGCode:
  writePathListToGCode(nosePaths(), "Oval")
  #writePathListToGCode(recessPaths(), "Recess")
  #xyz = xytoxyz([0,0])
  #cutVant += xyz[2]
  #writePathListToGCode(pointPaths(), "Jaca")

boolPlot = False
if boolPlot:
  lines = [np.array([xytoxyz((x,y)) for x in np.arange(-41,41,0.4)])
    for y in np.arange(-32,45,0.4)]
  lines = [l[l[:,2]>0,:] for l in lines]  #  trim beyond piece
  lines = [np.concatenate((l,l*[[-1,1,-1]])) for l in lines]
  ax = plt.axes(projection='3d')
  for l in lines:
    ax.plot3D(l[:,0],l[:,1],l[:,2],'blue')

  # Trying to set aspect ratio to 1:1:1:
  # plt.gca().set_aspect('equal')  #  should work but doesn't.
  ax.set_box_aspect([ub-lb for (lb,ub) in (getattr(ax, f'get_{a}lim')()
    for a in 'xyz')])  #  https://github.com/matplotlib/matplotlib/issues/17172
  plt.show()

boolMeasure = True
if boolMeasure:
  print([np.sum(np.linalg.norm(np.diff(path,axis=0),axis=1))
    for path in procPaths])
# 170.9
#  17.4
#  14.4
#   8.9
#  15.3
#  14.1
#  12.9
#  28.4

# Add jitter to paths, & then plot in random colours, to check for overlap:
boolCheckPaths = False
if boolCheckPaths:
  fig = plt.figure()
  for p in procPaths:
    rp = p + np.random.normal(scale=0.2, size=p.shape)
    plt.plot(rp[:,0], rp[:,1], color=0.8*np.random.rand(3))
  plt.show()
