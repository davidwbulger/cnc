# Produce necklace. Two gcode programs are needed:
#  cut oval contour for front of piece
#  cut path with v-bit
#  cut oval contour for back of piece (1st gcode program again)

import cnc as cnc
import numpy as np
import matplotlib.pyplot as plt

# Parameters:
ovWidth = 99
ovHeight = 77
ovDepth = 9
vWaste = 0.5  #  extra depth to be cut away; orig depth should be oD+2*vW
imWidth = 77
imOffset = np.array([1.5,-3.5])
rad = 3  #  radius of ball nose
graveDepth = 1.5  #  excessive?
safeHt = 3

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
procPaths = [imOffset+path*imWidth/(maxs[0]-mins[0]) for path in procPaths]

# Now how to control ball end to cut the oval? The equation for the oval should
# be
# ((2x/oW)^2+(2y/oH)^2)^2+(2z/oD)^4=1.
# Given x & y, want to choose z such that the distance from (x,y,z) to the
# surface is r.
# For parallel plane-intersecting paths (of constant y, say), sample a grid of
# surface points' heights AND GRADIENTS, and calculate position of ball centre
# (and thus of nose) for ball to be tangent.

def bitPos(x,y,bitrad):
  # Inputs x & y coordinates for a point on the oval's surface.
  # Outputs the corresponding (x,y,z) for the nose of the bit.
  e = 1-((2*x/ovWidth)**2+(2*y/ovHeight)**2)**2
  if e > 0:
    z = e**0.25*ovDepth/2
  else:
    z = 0;
  if bitrad>0:  #  cutting the oval surface
    grad = np.array([x/ovWidth**2*((x/ovWidth)**2+(y/ovHeight)**2),
      y/ovHeight**2*((x/ovWidth)**2+(y/ovHeight)**2), z**3/ovDepth**4])
    grad = grad / np.linalg.norm(grad)
    return np.array([x,y,z])+rad*grad+np.array(
      [0,0,-bitrad-0.5*ovDepth-vWaste])
  else:  #  engraving step!
    return np.array([x,y,z])+np.array([0,0,-0.5*ovDepth-graveDepth])

def tanPath():
  # Returns the path in 2D of the tangent point to cut. Only a function for the
  # sake of tidiness.
  stepangle = 0.01  #  radians; about 0deg34'23"
  offset = 0.8
  N1 = int(np.pi*ovWidth/(stepangle*offset))  #  num nodes in path's 1st part
  N2 = int(2*np.pi/stepangle)  #  num nodes in path's 2nd part (outer circumf)
  r = np.concatenate((1-np.linspace(1,0,N1)**1.5, np.ones(N2)))
  th = stepangle*np.arange(N1+N2)
  return np.vstack((0.5*ovWidth*r*np.cos(th), 0.5*ovHeight*r*np.sin(th))).T

def nosePath():
  # Returns the path in 3D of the nose of the bit while cutting the oval
  # surface.
  return np.array([bitPos(x,y,rad) for (x,y) in tanPath()])

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
  # cnc.ToolPath(nodes, taxis).PathToGCode(1200, progName)
  tp = cnc.ToolPath(nodes, taxis)
  tp.PathToGCode(1200, progName)

writePathListToGCode([nosePath()], "Oval")
writePathListToGCode(pointPaths(), "Jaca")

boolPlot = False
if boolPlot:
  for path in procPaths[0:]:
    plt.plot(path[:,0], path[:,1], color="blue")
  phi = np.linspace(0,2*np.pi,361)
  plt.plot(0.5*ovWidth*np.cos(phi), 0.5*ovHeight*np.sin(phi), color="black")
  xy = tanPath()
  plt.plot(xy[:,0], xy[:,1], color="grey")
  plt.show()
