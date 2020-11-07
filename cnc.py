# CNC utilities module.
# davidwbulger@gmail.com, 15 August 2020

# Assume Cartesian coordinates and millimetres.

import numpy as np
import struct
import re
import os
import numbers

######################################################################################################################################
####    DEFINE CLASSES:    ###########################################################################################################
######################################################################################################################################

class polyTri:
  def __init__(self,fname):
    # Read fname (an STL file) to build a polyTri structure
    fid = open(fname, "rb")
    if fid.read(5) == b"solid":
      # Read an ASCII format STL file:
      fid.close()
      fid = open(fname, "r")
      line = fid.readline()
      triList = []
      while line:
        if re.match("^\s*outer loop", line):
          triList.append(np.vstack([np.fromstring(re.sub("\s*vertex\s*","",fid.readline()),dtype=float,sep=' ') for ix in [0,1,2]]))
        line = fid.readline()
      fid.close()
      self.facets = np.array(triList,ndmin=3)
    else:
      # Read a binary format STL file:
      fid.read(75)  #  discard remainder of header
      numTri = struct.unpack('<I', fid.read(4))[0]  #  little-endian UINT32
      facets = np.zeros((numTri,3,3))
      for tx in range(numTri):
        fid.read(12)  #  discard normal vector
        facets[tx,:,:] = np.array([[struct.unpack('<f', fid.read(4))[0] for cx in range(3)] for rx in range(3)]) 
        fid.read(2)  #  discard 'attribute byte count'
      fid.close()
      self.facets = facets

  def __str__(self):
    return(self.facets.__str__())

######################################################################################################################################

class ToolPath:
  # Class to represent a sequence of tool motions.
  # Attributes are:
  #   nodes: a 3xN array of node coordinates.
  #   cutx: a 2xn array of cut indices (so kth cut follows .nodes[cutx[0,k]:cutx[1,k]]).
  #   rahe: a length n vector of heights at which to do rapid movements before each cut.

  def __init__(self, nodes, cutx, rahe):
    self.nodes = nodes
    self.cutx = cutx
    self.rahe = rahe

  def __str__(self):
    cutlen = sum([np.linalg.norm(self.nodes[:,j+1]-self.nodes[:,j]) for cx in self.cutx.T for j in range(cx[0],cx[1]-1)])
    m = np.min(self.nodes,1)
    M = np.max(self.nodes,1)
    return(f"ToolPath object with {self.cutx.shape[1]} cuts and {self.nodes.shape[1]} nodes in total."
      " Total cut length %.2f within box [%.2f,%.2f]x[%.2f,%.2f]x[%.2f,%.2f]." % (cutlen,m[0],M[0],m[1],M[1],m[2],M[2]))

  def afxform(self, A):
    # Applies the affine transformation described by the 4x4 matrix A to the nodes of self. Returns a copy.
    # Obviously this will need to be rethought if A[2,0:3] isn't [0,0,1].
    return ToolPath(np.matmul(A,np.vstack((self.nodes,[1]*self.nodes.shape[1])))[0:3,:], self.cutx, self.rahe)

  def PathToGCode(self,feedrate,fname):
    # Create a g-code file corresponding to a ToolPath object.
    # self.rahe[0] = Absolute Safe Height; assume no obstacles at that height anywhere
  
    # CONSTANTS:
    atol = 0.01  #  10 microns; well below machine precision
    coordnames = "XYZ"
  
    progname = (os.path.splitext(fname)[0]+"xxxxx")[:5].upper()
    # There are many popular file extensions for g-code! Assume '.gcode' if not provided:
    # (But some others are ".nc", ".tap", ".mpt", ".mpf")
    if not "." in fname: fname += ".gcode"
    fidout = open(fname, 'w')
  
    ## HEADER:
    fidout.write("%\nO" + progname + "\nG17 G21 G40 G49 G80 G90\n")
    fidout.write("G00 G54 X0. Y0.\nG43 H1 Z%.2f\n" % self.rahe[0])
  
    ## PATHS:
    curpos = np.array([0,0,self.rahe[0]])  #  current position, as a row
    setfeedrate = False  #  it's not set yet. The machine should be told the feedrate on the first cut instruction.
    motion = None

    for (cx,saht) in zip(self.cutx.T,self.rahe):
      # DO RAPID MOVE UP TO SAfeHeighT, TO NEXT (x,y), & then down to (z):
      newpos = self.nodes[:,cx[0]]
      if curpos[2]<saht:
        if not motion==0:
          fidout.write("G00 ")
          motion = 0
        fidout.write("Z%.2f\n" % saht)
      if not np.allclose(curpos[0:2], newpos[0:2], rtol=0, atol=atol):
        if not motion==0:
          fidout.write("G00 ")
          motion = 0
        fidout.write("X%.2f Y%.2f\n" % tuple(newpos[0:2]))
      for newpos in self.nodes[:,cx[0]:cx[1]].T:
        if not motion==1:
          fidout.write("G01 ")
          motion = 1
        fidout.write("X%.2f Y%.2f Z%.2f" % tuple(newpos))
        if not setfeedrate:
          fidout.write(" F%d" % feedrate)
          setfeedrate = True
        fidout.write("\n")
      curpos = newpos

    ## FOOTER:
    fidout.write("G00 Z%.2f\nG91 G28 X0. Y0.\nG90\nM30\n%%\n" % self.rahe[0])
    fidout.close()

  def plot(self, ax, color='black', linewidth=1):
    # Plot a ToolPath on the 3D axis "ax", in the given color.
    for cx in self.cutx.T:
      ax.plot(*self.nodes[:,cx[0]:cx[1]], color, linewidth=linewidth)

def compileToolPath(paths, hts):
  # Combine compatible lists of 3D paths and heights into a ToolPath object.
  nodes = np.hstack(paths)
  cumpathlengths = np.cumsum([path.shape[1] for path in paths])
  cutx = np.vstack((np.insert(cumpathlengths[0:-1],0,0), cumpathlengths))
  rahe = np.array(hts)
  return ToolPath(nodes, cutx, rahe)

def catToolPaths(TPList):
  # Joins the ToolPaths in TPList.
  indexoffsets = np.cumsum([0]+[tp.nodes.shape[1] for tp in TPList[0:-1]])
  return ToolPath(np.hstack([tp.nodes for tp in TPList]),
    np.hstack([tp.cutx + io for (tp,io) in zip(TPList, indexoffsets)]),
    np.concatenate([tp.rahe for tp in TPList]))
    
######################################################################################################################################

class PathGrid:
  # A PathGrid is a system of paths, each of constant y, with y increasing. Within each path, z varies as x increases. This definition
  # is not remotely isotropic, but the idea is to set up a PathGrid as a sequence of cuts, and then reorient it as necessary.

  def __init__(self, y, xz):
    # Create a PathGrid with ordinates 'y' and corresponding paths 'xz.'

    # Error checking:
    if not isinstance(y, np.ndarray):
      raise PathGridError("PathGrid attribute y must be a Numpy NDArray.")
    if not all(isinstance(xzp, np.ndarray) for xzp in xz):
      raise PathGridError("PathGrid attribute xz must be a list of Numpy NDArrays.")
    if any (np.diff(y) <= 0):
      raise PathGridError("PathGrid attribute y must be strictly increasing.")
    if any(any(np.diff(xzp[0,:])<0) for xzp in xz):
      raise PathGridError("Each 2D path in PathGrid attribute xz must be non-decreasing.")
    if y.shape[0] != len(xz):
      raise PathGridError("PathGrid attributes y and xz must have the same length.")

    # Okay, it's valid.
    self.y = y
    self.xz = xz

  def __str__(self):
    return(f"PathGrid object with {len(self.y)} paths and {self.nodeCount()} nodes in total.")

  def maxoxy(self):
    # returns the scalar maximum of z over x & y
    return max([max(xzp[1,:]) for xzp in self.xz])

  def pospar(self):
    return PathGrid(self.y, [pospar(xzp) for xzp in self.xz])

  def whereBelow(self, other):
    # Trims a PathGrid to only the parts where self is below other.
    # Return structure is not a PathGrid; it's a list (one entry per y value) of lists (one entry per disconnected interval) of xz
    # cutting paths with fixed y.
    dg = (other-self).pospar()  #  "difference grid"
    segl = [np.hstack((
      np.argwhere(np.diff(np.insert(np.logical_or(xzp[1,:-1]>1e-3, xzp[1,1:]>1e-3),0,False).astype(int))==1),
      2 + np.argwhere(np.diff(np.append(np.logical_or(xzp[1,:-1],xzp[1,1:]), False).astype(int))==-1)
      )) for xzp in dg.xz]
    return [[np.concatenate((
          np.array([[xzd[0,se[0]]],[np.interp(xzd[0,se[0]],xzs[0,:],xzs[1,:])]]),
          xzs[:,np.logical_and(xzs[0,:]>xzd[0,se[0]], xzs[0,:]<xzd[0,se[1]-1])],
          np.array([[xzd[0,se[1]-1]],[np.interp(xzd[0,se[1]-1],xzs[0,:],xzs[1,:])]])),1)
        for se in seg]
      for (yp,seg,xzd,xzs) in zip(self.y,segl,dg.xz,self.xz)]

  def maxabs(self):
    return max([max(abs(xzp[1,:])) for xzp in self.xz])

  def nodeCount(self):
    return sum([xzp.shape[1] for xzp in self.xz])

  # Now overload arithmetic operators. Generally these apply to the z component. The y attribute won't be altered, and a copy isn't
  # necessary; there's no foreseeable use-case of a PathGrid being constructed and used in calculations, and then having its y values
  # changed.

  def __neg__(self):
    # return PathGrid(self.y, [np.array([[1],[-1]]) * xzp for xzp in self.xz])
    return self * (-1)

  def __mul__(self, other):
    if isinstance(other, numbers.Number):
      return PathGrid(self.y, [np.array([[1],[other]]) * xzp for xzp in self.xz])
    raise PathGridError("Multiplier of PathGrid object must be a scalar.")

  def __rmul__(self, other):
    return self * other

  def __add__(self, other):
    if isinstance(other, numbers.Number):
      return PathGrid(self.y, [np.array([[0],[other]]) + xzp for xzp in self.xz])
    else:
      if not isinstance(other, PathGrid):
        raise PathGridError("Addend to a PathGrid must be another PathGrid or a scalar.")
      if not np.array_equal(self.y, other.y):
        raise PathGridError("Two PathGrid addends must have the same y attributes.")
    return PathGrid(self.y, [addPLFs(sxzp,oxzp) for (sxzp,oxzp) in zip(self.xz, other.xz)])

  def __sub__(self, other):
    return self + (-other)

  def __rsub__(self, other):
    return -self + other

  def castToMold(self, ballrad, ltol, floor):
    # Converts "self," describing an actual cut surface, to a corresponding "mold," describing the lowest surface the drillbit's
    # spherical head's centre can traverse. Both are PathGrids.
    # The drillbit's head's radius is ballrad. The lateral tolerance for path equality is ltol. Don't cut below floor.
    # If no floor is needed, set floor to np.NINF.
  
    mxz = [None for _ in self.xz]
    for j,(yp,xzp) in enumerate(zip(self.y, self.xz)):
      # j & y describe the cross-section that the drillbit centre traverses.
      # Initialise a very fine discretisation of this cross-section:
      xgrid = np.linspace(xzp[0,0], xzp[0,-1], 2+int((xzp[0,-1]-xzp[0,0])/ltol))
      zgrid = floor * np.ones(xgrid.shape)
      # Now "push up" the values in zgrid to sit on spheres & cylinders around neighbouring paths:
      for alty,altxz in zip(self.y, self.xz):
        if np.abs(alty-yp) < ballrad:
          # The neighbouring path's avoidance zone intersects our xz plane.
          xzrad = np.sqrt(ballrad**2-(alty-yp)**2)  #  reduced radius of avoidance zone
          # First avoid spheres:
          for node in altxz.T:
            vex = np.abs(xgrid - node[0]) < xzrad  #  indeX of points coincident with sphere under VErtical projection
            zgrid[vex] = np.maximum(zgrid[vex], node[1] + np.sqrt(xzrad**2 - (xgrid[vex]-node[0])**2))
          # Now avoid cylinders:
          for nx in range(1,altxz.shape[1]):
            segwid = altxz[0,nx]-altxz[0,nx-1]
            if segwid > ltol:
              seglen = np.linalg.norm(altxz[:,nx]-altxz[:,nx-1])
              segslo = (altxz[1,nx]-altxz[1,nx-1])/segwid
              edgex = altxz[0,nx-1]-(altxz[1,nx]-altxz[1,nx-1])*xzrad/seglen
              edgez = altxz[1,nx-1]+(altxz[0,nx]-altxz[0,nx-1])*xzrad/seglen
              vex = np.logical_and(xgrid>edgex, xgrid<edgex+segwid)
              zgrid[vex] = np.maximum(zgrid[vex], edgez + segslo * (xgrid[vex]-edgex))
      # Now reduce load by approximating zgrid more efficiently:
      mxz[j] = fitPWL(np.vstack([xgrid,zgrid]),ltol)
    return PathGrid(self.y, mxz)
  
  def roundJoint(self,ballrad,ltol):
    # The idea is to minimally modify the PathGrid "self" so as to be tangible from above & below with a ball-end drill bit with end
    # radius "ballrad."
  
    # This works best under the condition that all parts of the original surface are tangible from above OR below.
    # Under that condition, we can round the concave & convex extremities in either order with identical results.
    return (-((-(self.castToMold(ballrad,ltol,np.NINF))).castToMold(2*ballrad,ltol,np.NINF))).castToMold(ballrad,ltol,np.NINF)
  
  def plot(self, ax, color='black'):
    # Plot a PathGrid on the 3D axis "ax", in the given color.
    for (yp,xzp) in zip(self.y, self.xz):
      ax.plot(xzp[0,:], yp*np.ones(xzp.shape[1]), xzp[1,:], color, linewidth=1)

  def pacePathGrid(self, shpaid, abssh, cude):
    # "Pace" a target PathGrid (self), by creating a ToolPath that works its way down to the target depth from an initial stockheight
    # (described by shpaid---a PathGrid, or a scalar to indicate starting with a flat block), cutting no more than the prescribed
    # cutting depth (cude) on each path. As close as possible to a constant cut depth is left for the final cut, with the idea that it
    # will result in an evener finish. The "absolute stock height" (abssh) is a height at which rapid movement is always safe, used
    # between passes.
  
    wath = shpaid - self  #  a PathGrid representing the waste depth
    maxwaste = wath.maxoxy()
    numcuts = max(0, int(np.ceil(maxwaste/cude)))
    if numcuts < 1:
      raise PathGridError("No cutting required in pacePathGrid.")
    cutlevels = [maxwaste - cude*l for l in range(1,numcuts)] + [0]  #  includes target cut at end

    # Here we want to build the ToolPath in a loop. Unfortunately (though understandably), iteratively appending to numpy arrays
    # would create new copies at each iteration, and be laughably inefficient. But it's hard to predict how much memory to allocate
    # for the arrays. Therefore, in this loop we'll build two lists: a list of single-cut paths, and a same-length list of safe
    # heights for rapid motion between them. Then use compileToolPath to convert that to a ToolPath.
    scpaths = []
    safehts = [abssh]
    steptar = shpaid
    curpos = None  #  so it knows not to add a rapid move before first cut
    for cule in cutlevels:
      startOfLevel = True  #  so it knows to avoid at curheight rather than steptar
      curheight = steptar  #  set current stock height to target height at previous step
      steptar = shpaid - (wath - cule).pospar()
      cutlist = steptar.whereBelow(curheight)  #  just the actual cutting needing doing

      # Reverse cut list in the y direction if necessary:
      if curpos != None and abs(curpos[1]-steptar.y[-1]) < abs(curpos[1]-steptar.y[0]):
        cutlist.reverse()
        stepy = reversed(self.y)
      else:
        stepy  = self.y

      for (yp, cutp) in zip(stepy, cutlist):
        if len(cutp):
          if curpos != None:
            if abs(curpos[0]-cutp[-1][0,-1]) < abs(curpos[0]-cutp[0][0,0]):
              cutp = [np.fliplr(cu) for cu in reversed(cutp)]  #  do these cuts with decreasing x
            if startOfLevel:
              safehts.append(maxcrop(curheight, curpos, [cutp[0][0,0],yp]))
              startOfLevel = False
            else:
              safehts.append(max(maxcrop(curheight,[curpos[0],yp], [cutp[0][0,0],yp]), maxcrop(steptar,curpos, [cutp[0][0,0],yp])))
          safehts += [maxcrop(curheight,[cutp[j-1][0,-1],yp],[cutp[j][0,0],yp]) for j in range(1,len(cutp))]
          scpaths += [np.vstack((cu[0,:], yp*np.ones(cu.shape[1]), cu[1,:])) for cu in cutp]
          curpos = [scpaths[-1][0,-1], yp]
    return compileToolPath(scpaths, safehts)

def addPLFs(a,b):
  # Inputs two lists of 2D numpy arrays of two rows each (x & z), increasing in x, describing each PLF (piece-wise linear function).
  # Outputs another PLF, in the same format, being a+b.
  # Assumes a & b start & end at the same x values.
  x = np.sort(np.unique(np.concatenate((a[0,:],b[0,:]))))
  z = np.interp(x, a[0,:], a[1,:]) + np.interp(x, b[0,:], b[1,:])
  return np.vstack([x,z])
  
def pospar(plf):
  # Inputs a piecewise linear function (described by a 2-row Numpy array, with the first row non-descending).
  # Outputs max(0,plf) in the same format & over the same range.

  px = plf[1,:]>0  #  indices of positive values
  csx = np.nonzero(np.diff(px))[0]  #  indices of segments crossing horizontal axis
  nodes = np.column_stack((
    [plf[:,px]] +
    list(((plf[0,cx]*plf[1,cx+1]-plf[0,cx+1]*plf[1,cx])/(plf[1,cx+1]-plf[1,cx]),0) for cx in csx) +
    [(np.array([[1],[0]])*plf[:,[0,-1]])[:, np.logical_not(px[[0,-1]])]]))
  return nodes[:,np.argsort(nodes[0,:])]
  
# This was a method of PathGrid, but I want it to work sensibly if the 'paid' argument is just a constant, too.
def maxcrop(paid, p0, p1):
  # Returns the greatest height in paid within the closed region [x[0],x[1]]*[y[0],y[1]]*[-inf,inf].

  if isinstance(paid, PathGrid):
    # The rectangle is delimited axially by the two points, which can be input as any two opposite corners.
    x = [min(p0[0],p1[0]),max(p0[0],p1[0])]
    y = [min(p0[1],p1[1]),max(p0[1],p1[1])]
  
    retval = np.NINF
    for (yp,xzp) in zip(paid.y,paid.xz):
      if yp>=y[0] and yp<=y[1]:
        for xe in x:  #  each endpoint, i.e., x[0] & then x[1]
          if xzp[0,0]<xe and xzp[0,-1]>xe:
            retval = max(retval, np.interp(xe, xzp[0,:], xzp[1,:]))
        retval = np.amax(xzp[1,np.logical_and(xzp[0,:]>=x[0], xzp[0,:]<=x[1])], initial=retval)
  else:
    retval = paid  #  which is in this case hopefully just a numeric value
  return retval

######################################################################################################################################

class Error(Exception):  #  Base class for exceptions in this module.
  pass

class PathGridError(Error):  #  Exception raised for errors in the input.
  # Attribute:
  #   message -- explanation of the error
  def __init__(self, message):
    self.message = message

class ToolPathError(Error):  #  Exception raised for errors in the input.
  # Attribute:
  #   message -- explanation of the error
  def __init__(self, message):
    self.message = message

######################################################################################################################################
####    DEFINE UTILITIES:    #########################################################################################################
######################################################################################################################################

# Affine transformation matrices ("ATMs") should probably be a class ... but that seems like a lot of reinventing, and also I may not
# actually need much functionality.

def rotAround(x,y,theta):
  # Returns an ATM fixing z and rotating x and y an angle theta around (x,y).
  c = np.cos(theta)
  s = np.sin(theta)
  return np.array([[c, -s, 0, (1-c)*x+s*y], [s, c, 0, -s*x+(1-c)*y], [0, 0, 1, 0], [0, 0, 0, 1]])

def xlate(x,y):
  # Returns an ATM fixing z and translating x and y by (x,y).
  return np.array([[1, 0, 0, x], [0, 1, 0, y], [0, 0, 1, 0], [0, 0, 0, 1]])

def hackaspect(ax):
  # add a bounding box to force the aspect ratio to be 1:1:1
  xl = ax.get_xlim()
  yl = ax.get_ylim()
  zl = ax.get_zlim()
  M = 0.5 * max(np.diff(xl),np.diff(yl),np.diff(zl))[0]  #  max radius
  xl = np.mean(xl) + np.array([-M,M])
  yl = np.mean(yl) + np.array([-M,M])
  zl = np.mean(zl) + np.array([-M,M])
  ax.plot(xl, [yl[0],yl[0]], [zl[0],zl[0]], 'white')
  ax.plot(xl, [yl[1],yl[1]], [zl[1],zl[1]], 'white')

######################################################################################################################################
# THE FUNCTIONS meanroundingup & fitPWL ARE ABOUT EFFICIENTLY TRACING A PATH WITH LINE SEGMENTS:

def meanroundingup(a,b):
  return(int(0.75+(a+b)/2))  #  mean of integers a & b, but rounds up if sum is odd.

def fitPWL(xy, ltol):
  # Inputs a 2xN array of plane points, assumed to be increasing in x, and a lateral tolerance, and outputs a piecewise-linear
  # approximation to the heights. Specifically, returns a subset of the entered nodes, minimising node count subject to preserving
  # the path to within a tolerance.

  # Use a greedy method to determine breakpoints. (Should then be refined, but maybe later.)
  N = xy.shape[1]
  keepNodes = [False]*N
  keepNodes[0] = keepNodes[-1] = True  #  definitely want the first & last points!
  fitto = 0  #  index of how much of the curve we've approximated so far
  while fitto < N-1:
    gotoatleast = fitto+1  #  we'll always be able to reach at least the next node
    gotoatmost = N-1
    while gotoatleast < gotoatmost:
      gototry = meanroundingup(gotoatleast,gotoatmost)
      # now need a test of collinearity that doesn't break if segments are vertical...
      start = xy[:,fitto,None]  #  'None' = no drop
      midpoints = xy[:,fitto+1:gototry]
      finish = xy[:,gototry,None]
      dotprods = np.matmul((finish-start).T,midpoints-start)
      straightEnough = np.all(dotprods > 0)
      straightEnough = straightEnough and np.all(dotprods < np.matmul((finish-start).T,finish-start))
      if straightEnough:
        crossprods = np.matmul(np.matmul((finish-start).T,np.array([[0,-1],[1,0]])),midpoints-start)
        straightEnough = np.all(np.abs(crossprods) < ltol*np.linalg.norm(finish-start)/2)
      if straightEnough:
        gotoatleast = gototry
      else:
        gotoatmost = gototry-1
    keepNodes[gotoatmost] = True
    fitto = gotoatmost
  return(xy[:,keepNodes])

