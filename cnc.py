# CNC utilities module.
# davidwbulger@gmail.com, started 5 August 2020

# Assume Cartesian coordinates and millimetres.

import numpy as np
import struct
import re
import os
import numbers
#import warnings
#warnings.filterwarnings("error")

############################################################################################################
####    DEFINE CLASSES:    #################################################################################
############################################################################################################

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
          triList.append(np.vstack([
            np.fromstring(re.sub("\s*vertex\s*","",fid.readline()),dtype=float,sep=' ') for ix in [0,1,2]]))
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
        facets[tx,:,:] = np.array([[struct.unpack('<f', fid.read(4))[0]
          for cx in range(3)] for rx in range(3)]) 
        fid.read(2)  #  discard 'attribute byte count'
      fid.close()
      self.facets = facets

  def __str__(self):
    return(self.facets.__str__())

############################################################################################################

class ToolPath:
  # Class to represent a sequence of tool motions.
  # Attributes are:
  #   nodes: a 3xN array of node coordinates.
  #   taxis: a 2xR array of taxis settings, so, e.g., a column [n,0] switches to mode G0 when moving from node n.
  #          Second row should just be 0s & 1s. Require taxis[0,0]=0 so we have a taxis mode at outset.

  def __init__(self, nodes, taxis):
    if nodes[:,[0,-1]].any():
      raise ToolPathError("First and last node must be the origin.")
    if taxis[0,0]!=0 or taxis[0,-1]>=nodes.shape[1]-1 or np.any(np.diff(taxis[0,:])<1):
      breakpoint()
      raise ToolPathError(
        "Taxis change nodes (first row of taxis) must increase from first node index, stopping before the last.")
    if not all(np.logical_or(taxis[1,:]==0, taxis[1,:]==1)):
      raise ToolPathError("Taxis changes (second row of taxis) must be binary.")
    self.nodes = nodes
    self.taxis = taxis

  def __str__(self):
    cuts = [nodes[:,taxis[0,j]:(taxis[0,j+1]+1)] for j in range(taxis.shape[1]-1) if taxis[1,j]==1]
    cutlen = sum([sum(np.linalg.norm(np.diff(cut,axis=1),axis=0)) for cut in cuts])
    raps = [nodes[:,taxis[0,j]:(taxis[0,j+1]+1)] for j in range(taxis.shape[1]-1) if taxis[1,j]==0]
    raplen = sum([sum(np.linalg.norm(np.diff(rap,axis=1),axis=0)) for rap in raps])
    m = np.min(self.nodes,1)
    M = np.max(self.nodes,1)
    return(f"ToolPath object with {len(cuts)} cuts and {self.nodes.shape[1]} nodes in total."
      " Total cut length %.2f and rapid movement length %.2f within box [%.2f,%.2f]x[%.2f,%.2f]x[%.2f,%.2f]." %
      (cutlen,raplen,m[0],M[0],m[1],M[1],m[2],M[2]))

  def afxform(self, A):
    # Applies the affine transformation described by the 4x4 matrix A to the nodes of self. Returns a copy.
    # Obviously this will need to be rethought if A[2,0:3] isn't [0,0,1].
    return ToolPath(np.matmul(A,np.vstack((self.nodes,[1]*self.nodes.shape[1])))[0:3,:],self.taxis)

  def PathToGCode(self,feedrate,fname):
    # Create a g-code file corresponding to a ToolPath object.
    # Note: this is currently tailored to use with Easel's controller, all coords relative to a starting
    # point manually jogged to.
  
    # CONSTANTS:
    atol = 0.01  #  10 microns; well below machine precision
    coordnames = "XYZ"
  
    progname = (os.path.splitext(fname)[0]+"xxxxx")[:5].upper()
    # There are many popular file extensions for g-code! Assume '.gcode' if not provided:
    # (But some others are ".nc", ".tap", ".mpt", ".mpf")
    if not "." in fname: fname += ".gcode"
    fidout = open(fname, 'w')
  
    ## HEADER:
    fidout.write("%\nO" + progname + "\nG17 G21 G40\n")
    # already at 0 # fidout.write("G00 X0. Y0.\n")
    # Easel automatically pulls bit up to about [0,0,5] before starting. Presumably it moves back to origin.
    # TEST THIS THOUGH!!!!!!!!!!!!!!!!!!
  
    ## PATHS:
    setfeedrate = False  #  it's not set yet. The machine should be told the rate on the 1st G01.
    for j in range(self.taxis.shape[1]-1):
      # taxis mode is taxis[1,j]
      # sequence of nodes whither to move is nodes[:,taxis[0,j]+1:taxis[0,j+1]+1]
      fidout.write(f"G0{self.taxis[1,j]} ")  #  set motion mode (cutting or rapid)
      for newpos in self.nodes[:,self.taxis[0,j]+1:self.taxis[0,j+1]+1].T: 
        fidout.write("X%.2f Y%.2f Z%.2f" % tuple(newpos))
        if self.taxis[1,j] and not setfeedrate:
          fidout.write(" F%d" % feedrate)
          setfeedrate = True
        fidout.write("\n")

    ## FOOTER:
    # fidout.write("G00 Z5.\nM30\n%%\n")  #  retract 5mm
    fidout.write("M30\n%%\n")
    fidout.close()

  def plot(self, ax, color='black', linewidth=1):
    # Plot a ToolPath on the 3D axis "ax", in the given color.
    for j in range(self.taxis.shape[1]-1):
      if self.taxis[1,j]:  #  only the cut paths, not the motions.
        ax.plot(*self.nodes[:,taxis[0,j]:taxis[0,j+1]+1], color, linewidth=linewidth)

def compileToolPath(paths, taxes):
  # Combine compatible lists of 3D paths and heights into a ToolPath object.
  # Each starts and stops from the origin (so they should be prepared individually as fairly large operations)
  nodes = np.hstack([np.zeros((3,1))]+paths)
  cumpathlengths = np.cumsum([path.shape[1] for path in paths])
  taxis = np.row_stack((np.hstack((0,cumpathlengths[:-1])),taxes))
  return ToolPath(nodes, taxis)

def catToolPaths(TPList):
  raise ToolPathError("This isn't yet updated for the new ToolPath definition.")
  # Joins the ToolPaths in TPList.
  indexoffsets = np.cumsum([0]+[tp.nodes.shape[1]-1 for tp in TPList[0:-1]])
  return ToolPath(np.hstack([tp.nodes[:,0:-1] for tp in TPList[0:-1]] + [TPList[-1].nodes]),
    np.vstack([np.hstack([tp.taxis[0,0:-1]+io for (tp,io) in zip(TPList[0:-1],indexoffsets)] +
      [TPList[-1].taxis[0,:]+indexoffsets[-1]]),
      np.hstack([tp.taxis[1,0:-1] for tp in TPList[0:-1]]+[TPList[-1].taxis[1,:]])]))

def thicknesser(xran, yran, zht, sht, offset, feedrate, fname):
  # Utility to plane a rectangle.
  # Cuts rectangle at height zht over rectangle [xran[0],xran[1]]x[yran[0],yran[1]].
  # Uses safe height sht.
  xwd = np.abs(np.diff(xran)[0])
  ywd = np.abs(np.diff(yran)[0])
  if True:  #  Condition ought to be "xwd<ywd:", but the Z & Y axes aren't quite perpendicular
    (xran,yran,xwd,ywd) = (yran,xran,ywd,xwd)
    swapped = True
  else:
    swapped = False
  numRows = 2 + int(np.floor(ywd/offset))
  y = np.linspace(yran[0],yran[1],num=numRows)
  nodes = np.vstack((np.pad(xran[[0,0,1,1]], (0,2*numRows-3), mode='wrap'),
    np.hstack((y[0],np.kron(y,[1,1]))), np.hstack((sht,np.repeat([zht],2*numRows)))))
  nodes = np.hstack((nodes, np.matmul(nodes[:,-1,None],np.array([[1,0]]))))
  nodes[2,-2:0] = sht
  #  print(nodes)  #  debugginh
  if swapped:
    nodes = nodes[[1,0,2],:]
  ToolPath(nodes, np.array([[0],[nodes.shape[1]]]), np.array([sht])).PathToGCode(feedrate, fname)

############################################################################################################

class PathGrid:
  # A PathGrid is a system of paths, each of constant y, with y increasing. Within each path, z varies as x
  # increases. This definition is not remotely isotropic, but the idea is to set up a PathGrid as a sequence
  # of parallel cuts, and then reorient it as necessary.

  # At some stage I should probably refactor this so that a PathGrid is a list of dictionaries rather than a
  # dictionary of arrays and lists.

  def __init__(self, y, xz):
    # Create a PathGrid with ordinates 'y' and corresponding paths 'xz.'

    # Error checking:
    if not isinstance(y, np.ndarray):
      raise PathGridError("PathGrid attribute y must be a Numpy NDArray.")
    if not all(isinstance(xzp, np.ndarray) for xzp in xz):
      raise PathGridError("PathGrid attribute xz must be a list of Numpy NDArrays.")
    if len(y)>1 and any (np.diff(y) <= 0):
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

  def maxoseg(self, segment):
    # returns scalar maximum z height over the segment from segment[:,0] to segment[:,1], interpolated
    if np.diff(segment[1,:]):
      try:
       return max([np.interp(np.interp(yp,segment[1,:],segment[0,:]),xzp[0,:],xzp[1,:])  #  self's hgt at intersec
         - np.interp(yp,segment[1,:],segment[2,:])  #  segment's height at intersection
         for (yp,xzp) in zip(self.y,self.xz) if yp>=min(segment[1,:]) and yp<=max(segment[1,:])])
      except:
       breakpoint()
    else:
      xzp = next(xzp for (yp,xzp) in zip(self.y,self.xz) if yp==segment[1,0])  #  assume existence & uniqueness
      # Find x values of nodes within segment and endpoint of segment:
      xChek=np.hstack((xzp[0,np.logical_and(xzp[0,:]>min(segment[0,:]),xzp[0,:]<max(segment[0,:]))],segment[0,:]))
      return np.max(np.interp(xChek,xzp[0,:],xzp[1,:])-np.interp(xChek,segment[0,:],segment[2,:]))

  def pospar(self):
    return PathGrid(self.y, [pospar(xzp) for xzp in self.xz])

  def min(self,other):
    return PathGrid(self.y, [minPLFs(xza,xzb) for (xza,xzb) in zip(self.xz,other.xz)])

  def whereBelow(self, other):
    # Trims a PathGrid to only the parts where self is below other.
    # Return structure is not a PathGrid; it's a list of continuous xyz cut-paths. (Each will have constant y.)
    dg = (other-self).pospar()  #  "difference grid"
    segl = [np.hstack((
      np.argwhere(np.diff(np.insert(np.logical_or(xzp[1,:-1]>1e-3,xzp[1,1:]>1e-3),0,False).astype(int))==1),
      2+np.argwhere(np.diff(np.append(np.logical_or(xzp[1,:-1]>1e-3,xzp[1,1:]>1e-3),False).astype(int))==-1)
      )) for xzp in dg.xz]
    xzll = [[yp,np.concatenate((
      np.array([[xzd[0,se[0]]],[np.interp(xzd[0,se[0]],xzs[0,:],xzs[1,:])]]),
      xzs[:,np.logical_and(xzs[0,:]>xzd[0,se[0]], xzs[0,:]<xzd[0,se[1]-1])],
      np.array([[xzd[0,se[1]-1]],[np.interp(xzd[0,se[1]-1],xzs[0,:],xzs[1,:])]])),1)]
      for (yp,seg,xzd,xzs) in zip(self.y,segl,dg.xz,self.xz) for se in seg]
    return ([np.row_stack((y_xz[1][0,:],np.full(y_xz[1].shape[1],y_xz[0]),y_xz[1][1,:])) for y_xz in xzll])

  def maxabs(self):
    return max([max(abs(xzp[1,:])) for xzp in self.xz])

  def nodeCount(self):
    return sum([xzp.shape[1] for xzp in self.xz])

  # Now overload arithmetic operators. Generally these apply to the z component. The y attribute won't be
  # altered, and a copy isn't necessary; there's no foreseeable use-case of a PathGrid being constructed and
  # used in calculations, and then having its y values changed.

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

  def castToMold(self, ballrad, ltol, floor=-10000):  #  (NINF was producing NINF-NINF warnings)
    # Converts "self," describing an actual cut surface, to a corresponding "mold," describing the lowest
    # surface the drillbit's spherical head's tip can traverse. Both are PathGrids. The drillbit's head's
    # radius is ballrad. The lateral tolerance for path equality is ltol. Don't cut below floor. If no floor
    # is needed, set floor to -10000.
  
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
            # indeX of points coincident with sphere under VErtical projection:
            vex = np.abs(xgrid - node[0]) < xzrad
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
      mxz[j] = fitPWL(np.vstack([xgrid,zgrid-ballrad]),ltol)  #  -ballrad since tracking tip, not centre
    return PathGrid(self.y, mxz)
  
  def moldToCast(self, ballrad, ltol):
    # Inverts the above function. Assumes no floor.
    return -((-self).castToMold(ballrad,ltol))

  def roundJoint(self,ballrad,ltol):
    # The idea is to minimally modify the PathGrid "self" so as to be tangible from above & below with a
    # ball-end drill bit with end radius "ballrad."
  
    # This works best under the condition that all parts of the original surface are tangible from above OR
    # below. Under that condition, we can round the concave & convex extremities in either order with
    # identical results.
    return self.castToMold(ballrad,ltol).moldToCast(2*ballrad,ltol).castToMold(ballrad,ltol)
  
  def plot(self, ax, color='black'):
    # Plot a PathGrid on the 3D axis "ax", in the given color.
    for (yp,xzp) in zip(self.y, self.xz):
      ax.plot(xzp[0,:], yp*np.ones(xzp.shape[1]), xzp[1,:], color, linewidth=1)

  def MultiToolGreedy(self, shpaid, toolSeq):
    # This creates a sequence of ToolPaths (one per tool in the provided sequence of ball-end drill bits, of
    # presumably decreasing radius), each working its way down from an initial stockheight, in general via a
    # sequence of not-very-deep cuts. Probably in practice the toolSeq will have length 1 (a single bit) or
    # 2 (a rough cut followed by a fine cut).
    # Inputs:
    #   self:        the target pathgrid
    #   shpaid:      initial stock height [PathGrid or scalar]
    #   toolSeq:     list of dictionaries defining (in cut order, coarse to fine):
    #     bitrad:      RADIUS (NOT DIAMETER) of bit, assumed to be ball-nose
    #     cude:        depth of wood it can remove in each pass
    #     ds:          ('downsample') offset for this tool will be ds times self's offset

    # Most of the logic here deals with the first tool in the list. The rest are handled via recursion.

    # Calculate this tool's target centre depth & resulting stock depth:
    # Minimum waste to leave for the final fine cut; may need tweaking:
    buffer = 0.75*toolSeq[-1]['cude'] if len(toolSeq)>1 else 0
    ballrad = toolSeq[0]['bitrad']
    target = (self+buffer).castToMold(ballrad, 0.1*ballrad)  #  target mold for this tool
    targetDS = target.downsample(toolSeq[0]['ds'])  #  downsampled for actual cutting
    toolFinalSH = targetDS.upsample(self).moldToCast(ballrad,0.1*ballrad)

    # There's probably a better way to deal with this, but at this point we'll convert shpaid to a PathGrid
    # if it's just a scalar. Note shpaid remains fullsampled.
    if isinstance(shpaid, numbers.Number):
      shpaid = PathGrid(self.y, [np.vstack((xzp[0,:],np.full(xzp.shape[1],shpaid))) for xzp in self.xz])
    shpaid = shpaid.castToMold(ballrad, 0.1*ballrad)  #  working in 'mold' view from now on.
    shpaidDS = shpaid.downsample(toolSeq[0]['ds'])

    # Now we want to build a ToolPath in a loop, but shouldn't do that directly, for memory-handling
    # reasons. So we build two lists:
    scpaths = []  #  a list of paths to be joined end-to-end, and
    taxes = []    #  a corresponding list of taxis modes (0 or 1)
    curpos = np.zeros((3,1))  #  just used to determine where to go next
    numcuts = int(np.ceil((shpaid-target).maxoxy()/toolSeq[0]['cude']))
    if numcuts<1:
      raise PathGridError("Wait, I'm confused. One tool's not cutting anything.")
    cutlevels = toolSeq[0]['cude'] * np.flip(np.arange(numcuts))
    PostSH = shpaidDS  #  effectively initialising PreSH [see 2 lines further]
    for cule in cutlevels:
      PreSH = PostSH
      PostSH = (shpaidDS+0.001).min(targetDS+cule)  #  added micron to clarify where PostSH is below shpaid
      cutxyz = PostSH.whereBelow(shpaidDS)  # segments needing cutting at this level
      # breakpoint()
      while len(cutxyz):
        dists = np.row_stack([[np.sum((curpos-seg[:,0,None])**2), np.sum((curpos-seg[:,-1,None])**2)]
          for seg in cutxyz])
        # Now find location in this two column array of the minimum; pop the rowth el of cutxyz; rev if 2nd col:
        locmin = np.unravel_index(dists.argmin(),dists.shape)  #  location of shortest (squared) distance
        try:
          nextpath = cutxyz.pop(locmin[0])  #  grab the closest path from the list
        except:
          breakpoint()
        if locmin[1]: nextpath = np.fliplr(nextpath)  #  reverse it if the end is closer than the start

        # Now we know where we're headed, but there may be stuff to avoid en route. Firstly, if it's close, let's
        # see whether we can go there in a single straight line:
        if dists[locmin] > 625 or PostSH.maxoseg(np.column_stack((curpos, nextpath[:,0]))) > 0:
          # There are obstacles.
          minZ = PreSH.maxoseg(np.array([[curpos[0,0],nextpath[0,0]],[curpos[1,0],nextpath[1,0]],[0,0]])) + 2
          # Rapidly rise to minZ, move x & y, then slowmo down to start of nextpath:
          scpaths.append(np.array([[curpos[0,0],nextpath[0,0]],[curpos[1,0],nextpath[1,0]],[minZ,minZ]]))
          taxes.append(0)  #  rapid mode
        scpaths.append(nextpath)
        taxes.append(1)  #  feed mode
        curpos = nextpath[:,-1,None]

    # Now also need to return to origin:
    taxes.append(0)
    if PostSH.maxoseg(np.column_stack((curpos, np.zeros((3,1))))) > -2:
      # There are obstacles.
      minZ = PostSH.maxoseg(np.array([[curpos[0,0],0],[curpos[1,0],0],[0,0]])) + 2
      scpaths.append(np.array([[curpos[0,0],0,0],[curpos[1,0],0,0],[minZ,minZ,0]]))
    else:
      scpaths.append(np.zeros((3,1)))

    # Return value, with recursive call if we aren't already looking at the finest tool:
    TPList = [compileToolPath(scpaths,taxes)]
    if len(toolSeq)>1:
      # 'Upsample'; determine current stock height from target depth:
      newStockHeight = target.upsample(self).moldToCast(ballrad,0.1*ballrad)
      # Call MultiToolPG recursively for remaining tools:
      TPList += self.MultiToolGreedy(toolFinalSH, toolSeq[1:])
    return TPList

  def downsample(self, ds):
    # Return a coarser PathGrid with every dsth row from self.
    if ds==1:
      return self
    else:
      retind = np.arange(int(ds/2),len(self.y),ds)
      return PathGrid(self.y[retind], [self.xz[k] for k in retind])

  def upsample(self, finer):
    # Return a finer PathGrid, every dsth row of which is from self. The new rows are at infinite height.
    # (Actually I've changed every infinity to 10 metres to dodge some "inf-inf" warnings.)
    moldxz = [np.vstack((xzp[0,:],np.full(xzp.shape[1],10000))) for xzp in finer.xz]
    for (yp,xzp) in zip(self.y, self.xz):
      yind = next(j for (j,ypf) in enumerate(finer.y) if ypf==yp)
      moldxz[yind] = xzp
    return PathGrid(finer.y, moldxz)

  def clipToRegion(self, reg):
    # Return a new PathGrid, just the parts of self that obey the condition in 'reg'.
    # Each region is assumed to be convex!!!

    # Go through the rows one by one, refine them to high resolution, mask them, and then renode. (Note that
    # the method used in pospar won't work because we don't have a linearity assumption in the region's
    # boundary.) For now I will hard-code the high resolution at 1mm, which should be fine if we always
    # overlap by at least that much.
    resol = 1
    newy = []
    newxz = []
    for (kx,(yp,xzp)) in enumerate(zip(self.y,self.xz)):
      x = np.sort(np.unique(np.concatenate((xzp[0,:]+resol,np.arange(xzp[0,0],xzp[0,-1],resol)))))
      x = x[reg(x,yp)]
      if len(x)>0:
        newy.append(yp)
        z = np.interp(x, xzp[0,:], xzp[1,:])
        newxz.append(fitPWL(np.vstack((x,z)),1e-3))  #  can afford fine tol, since mostly can reuse nodes

    return PathGrid(np.array(newy), newxz)

def addPLFs(a,b):
  # Inputs two lists of 2D numpy arrays of two rows each (x & z), increasing in x, describing each PLF
  # (piece-wise linear function).
  # Outputs another PLF, in the same format, being a+b.
  # Assumes a & b start & end at the same x values.
  x = np.sort(np.unique(np.concatenate((a[0,:],b[0,:]))))
  z = np.interp(x, a[0,:], a[1,:]) + np.interp(x, b[0,:], b[1,:])
  return np.vstack([x,z])
  
def subtractPLFs(a,b):
  # Inputs two lists of 2D numpy arrays of two rows each (x & z), increasing in x, describing each PLF
  # (piece-wise linear function).
  # Outputs another PLF, in the same format, being a-b.
  # Assumes a & b start & end at the same x values.
  x = np.sort(np.unique(np.concatenate((a[0,:],b[0,:]))))
  z = np.interp(x, a[0,:], a[1,:]) - np.interp(x, b[0,:], b[1,:])
  return np.vstack([x,z])
  
def pospar(plf):  #  "positive part"
  # Inputs a piecewise linear function (described by a 2-row Numpy array with the first row non-descending).
  # Outputs max(0,plf) in the same format & over the same range.

  px = plf[1,:]>0  #  indices of positive values
  csx = np.nonzero(np.diff(px))[0]  #  indices of segments crossing horizontal axis
  nodes = np.column_stack((
    [plf[:,px]] +
    list(((plf[0,cx]*plf[1,cx+1]-plf[0,cx+1]*plf[1,cx])/(plf[1,cx+1]-plf[1,cx]),0) for cx in csx) +
    [(np.array([[1],[0]])*plf[:,[0,-1]])[:, np.logical_not(px[[0,-1]])]]))
  return nodes[:,np.argsort(nodes[0,:])]

def minPLFs(a,b):  #  pointwise minimum of two piecewise linear functions
  # Inputs two piecewise linear functions (described by 2-row Numpy arrays with first row non-descending).
  # Outputs min(a,b) in the same format & over the same (assumed common) range.
  return subtractPLFs(a, pospar(subtractPLFs(a,b)))

############################################################################################################

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

############################################################################################################
####    DEFINE UTILITIES:    ###############################################################################
############################################################################################################

# Affine transformation matrices ("ATMs") should probably be a class ... but that seems like a lot of
# reinventing, and also I may not actually need much functionality.

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

############################################################################################################
# THE FUNCTIONS meanroundingup & fitPWL ARE ABOUT EFFICIENTLY TRACING A PATH WITH LINE SEGMENTS:

def meanroundingup(a,b):
  return(int(0.75+(a+b)/2))  #  mean of integers a & b, but rounds up if sum is odd.

def fitPWL(xy, ltol):
  # Inputs a 2xN array of plane points, assumed to be increasing in x, and a lateral tolerance, and outputs
  # a piecewise-linear approximation to the heights. Specifically, returns a subset of the entered nodes,
  # minimising node count subject to preserving the path to within a tolerance.

  # Use a greedy method to determine breakpoints. (Should then be refined, but maybe later.)
  N = xy.shape[1]
  if N<2: return xy
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
  return xy[:,keepNodes]
