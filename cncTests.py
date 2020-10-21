# Tests for the module "cnc"

import cnc
import numpy as np
import matplotlib.pyplot as plt

######################################################################################################################################
# DEFINE UTILITIES:

def plotseglist(seglist):
  for seg in seglist:
    xy = seg['xy']
    K = seg['K']
    x = np.linspace(xy[0,0],xy[0,-1],num=200)
    y = cnc.calcArc(xy,K,x)
    plt.plot(x,y,color='red')
  xy = np.concatenate([seglist[0]['xy'][:,[0]],np.concatenate([seg['xy'][:,[1]] for seg in seglist],axis=1)],axis=1)
  plt.plot(xy[0],xy[1],'ro')

######################################################################################################################################
# DEFINE TESTS:

def test_constructor():
  # A function to test the constructor of the PolyTri class. It has two modes of operation;
  # it can load an ASCII-format .stl file:
  eph = cnc.PolyTri('tet.stl')
  print(eph)
  print(eph.facets.shape)
  print(eph.facets.strides)

  # or it can load a binary-format .stl file:
  eph = cnc.PolyTri('Tetra_v1.stl')
  print(eph)
  print(eph.facets.shape)
  print(eph.facets.strides)

def test_findK():
  K = cnc.findK(np.array([[0,500,1000],[0,0.5,0]]),0.1)
  print(K)

def test_fitTL():
  # Create a random polynomial, & then approximate it with arcs:
  x = np.linspace(0,1,num=201)
  rts = np.random.rand(5,1)
  y = np.prod(x-rts, axis=0)
  y = y * 0.6 / np.max(np.abs(y))

  xy = np.concatenate(([x],[y]))
  vtol = 0.0001
  seglist = cnc.fitTL(xy,vtol)
  print(len(seglist))

  plt.plot(x,y,color='black')
  plotseglist(seglist)
  plt.axis('equal')
  plt.show()

def test_fitPWL():
  # Create a random polynomial, & then approximate it with arcs:
  x = np.linspace(0,1,num=2001)
  rts = np.random.rand(5,1)
  y = np.prod(x-rts, axis=0)
  y = y * 0.6 / np.max(np.abs(y))

  xy = np.concatenate(([x],[y]))
  ltol = 0.001
  xyout = cnc.fitPWL(xy,ltol)
  print(xyout.shape[1])

  plt.plot(x,y,color='black')
  plt.plot(xyout[0],xyout[1],'ro')
  plt.axis('equal')
  plt.show()

def test_PathsToGCode():
  # Create a star-shape gcode file to test PathsToGCode:
  phi = np.array([np.linspace(0,2*np.pi,11)]).T
  r = np.cos(2.5*phi)**2 + np.sqrt((7-3*np.sqrt(5))/2) * np.sin(2.5*phi)**2
  d = 3*np.sin(phi)**2
  verts = np.concatenate((40*r*np.cos(phi), 40*r*np.sin(phi), d), axis=1)
  path = {'fr':300, 'itin':verts}
  cnc.PathsToGCode([path], "Star.gcode")

def test_avoid(actuallyLSZ=False):
  # Tests "avoid" & plots the result, unless the argument is true: then tests "lsz" instead
  rng = np.random.default_rng()
  ry = rng.random()*10-5
  segment = np.array([[rng.random()*10-5,ry,rng.random()*10-5], [rng.random()*10-5,ry,rng.random()*10-5]])
  xy = np.array([rng.random()*10-5,rng.random()*10-5])
  ballrad = 2+2*rng.random()
  if actuallyLSZ:
    avret = cnc.lsz(segment, xy, ballrad)
  else:
    avret = cnc.avoid(segment, xy, ballrad)
  print(avret)

  # Also plot it!
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d') 
  ax.plot(segment[:,0], segment[:,1], segment[:,2], 'b', linewidth=3)
  phi = np.linspace(0,2*np.pi,25)
  cylx = np.outer([1,1], xy[0] + ballrad * np.cos(phi))
  cyly = np.outer([1,1], xy[1] + ballrad * np.sin(phi))
  cylz = np.outer([-5,5], np.ones(25))
  ax.plot_wireframe(cylx, cyly, cylz, color='red', linewidth=0.5)
  if np.all(np.isfinite(avret)):
    phi = np.outer([1]*13, phi)
    theta = np.outer(np.linspace(0,np.pi,13), [1]*25)
    sphx = ballrad * np.sin(theta) * np.cos(phi) + xy[0]
    sphy = ballrad * np.sin(theta) * np.sin(phi) + xy[1]
    sphz = ballrad * np.cos(theta)
    if actuallyLSZ:
      ax.plot_wireframe(sphx, sphy, sphz + avret, color='green', linewidth=1)
    else:
      ax.plot_wireframe(sphx, sphy, sphz + avret[0], color='black', linewidth=1)
      ax.plot_wireframe(sphx, sphy, sphz + avret[1], color='green', linewidth=1)
  ax.set_xlim(-5,5)
  ax.set_ylim(-5,5)
  ax.set_zlim(-5,5)
  ax.set_box_aspect((10,10,10))
  plt.show()

def makePaid():
  # A test function to make a "PathGrid" object for testing castToMold, widen & any similar functions yet to come.
  y = np.linspace(-0.05,2.05,22)
  xz = [np.array([[0,4],[1,1]])] + [np.array([[0,0,1,2,2,3,4],[1,2,2,1,2,2,1]])]*10 + \
    [np.array([[0,1,2,2,3,4,4],[1,2,2,1,2,2,1]])]*10 + [np.array([[0,4],[1,1]])]
  return {'y':y, 'xz':xz}

def testPaidUtes():
  paid = makePaid()
  # mold = cnc.castToMold(paid, 0.4, 0.002, 0)
  mold = cnc.roundJoint(paid, 0.4, 0.001)

  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d') 
  cnc.plot(paid, ax, 'red')
  cnc.plot(mold, ax, 'blue')
  plt.show()

def test_whereBelow():
  y = np.linspace(0,1,26)
  x = np.linspace(-1,1,101)
  pga = cnc.PathGrid(y, [np.vstack((x,yp*(1-x**2))) for yp in y])
  pgb = cnc.PathGrid(y, [np.vstack((x,1-yp*(1-x**2))) for yp in y])

  cutlist = pga.whereBelow(pgb)

  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d') 
  for (yp,clp) in zip(y,cutlist):
    for cl in clp:
      ax.plot(cl[0,:],y*np.ones(cl.shape[1]),cl[1,:],'black',linewidth=2)
  pga.plot(ax, 'red')
  pgb.plot(ax, 'blue')
  plt.show()

######################################################################################################################################
# RUN SOMETHING:

# test_constructor()
# test_findK()
# test_fitTL()
# test_PathsToGCode()
# test_avoid(True)
# test_fitPWL()
# testPaidUtes()

# 20 October 2020: the above was all quite a while ago & mightn't still work.
test_whereBelow()
