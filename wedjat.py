from matplotlib import pyplot as plt
import numpy as np

from scipy import interpolate
"""
Image tool ideas for inlays:
  Given a B&W raster image:
    compute the boundaries,
    compute offset boundaries
  Repeat until no change:
    replace black with union of circles entirely contained within black
    replace white with union of circles entirely contained within white
"""

factor = 2
LW = 24*factor  #  linewidth
(fig,ax) = plt.subplots(figsize=[6.4*factor,4.8*factor])

# Iris
theta = np.linspace(0,2*np.pi,361)
# ax.plot(np.cos(theta), np.sin(theta), color='r', linewidth=1)
ax.fill(np.cos(theta), np.sin(theta), color='grey')

# Lids
# First define exact, angular path, before smoothing:
angle = 7/12*np.pi
lowth = 0.17  #  in units of iris radius
width = 3.2
tail = 0.8
nc = 17  #  node count for splines
theta = np.linspace(np.pi/2-angle/2,np.pi/2+angle/2,2*nc,endpoint=False)[1:]
lidrad = width/(2*np.sin(angle/2))
uplid = np.array([lidrad*np.cos(theta),lidrad*(np.sin(theta)-np.cos(angle/2))])
# Firstly erase top part of iris:
ax.fill(np.concatenate((uplid[0],uplid[0][::-1])),
  np.concatenate((uplid[1]-lowth,uplid[1][::-1]-lowth+0.4)), color='white')
xy = np.concatenate((
  np.array([np.linspace(width/2+tail,width/2,nc,endpoint=False),
    np.full(nc,-lowth)]),
  uplid - [[0],[lowth]],
  uplid[:,::-1]*[[1],[-1]] - [[0],[lowth]],
  [[width/2],[-lowth]]
  ), axis=1)
(tck,u) = interpolate.splprep(xy, s=0)
unew = np.arange(0, 1.001, 0.001)
draw = interpolate.splev(unew, tck)
ax.plot(draw[0], draw[1], color='k', linewidth=LW, solid_capstyle='round'
  )[0].set_antialiased(False)

# Brow
lift = 0.6
brow = uplid*(lift+lidrad)/lidrad
brow[1,:] += lift - np.max(brow[1,:]) + np.max(uplid[1,:]) - lowth
corner = np.array([[width/2],[-lowth]]) +lift*np.array([[np.tan(angle/4)],[1]])
brow = brow[:,(brow[0,:]>-width/2) & (brow[1,:]>corner[1])]
xy2 = np.concatenate((
  np.array([np.linspace(width/2+tail,corner[0,0],nc,endpoint=False),
  np.full(nc,lift-lowth)]),
  brow
  ), axis=1)
(tck,u) = interpolate.splprep(xy2, s=0)
unew = np.arange(0, 1.001, 0.001)
draw = interpolate.splev(unew, tck)
ax.plot(draw[0], draw[1], color='k', linewidth=LW, solid_capstyle='round'
  )[0].set_antialiased(False)

# print(np.concatenate((xy[:,[0]],xy2[:,[0]]),axis=1))
# print([np.max(xy[1,:]),np.max(xy2[1,:])])

LLang = np.pi*4/3    #  240deg
RLang = np.pi*17/12  #  255deg

pt = xy[:,[np.argmax([np.cos(LLang),np.sin(LLang)]@xy)]]
LL = pt + [[0,0,-0.1],[0,-0.9,-0.5]]
ax.plot(LL[0], LL[1], color='k', linewidth=LW, solid_capstyle='round'
  )[0].set_antialiased(False)

pt = xy[:,[np.argmax([np.cos(RLang),np.sin(RLang)]@xy)]]
floor = np.min(LL[1,:])
centre = LL[1,-1]
theta = np.linspace(np.pi*5/6,np.pi*8,1000)
RL = np.array([np.cos(theta), np.sin(theta)])
RL *= np.exp(-0.75*theta)
RL[1,:] *= (floor-centre)/np.min(RL[1,:])
RL[1,:] += centre - RL[1,-1]
RL = RL[:,RL[1,:]<=pt[1]]
RL[0,:] *= (width/2+tail-pt[0])/(np.max(RL[0,:])-np.min(RL[0,:]))
RL[0,:] += pt[0]-RL[0,0]
ax.plot(RL[0], RL[1], color='k', linewidth=LW, solid_capstyle='round'
  )[0].set_antialiased(False)

ax.axis('equal')
ax.axis('off')
plt.savefig("wedjat.png")
plt.show()
