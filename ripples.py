# create g-code for test run

import cnc
import numpy as np
import matplotlib.pyplot as plt

# PARAMETERS:
radius = 30
depth = 2
dpc = 0.5  #  depth per cut
frate = 150
offset = 0.7
nripples = 11.5  #  should be half-integer

# ADJUSTMENTS:
nrows = int(np.ceil(2*radius/offset)) - 1

# CREATE SINGLE GRID:
y = np.linspace(offset-radius, radius-offset, num=nrows)
x = [np.arange(-np.sqrt(radius*radius-yp*yp), np.sqrt(radius*radius-yp*yp), offset) for yp in y]
# print([len(xp) for xp in x])
dist = [np.sqrt((xp-radius)**2+yp**2)*nripples*np.pi/radius for (xp,yp) in zip(x,y)]  #  in radians of ripples
z = [(np.cos(distp)-1)*depth/2 for distp in dist]
# z = [(np.cos(distp)-1)*depth/2 + depth-dpc for distp in dist]
xz = [np.vstack((xp,zp)) for (xp,zp) in zip(x,z)]
taPG = cnc.PathGrid(y,xz)

# CONVERT TO TOOLPATH:
# Note that, currently, PacePathGrid is the only way to convert a PathGrid to a ToolPath. This code 'manually' handles the layers,
# so we'll start by pretending that the safe cut per depth is huge to get the target ToolPath:
taTP = taPG.pacePathGrid(0, 2, depth+1)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d') 
taTP.plot(ax, 'red')
cnc.hackaspect(ax)
plt.show()

# LAYERS:
# This would automate the whole layers thing:
# layers = cnc.catToolPaths([taTP.afxform(np.array([[1,0,0,0],[0,1,0,0],[0,0,1,zadj],[0,0,0,1]])) for zadj in np.flip(np.arange(0,depth,dpc))])
# layers.PathToGCode(frate, "ripples.gcode")

# But instead I'll just do the one cut & manually jog & repeat.
taTP.PathToGCode(frate, "ripples.gcode")
