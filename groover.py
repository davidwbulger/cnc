import cnc
import numpy as np

# Add grooves to the iris/album.

radius = 7
depth = 0.3
offset = 0.21387  #  mm per revolution
tracks = np.array([[23.3, 30.5], [31.4, 35.9], [36.8, 45.8], [46.7, 51.2],
  [52.1, 56.1]])

rpm = 2*np.pi/offset  #  radians per millimetre
angleRanges = [np.arange(d[0]*rpm, d[1]*rpm, 0.01) for d in tracks]
paths = [np.array([np.cos(ar)*offset*ar, -np.sin(ar)*offset*ar, 0*ar-depth]).T
  for ar in angleRanges]
cnc.cutSequence(paths, 1000, 3, "grooves.gcode")