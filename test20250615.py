# Just a simple test. Goals are:
  # start to remind myself how to use my cnc.py library
  # test set-up after
    # a two-year hiatus
    # migration to gSender
# This just cuts a square, in two passes. My intention is just to "cut air."

import cnc
import numpy as np

(x,y) = np.array([[0, 0], [0,40], [40, 40], [40, 0], [0,0]]).T

# cnc.cutPath(x, y, 3, 2, 300, 3, "testSquare.gcode")
cnc.cutPath(x, y, 13, 2, 300, -7, "testSquare.gcode", 1.5)
# START FROM 10mm ABOVE SURFACE!!!
