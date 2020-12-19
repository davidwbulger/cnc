# Load an stl file designed for 3D printing (downloaded from the internet) and manipulate it to produce code for
# CNC carving.

import cnc
import numpy as np
import matplotlib.pyplot as plt

##################################################################################################################
# PARAMETERS:

##################################################################################################################
# LOAD AND MANIPULATE:

pt = cnc.polyTri("./Full_Bus.stl")
pg = pt.toPG(0.1)

##################################################################################################################
# PLOT:

fig = plt.figure()
ax = fig.add_subplot(projection="3d")
pg.plot(ax)
mng = plt.get_current_fig_manager()
mng.window.state("zoomed")
plt.show()

##################################################################################################################
# WRITE THE GCODE:

