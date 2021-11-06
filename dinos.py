# Load the .stl file for Merry Cretaceousmas exported from Blender, and create
# Gcode to cut it.

import cnc

pt = cnc.polyTri("./Dinos.stl")

# Convert from m to mm:
pt = pt.afxform([[1e3,0,0,0],[0,1e3,0,0],[0,0,1e3,0],[0,0,0,1]])
print(pt)

pg = pt.toPG(0.3)

tp = pg.SingleToolNoOpt(1)
tp.PathToGCode(1000, "dinos.gcode")
