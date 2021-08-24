# Load the Monty stl exported from Blender, and convert it to g-code.

import cnc
import numpy as np

pt = cnc.polyTri("./Monty.stl")

# Blender has measured in metres, whereas we're using mm, so need to multiply
# everything by 1000:
pt = pt.afxform(np.diag([1000,1000,1000,1]))
# print(pt.bbox())

# Do the rough cut:
pg = pt.toPG(0.6)
tp = pg.MultiToolGreedy(0, [dict(bitrad=1,cude=2,ds=1)])[0]
tp.PathToGCode(1200, "monty_rough.gcode")

# Do the fine cut:
pg = pt.toPG(0.3)
tp = pg.SingleToolNoOpt(0.75)
tp.PathToGCode(1200, "monty_fine.gcode")
