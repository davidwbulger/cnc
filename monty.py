# Load the Monty stl exported from Blender, and convert it to g-code.

import cnc
import numpy as np

pt = cnc.polyTri("./Monty.stl")

# Blender has measured in metres, whereas we're using mm, so need to multiply
# everything by 1000:
pt = pt.afxform(np.diag([1000,1000,1000,1]))
# print(pt.bbox())

if False:  #  This was all good, but we're just doing spots now.
  # Do the rough cut:
  pg = pt.toPG(1.0)
  tp = pg.MultiToolGreedy(0, [dict(bitrad=3.0,cude=3.0,ds=1)])[0]
  # tp = pg.SingleToolNoOpt(3.0)
  tp.PathToGCode(1200, "monty_rough.gcode")

  # Do spots!
  pg = pt.toPG(0.3, yrange=[22,68])
  tp = pg.SingleToolNoOpt(1.0)
  tp.PathToGCode(1200, "monty_spot.gcode")
  
  pg = pt.toPG(0.3, xrange=[170,215], yrange=[24,36])
  tp = pg.SingleToolNoOpt(1.0)
  tp.PathToGCode(1200, "monty_mast.gcode")

# Do the fine cut:
pg = pt.toPG(0.3)
tp = pg.SingleToolNoOpt(1.0)
tp.PathToGCode(1800, "monty_fine.gcode")
