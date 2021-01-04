# create g-code for clamps for dode

import cnc
import sys
import numpy as np

# This creates a gcode file to cut a clamping jig layer (presumably out of mdf).
if len(sys.argv) != 3:  #  number of parameters, + 1 for "dodeClamps.py"
  print("Usage: python dodeClamps.py rf mdfDepth hiZ ballrad offset")
  print("  rf       = planar exradius of exterior faces")
  print("  mdfDepth = depth of the stock we're cutting the jig layer from")
  print("  hiZ      = upper end of the range of Z coordinates this layer will clamp")
  print("  ballrad  = radius of ball-end bit cutting the clamp")
  print("  offset   = offset between cuts")
  print("If mdfDepth>hiZ then it's an end clamp, so cut bottom as well as sides."
else:
  do the thing

##################################################################################################################
#####     APPENDIX A:  VERTICES' COORDINATES     #################################################################
##################################################################################################################
"""
This script creates gcode files for cutting jigs for a dodecahedral box. The size parameter we're using, rf,
is defined as the planar exradius of the exterior faces, i.e., the distance from the centre of a face to a vertex
of that face. To obtain coordinates of the exterior vertices, multiply the following vectors by rf/2, and note
that phi denotes the "golden ratio" (1+sqrt(5))/2:

BOTTOM FACE:
          phi                 -phi+1              -2                  -phi+1              phi
          sqrt(3-phi)         sqrt(phi+2)         0                   -sqrt(phi+2)        -sqrt(3-phi)
          -phi-1              -phi-1              -phi-1              -phi-1              -phi-1

'TROPIC OF CAPRICORN':
          phi+1               -1                  -2*phi              -1                  phi+1
          sqrt(phi+2)         sqrt(4*phi+3)       0                   -sqrt(4*phi+3)      -sqrt(phi+2)
          -phi+1              -phi+1              -phi+1              -phi+1              -phi+1

'TROPIC OF CANCER':
2*phi               1                   -phi-1              -phi-1              1
0                   sqrt(4*phi+3)       sqrt(phi+2)         -sqrt(phi+2)        -sqrt(4*phi+3)
phi-1               phi-1               phi-1               phi-1               phi-1

TOP FACE:
2                   phi-1               -phi                -phi                phi-1
0                   sqrt(phi+2)         sqrt(3-phi)         -sqrt(3-phi)        -sqrt(phi+2)
phi+1               phi+1               phi+1               phi+1               phi+1
"""
