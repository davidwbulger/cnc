# create g-code for clamps for dode

import cnc
import numpy as np

# This creates a gcode file to cut a clamping jig layer (presumably out of mdf).

# PARAMETERS:
cude = 4  #  safe cutting depth for a single pass
offset = 2
feedrate = 800
stockdepth = 20
rf = 67
parity = 0


# 

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
