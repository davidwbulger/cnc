import matplotlib.pyplot as plt
import numpy as np
import re

# Ad-hoc script to extract and plot points from gcode.

pat = re.compile("Y([-\d\.]+) Z([-\d\.]+)")
yz = np.array([list(map(float,m.group(1,2)))
  for l in open("steamback1.8.gcode") for m in re.finditer(pat,l)])

(fig,ax) = plt.subplots()
ax.scatter(yz[:,0], yz[:,1], color='r')
ax.set_xlabel('Y')
ax.set_ylabel('Z')
ax.set_title('Range of Z per Y')
plt.show()
