# cnc
Tools for controlling a CNC machine

A python library & associated code for controlling a grbl CNC machine. I'm not at this stage expecting any of this to be of interest to anyone else, but who knows. Essentially, I want more direct control and legal freedom than available tools provide, so I am creating this python framework to produce gcode.

Files:
* cnc.py: the library.
* cncTests.py: some unit tests. I've been using this during development and have of course redesigned elements of cnc.py during that development, so most of these tests might not currently work.
* smallDode.py: a script that uses the library to create gcode to cut a certain shape.
* Star.gcode: some earlier experiment. I think this was manually written gcode to cut a 2D star.
* tet.scad: an OpenSCAD program to create a tetrahedron.
* tet.stl: the tetrahedron made by tet.scad, output as an "stl" file (i.e., triangulated surface).
