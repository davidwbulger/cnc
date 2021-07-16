# create Blender model for big dode

import numpy as np
import bpy
import bmesh
C = bpy.context
D = bpy.data

# Working toward code for big dodecahedron. Currently at the design stage.

##  PARAMETERS  ###############################################################

# # COME BACK TO ALL OF THIS LATER:
# # Parameters for tenon & mortise shape:
# edgePropJoined = 0.9  #  max proportion of edge used for finger joint
# mitreMargin = 0.001  #  vertically, for simplicity
# cutsPerTooth = 6
# if cutsPerTooth%4 != 2:
#   raise ValueError("cutsPerTooth must be twice an odd number.")
# # This should be 2ce the bit radius, but really the bit carves slightly more:
# mortiseWidth = 0.0023
# toothWidth = 0.0021
# # Make the tooth "shorter" (in the direction normal to the "joint plane"
# # described by xzflat) by this amount:
# brachidont = 0.0005
# toothRoundingFactor = 1.9  #  Dunno why 1 doesn't work; this seems better.

# Main programs (SDFace#.gcode) stuff:
boxheight = 0.400  #  m
thm = 0.013  #  thickness of main exterior boards
thp = 0.009  #  thickness of plywood boards
thg = 0.002  #  thickness of gaps, e.g., around drawers (smallest straight bit)
mtd = 0.0055  #  depth of T-joint cuts into main boards
ptd = 0.0035  #  depth of T-joint cuts into plywood
ihwnd = 0.025  #  inner half-width of necklace drawer (gap between vert plies)

##  CALCULATED VALUES  ########################################################

phi = (1+np.sqrt(5))/2  #  the golden ratio
rf = boxheight/(2*phi+2)  #  planar exradius of face at outside

# heights of horizontal plies, at their tops:
capri = boxheight*(1.5-phi)  #  height of capricorn vertex
# hpts = np.concatenate((np.linspace(thm-0.5*boxheight,capri,2,False)[1:],
#   np.linspace(capri,0.5*boxheight-thm+thp,3,False)))
hpts = np.concatenate((np.linspace(thm-0.5*boxheight,capri,3,True)[1:],
  np.linspace(-capri,0.5*boxheight-thm+thp,2,False)))

##  UTILITIES  ################################################################

def MakeAndLinkMesh(name,meshdata):
  the_ob = MakeMesh(name,meshdata)
  C.collection.objects.link(the_ob)
  return the_ob

def MakeMesh(name,meshdata):
  the_mesh = D.meshes.new(name)
  the_mesh.from_pydata([tuple(v) for v in meshdata[0]],meshdata[1],meshdata[2])
  the_mesh.update()
  the_ob = D.objects.new(name, the_mesh)
  return the_ob

def BoolEq(mob1, mob2, boolop):
  # alter the mesh object mob1 to be mob1 bop mob2, for some operator bop.
  # (Thus named to suggest C's "assignment operators.")
  # Options for boolop are INTERSECT, UNION and DIFFERENCE.
  mod = mob1.modifiers.new(type="BOOLEAN", name="boolop")
  mod.object = mob2
  mod.operation = boolop
  C.view_layer.objects.active = mob1
  bpy.ops.object.modifier_apply(modifier=mod.name)

def DiffEq(mob1, mob2):
  # alter the mesh object mob1 to be mob1-mob2, "-" indicating boolean diff.
  BoolEq(mob1, mob2, "DIFFERENCE")

def AndEq(mob1, mob2):
  # alter the mesh object mob1 to be mob1^mob2, "^" indicating intersection.
  BoolEq(mob1, mob2, "INTERSECT")

def lambdaJoint(stileob, railob, margin, cude):
  # Helps to construct a "lambda" joint (an oblique T-joint, named for the
  # shape of a lowercase lambda).
  # stileob: mesh object representing the "stile" (receiving piece)
  # railob: mesh object representing the "rail" (inserted piece)
  # margin: width uncut around the edge of each stile face
  # cude: cut depth

  # This probably indicates a poor understanding of the available
  # functionality, but a lot of the logic in this function is simply pushing
  # shapes back & forth between "mesh" and "bmesh" format, since I know how to
  # iterate through the faces of, and extrude from, a bmesh, and I know how to
  # apply Boolean operations to a mesh.

  bmstile=bmesh.new()
  bmstile.from_mesh(stileob.data)
  #for sf in bmstile.faces[19:20]:
  for sf in bmstile.faces:
    # Obtain this face as a separate Bmesh object:
    mface = D.meshes.new("face")
    L = len(sf.verts)
    mface.from_pydata([v.co for v in sf.verts],
      [(j,j+1) for j in range(L-1)] + [(L-1,0)], [tuple(range(L))])
    bmface = bmesh.new()
    bmface.from_mesh(mface)
    origverts = set(bmface.verts)
  
    # Inset face:
    infa=bmesh.ops.inset_individual(bmface,faces=bmface.faces,thickness=margin)
  
    # Get just the interior, inset face, and convert it to a mesh:
    insfa = next(f for f in bmface.faces if not (set(f.verts) & origverts))
    insface = D.meshes.new("face")
    L = len(insfa.verts)
    insface.from_pydata([v.co for v in insfa.verts],
      [(j,j+1) for j in range(L-1)] + [(L-1,0)], [tuple(range(L))])
    cutob = D.objects.new("Inset", insface)
    C.collection.objects.link(cutob)
  
    # Now we want the intersection of cutob with railob. Unfortunately, if
    # there is a direct way in Blender of calculating the intersection of a
    # "watertight manifold" and a face, I haven't found it. So we'll extrude
    # cutob to make it 3D, then intersect, then discard all but the subfaces
    # of mface.
    cutbm = bmesh.new()
    cutbm.from_mesh(cutob.data)
    r = bmesh.ops.extrude_face_region(cutbm, geom=cutbm.faces)
    verts = [e for e in r['geom'] if isinstance(e, bmesh.types.BMVert)]
    bmface.faces.ensure_lookup_table()
    bmface.verts.ensure_lookup_table()
    bmesh.ops.translate(cutbm, vec=-0.1*cude*bmface.faces[0].normal,
      verts=verts)
    bmesh.ops.recalc_face_normals(cutbm, faces=cutbm.faces) # shouldn't need!
    cutbm.to_mesh(cutob.data)

    mod = cutob.modifiers.new(type="BOOLEAN", name="boolop")
    mod.object = railob
    mod.operation = "INTERSECT"
    C.view_layer.objects.active = cutob
    bpy.ops.object.modifier_apply(modifier=mod.name)
  
    cutbm = bmesh.new()  #  needed to get a clean slate
    cutbm.from_mesh(cutob.data)
    bmesh.ops.delete(cutbm,
      geom=[v for v in cutbm.verts if abs((v.co-bmface.verts[0].co) @
      bmface.faces[0].normal)>1e-4])
    bmesh.ops.translate(cutbm, vec=0.02*cude*bmface.faces[0].normal,
      verts=cutbm.verts)
    r = bmesh.ops.extrude_face_region(cutbm, geom=cutbm.faces)
    verts = [e for e in r['geom'] if isinstance(e, bmesh.types.BMVert)]
    bmesh.ops.translate(cutbm, vec=-1.02*cude*bmface.faces[0].normal,
      verts=verts)
    bmesh.ops.recalc_face_normals(cutbm, faces=cutbm.faces) # shouldn't need!
    cutbm.to_mesh(cutob.data)
  
    # Intersect this with Rail, and then remove that intersection from Stile:
    mod = cutob.modifiers.new(type="BOOLEAN", name="boolop")
    mod.object = railob
    mod.operation = "INTERSECT"
    C.view_layer.objects.active = cutob
    bpy.ops.object.modifier_apply(modifier=mod.name)
    mod = stileob.modifiers.new(type="BOOLEAN", name="boolop")
    mod.object = cutob
    mod.operation = "DIFFERENCE"
    C.view_layer.objects.active = stileob
    bpy.ops.object.modifier_apply(modifier=mod.name)
    D.objects.remove(cutob, do_unlink=True)
    
  # Now cut any remaining overlap away from the Rail:
  mod = railob.modifiers.new(type="BOOLEAN", name="boolop")
  mod.object = stileob
  mod.operation = "DIFFERENCE"
  mod.use_self=True
  C.view_layer.objects.active = railob
  bpy.ops.object.modifier_apply(modifier=mod.name)
    
def RegDodeGeom():
  # Determine the vertex coordinates and edge and face structures for a regular
  # dodecahedron. Scaled to have a height of 1.

  # Vertices:
  V = np.array([[np.sqrt(3-phi), phi, -phi-1],  #  bottom face, back R vertex
    [np.sqrt(phi+2), -phi+1, -phi-1],           #  bottom face, front R vertex
    [0, -2, -phi-1],                            #  bottom face, front C vertex
    [-np.sqrt(phi+2), -phi+1, -phi-1],          #  bottom face, front L vertex
    [-np.sqrt(3-phi), phi, -phi-1],             #  bottom face, back L vertex
    [np.sqrt(phi+2), phi+1, -phi+1],            #  capricorn, back R vertex
    [np.sqrt(4*phi+3), -1, -phi+1],             #  capricorn, front R vertex
    [0, -2*phi, -phi+1],                        #  capricorn, front C vertex
    [-np.sqrt(4*phi+3), -1, -phi+1],            #  capricorn, front L vertex
    [-np.sqrt(phi+2), phi+1, -phi+1],           #  capricorn, back L vertex
    [-np.sqrt(phi+2), -phi-1, phi-1],           #  cancer, front L vertex
    [-np.sqrt(4*phi+3), 1, phi-1],              #  cancer, back L vertex
    [0, 2*phi, phi-1],                          #  cancer, back C vertex
    [np.sqrt(4*phi+3), 1, phi-1],               #  cancer, back R vertex
    [np.sqrt(phi+2), -phi-1, phi-1],            #  cancer, front R vertex
    [-np.sqrt(3-phi), -phi, phi+1],             #  top face, front L vertex
    [-np.sqrt(phi+2), phi-1, phi+1],            #  top face, back L vertex
    [0, 2, phi+1],                              #  top face, back C vertex
    [np.sqrt(phi+2), phi-1, phi+1],             #  top face, back R vertex
    [np.sqrt(3-phi), -phi, phi+1]])             #  top face, front R vertex
  V = V*(1-0.5*phi)  #  normalise to height of 1
  
  # Edges and faces (won't change, so output as tuples):
  E = [(0,1),(1,2),(2,3),(3,4),(4,0),           #  bottom face
    (0,5),(1,6),(2,7),(3,8),(4,9),              #  lower longitudes
    (12,5),(5,13),(13,6),(6,14),(14,7),         #  right tropicals
    (7,10),(10,8),(8,11),(11,9),(9,12),         #  left tropicals
    (10,15),(11,16),(12,17),(13,18),(14,19),    #  upper longitudes
    (15,16),(16,17),(17,18),(18,19),(19,15)]    #  top face
  F = [(0,1,2,3,4),           #  bottom face
    (0,5,13,6,1),             #  lower back R face
    (1,6,14,7,2),             #  lower front R face
    (2,7,10,8,3),             #  lower front L face
    (3,8,11,9,4),             #  lower back L face
    (4,9,12,5,0),             #  lower back C face
    (5,12,17,18,13),          #  upper back R face
    (6,13,18,19,14),          #  upper front R face
    (7,14,19,15,10),          #  upper front C face
    (8,10,15,16,11),          #  upper front L face
    (9,11,16,17,12),          #  upper back L face
    (15,19,18,17,16)]         #  top face
  
  return (V,E,F)

###############################################################################
##  START CREATING THE SEPARATE OBJECTS  ######################################
###############################################################################

##  REMOVE EVERYTHING  ########################################################
for c in D.collections:
  D.collections.remove(c)
for c in D.objects:
  D.objects.remove(c)

##  CREATE OUTER & INNER SURFACES OF DODECAHEDRAL BOX  ########################

(V,E,F) = RegDodeGeom()  #  coordinates & connectivity of unit height dodec
doobj = MakeAndLinkMesh("DodeOut", (boxheight*V,E,F))  #  outer surface
diobj = MakeMesh("DodeIn", ((boxheight-2*thm)*V,E,F))  #  inner surface
DiffEq(doobj, diobj)  #  hollow out larger dode by removing smaller one
D.objects.remove(diobj, do_unlink=True)
doobj.hide_set(True)

##  CREATE VERTICAL AND HORIZONTAL PLYWOOD STRUTS  ############################
C.view_layer.active_layer_collection = C.layer_collection
toJoin = []
for m in [-1,1]:
  bpy.ops.mesh.primitive_cube_add(1, location=(m*(ihwnd+0.5*thp),0,0),
    scale=(thp,1,1))
  toJoin.append(C.object)
bpy.ops.object.select_all(action="DESELECT")
for ob in toJoin: ob.select_set(True)
bpy.ops.object.join()
vframob = C.object

toJoin = []
for h in hpts:
  bpy.ops.mesh.primitive_cube_add(1,
    location=(0.5+ihwnd+thp-ptd,0,h-0.5*thp), scale=(1,1,thp))
  toJoin.append(C.object)
  bpy.ops.mesh.primitive_cube_add(1,
    location=(-0.5-ihwnd-thp+ptd,0,h-0.5*thp), scale=(1,1,thp))
  toJoin.append(C.object)
bpy.ops.object.select_all(action="DESELECT")
for ob in toJoin: ob.select_set(True)
bpy.ops.object.join()
hframob = C.object

# Notch the horizontal plies to make a reinforced central vertical column:
bpy.ops.mesh.primitive_cube_add(1, location=(0,0,0),
  scale=(2*(ihwnd+thp),0.3*boxheight,boxheight))
centralStrutSpace = C.object
DiffEq(hframob, centralStrutSpace)
D.objects.remove(centralStrutSpace, do_unlink=True)

# Notch them again around the prime meridian:
bpy.ops.mesh.primitive_cube_add(1, location=(0,0,0),
  scale=(2*(ihwnd+thp),boxheight*2/phi,boxheight))
notch = C.object
coredode = MakeMesh("DodeIn", (0.875*boxheight*V,E,F))
DiffEq(notch,coredode)
DiffEq(hframob, notch)
D.objects.remove(coredode, do_unlink=True)
D.objects.remove(notch, do_unlink=True)

# Insert horizontal frame parts into vertical frame parts:
DiffEq(vframob, hframob)

# Remove frame parts outside of bigdode:
framedode = MakeMesh("DodeFrame", ((boxheight-2*(thm-mtd))*V,E,F))
AndEq(hframob,framedode)
AndEq(vframob,framedode)
D.objects.remove(framedode, do_unlink=True)

lambdaJoint(doobj, hframob, 0.0125, mtd)
# lambdaJoint(doobj, vframob, 0.0125, mtd)
# DiffEq(doobj, hframob)
# DiffEq(doobj, vframob)

