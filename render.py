# SPDX-License-Identifier: MIT
# === CONTEXT ===
# Purpose:     Headless preview renders of the generated STL (no GUI available here).
# Decisions:   Cycles on CPU. Eevee/Workbench need a GL context, which `blender -b`
#              on this host does not provide.
# Constraints: Blender 5.x API (bpy.ops.wm.stl_import); tested with 5.2.2. Both paths must
#              be absolute: with a relative output prefix, `-b` fails to save
#              ("cannot save: No such file or directory").
# Related:     Renders the output of orchid_pot.py.
# Usage:       blender -b --factory-startup --python render.py -- /abs/orchid_pot.stl /abs/preview
# === END CONTEXT ===
import math
import sys

import bpy
from mathutils import Vector

stl, prefix = sys.argv[sys.argv.index("--") + 1:][:2]

bpy.ops.wm.read_factory_settings(use_empty=True)
bpy.ops.wm.stl_import(filepath=stl)
obj = bpy.context.selected_objects[0]
bpy.ops.object.shade_smooth()

mat = bpy.data.materials.new("petg")
mat.use_nodes = True
bsdf = mat.node_tree.nodes["Principled BSDF"]
bsdf.inputs["Base Color"].default_value = (0.10, 0.45, 0.25, 1)
bsdf.inputs["Roughness"].default_value = 0.35
obj.data.materials.append(mat)

scn = bpy.context.scene
scn.render.engine = "CYCLES"
scn.cycles.device = "CPU"
scn.cycles.samples = 64
scn.cycles.use_denoising = True
scn.render.resolution_x = scn.render.resolution_y = 1000
scn.world = bpy.data.worlds.new("w")
scn.world.use_nodes = True
scn.world.node_tree.nodes["Background"].inputs["Color"].default_value = (0.9, 0.9, 0.88, 1)
scn.world.node_tree.nodes["Background"].inputs["Strength"].default_value = 0.6

sun = bpy.data.objects.new("sun", bpy.data.lights.new("sun", "SUN"))
sun.data.energy = 4.0
sun.rotation_euler = (math.radians(40), 0, math.radians(30))
scn.collection.objects.link(sun)

cam = bpy.data.objects.new("cam", bpy.data.cameras.new("cam"))
cam.data.clip_end = 10000
scn.collection.objects.link(cam)
scn.camera = cam
target = Vector((0, 0, 50))

# lens > 0: perspective focal length [mm]; lens < 0: orthographic, -lens = ortho scale
views = {
    "persp": ((330, -330, 260), 70, target),
    "side": ((0, -600, 55), -175, Vector((0, 0, 55))),
    "top": ((0, -1, 600), -175, Vector((0, 0, 0))),
    "detail": ((150, -150, 110), 35, Vector((40, -40, 60))),
}
for name, (loc, lens, tgt) in views.items():
    cam.location = loc
    cam.data.type = "PERSP" if lens > 0 else "ORTHO"
    cam.data.lens = max(lens, 1)
    cam.data.ortho_scale = -lens if lens < 0 else 1
    cam.rotation_euler = (tgt - cam.location).to_track_quat("-Z", "Y").to_euler()
    scn.render.filepath = f"{prefix}_{name}.png"
    bpy.ops.render.render(write_still=True)
