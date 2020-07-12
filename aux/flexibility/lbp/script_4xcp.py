from pymol.cgo import *
from pymol import cmd

cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_transparency", 0.5)
cmd.set("sphere_transparency", 0.7)
cmd.set("ray_trace_mode",  0)
cmd.set("two_sided_lighting", "on")
cmd.set("reflect", 0.9)
cmd.set("ambient", 0.1)
cmd.set('''ray_opaque_background''', '''off''')

cmd.load("avg_4xcp.pdb")
cmd.load("cav_4xcp_1.pdb")
cmd.load("ch_4xcp.pdb")
cmd.color("deepsalmon", "avg_4xcp")
cmd.color("forest", "cav_4xcp_1")
cmd.color("blue", "ch_4xcp")

cmd.hide("spheres")
cmd.hide("sticks")
cmd.show("nonbonded")
cmd.show("lines", "ch_4xcp")

cmd.set_view (\
  '''0.483242571,    0.869933009,   -0.098409966,\
    -0.440265745,    0.338637143,    0.831549525,\
     0.756720006,   -0.358518451,    0.546650290,\
    -0.000332687,   -0.000202216, -108.338195801,\
    33.437988281,   35.142623901,   33.455139160,\
    51.865219116,  164.737243652,  -20.000000000''')

cmd.png("4xcp.png", width=500, height=400, dpi=200, ray=1)

cmd.quit()
