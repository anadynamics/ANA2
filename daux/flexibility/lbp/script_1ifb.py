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

cmd.load("avg_1ifb.pdb")
cmd.load("cav_1ifb_1.pdb")
cmd.load("ch_1ifb.pdb")
cmd.color("wheat", "avg_1ifb")
cmd.color("forest", "cav_1ifb_1")
cmd.color("blue", "ch_1ifb")

cmd.hide("spheres")
cmd.hide("sticks")
cmd.show("nonbonded")
cmd.show("lines", "ch_1ifb")

cmd.set_view (\
  '''0.783268511,    0.174335986,   -0.596722603,\
     0.145346895,   -0.984611809,   -0.096872225,\
    -0.604435563,   -0.010855785,   -0.796566486,\
    -0.000247055,    0.000777386,  -98.827529907,\
    34.539207458,   29.541465759,   35.226402283,\
    -1.131647825,  198.827590942,  -20.000000000''')

cmd.png("1ifb.png", width=500, height=400, dpi=200, ray=1)

cmd.quit()
