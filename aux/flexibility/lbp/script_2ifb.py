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

cmd.load("avg_2ifb.pdb")
cmd.load("cav_2ifb_1.pdb")
cmd.load("ch_2ifb.pdb")
cmd.color("deepsalmon", "avg_2ifb")
cmd.color("forest", "cav_2ifb_1")
cmd.color("blue", "ch_2ifb")

cmd.hide("spheres")
cmd.hide("sticks")
cmd.show("nonbonded")
cmd.show("lines", "ch_2ifb")

cmd.set_view (\
  '''0.050421745,   -0.859337866,    0.508910358,\
    -0.963399172,    0.092455946,    0.251580775,\
    -0.263245493,   -0.502969980,   -0.823232293,\
     0.000638990,    0.000308108,  -99.765609741,\
    27.708494186,   29.903253555,   29.054700851,\
    70.827804565,  128.748825073,  -20.000000000''')

cmd.png("2ifb.png", width=500, height=400, dpi=200, ray=1)

cmd.quit()
