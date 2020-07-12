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

cmd.load("lc8.pdb")
cmd.color("sand", "lc8")

cmd.load("ch_ecf.pdb")
cmd.color("blue", "ch_ecf")

cmd.load("ch_edf.pdb")
cmd.color("blue", "ch_edf")

cmd.load("ecf_1.pdb")
cmd.color("forest", "ecf_1")

cmd.load("edf_1.pdb")
cmd.color("forest", "edf_1")

cmd.hide("spheres")
cmd.hide("sticks")
cmd.show("nonbonded")
cmd.show("lines", "ch_ecf")
cmd.show("lines", "ch_edf")

cmd.set_view (\
  '''0.030933606,    0.584483027,    0.810812235,\
    -0.455680430,   -0.713744700,    0.531896412,\
     0.889600039,   -0.385925800,    0.244258434,\
     0.000140842,   -0.000561976, -111.916397095,\
   -57.852790833,  -77.408004761,  -41.355598450,\
  -151.198226929,  375.676757812,  -20.000000000''')

cmd.png("lc8.png", width=500, height=400, dpi=200, ray=1)

cmd.quit()
