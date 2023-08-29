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

cmd.load("tctex.pdb")
cmd.color("darksalmon", "tctex")

cmd.load("ch_acb.pdb")
cmd.color("blue", "ch_acb")

cmd.load("ch_adb.pdb")
cmd.color("blue", "ch_adb")

cmd.load("acb_1.pdb")
cmd.color("forest", "acb_1")

cmd.load("adb_1.pdb")
cmd.color("forest", "adb_1")

cmd.hide("spheres")
cmd.hide("sticks")
cmd.show("nonbonded")
cmd.show("lines", "ch_acb")
cmd.show("lines", "ch_adb")

cmd.set_view (\
  '''0.559335649,    0.741895616,    0.369766444,\
     0.543353379,    0.008738876,   -0.839447320,\
    -0.626014233,    0.670448303,   -0.398220122,\
    -0.000068981,    0.000695011, -109.204299927,\
   -40.231609344,  -27.832706451,  -30.915164948,\
  -147.223098755,  365.669708252,  -20.000000000''')

cmd.png("tctex.png", width=500, height=400, dpi=200, ray=1)

cmd.quit()
