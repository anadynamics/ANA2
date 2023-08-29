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

cmd.load("avg_4uet.pdb")
cmd.load("cav_4uet_1.pdb")
cmd.load("ch_4uet.pdb")
cmd.color("wheat", "avg_4uet")
cmd.color("forest", "cav_4uet_1")
cmd.color("blue", "ch_4uet")

cmd.hide("spheres")
cmd.hide("sticks")
cmd.show("nonbonded")
cmd.show("lines", "ch_4uet")

cmd.set_view (\
  '''0.669605851,    0.610031366,   -0.423654735,\
    -0.187055930,    0.690543890,    0.698672831,\
     0.718765140,   -0.388591349,    0.576507986,\
    -0.000425146,   -0.000188941, -109.491302490,\
    32.991695404,   34.373760223,   34.586250305,\
    80.516609192,  138.437667847,  -20.000000000''')

cmd.png("4uet.png", width=500, height=400, dpi=200, ray=1)

cmd.quit()
