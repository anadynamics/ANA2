from pymol.cgo import *
from pymol import cmd

cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_transparency", 0.7)
cmd.set("ray_trace_mode",  0)
cmd.set("two_sided_lighting", "on")
cmd.set("reflect", 0.9)
cmd.set("ambient", 0.1)
cmd.set('''ray_opaque_background''', '''off''')

cmd.load("4xcp_hv.pdb")
cmd.color("red", "4xcp_hv")
cmd.load("hv_cavity_1.pdb")
cmd.color("salmon", "hv_cavity_1")

cmd.set_view (\
  '''0.456030250,    0.888168931,    0.056359597,\
    -0.276825786,    0.081376061,    0.957462490,\
     0.845804632,   -0.452238709,    0.282975733,\
    -0.000969257,   -0.000390973, -118.043930054,\
    33.933498383,   32.315837860,   33.869487762,\
    83.247825623,  152.777328491,  -20.000000000''')

cmd.png("dynamics_hv.png", width=500, height=400, dpi=100, ray=1)

cmd.quit()
