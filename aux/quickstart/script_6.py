from pymol.cgo import *
from pymol import cmd

cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_transparency", 0.7)
cmd.set("ray_trace_mode",  0)
cmd.set("two_sided_lighting", "on")
cmd.set("reflect", 0.9)
cmd.set("ambient", 0.1)
cmd.set('''ray_opaque_background''', '''off''')

cmd.load("input_pdb.pdb")
cmd.load("4_cavity_1.pdb")
cmd.set_view (\
  '''0.237595469,   -0.169580087,    0.956439614,\
     0.387604266,   -0.886297047,   -0.253433615,\
     0.890669584,    0.430938631,   -0.144854113,\
    -0.000374220,    0.000006055, -100.353157043,\
    31.699275970,   32.316436768,   33.389636993,\
    61.715843201,  138.970840454,  -20.000000000''')
cmd.png("quickstart_6.png", width=400, height=300, dpi=100, ray=1)
