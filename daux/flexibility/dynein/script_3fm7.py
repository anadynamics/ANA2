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

cmd.load("h3fm7.pdb")
cmd.color("darksalmon", "chain A")
cmd.color("darksalmon", "chain B")
cmd.color("marine", "chain C")
cmd.color("marine", "chain D")
cmd.color("wheat", "chain E")
cmd.color("wheat", "chain F")

cmd.set_view (\
  '''-0.111517593,   -0.977693856,   -0.177933186,\
    -0.873985291,    0.011277429,    0.485806793,\
    -0.472963661,    0.209689617,   -0.855753779,\
    -0.001382068,   -0.000694713, -115.162315369,\
   -43.635929108,  -42.908710480,  -44.754673004,\
  -148.208801270,  378.666076660,  -20.000000000''')

cmd.png("3fm7.png", width=600, height=300, dpi=200, ray=1)

cmd.quit()
