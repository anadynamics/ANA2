from pymol import cmd,stored
stored.list = []
cmd.iterate("(sele and name ca)", "stored.list.append(resi)")
print stored.list
