from pymol import cmd,stored
stored.list=[]
cmd.iterate("(sele)","stored.list.append(index)")
print stored.list
