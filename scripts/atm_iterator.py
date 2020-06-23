from pymol import cmd,stored
stored.list=[]
cmd.iterate("(wachin)","stored.list.append(index)")
print stored.list
