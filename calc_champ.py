import os
import sys
import glob
#template="./mazesim --eval --seed ./mazeres/nosticky/%s%d_rtgen_weakfirst.bst > out.txt"

#template="./mazesim -m medium_maze_list.txt --eval --seed %s > out.txt"
template="./mazesim -s maze_mut.ne -m hard_maze_list.txt --eval --seed %s > out.txt"

def read_in(fn):
 global template
 cmd = template % (fn)
 print cmd
 os.system(cmd)
 a=map(float,open("out.txt").read().split("\n")[-2].split())
 b=open("out.txt").read().split("\n")[-4]
 return a

ev=[]
for k in glob.glob("res/hard_champ/fr*novmuthard*.bst"):
 a=read_in(k)
 ev.append(a[0])

print ev
print len(ev)
print sum(ev)/len(ev)
