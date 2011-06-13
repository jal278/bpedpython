import os
import sys
template="./mazesim -s %s --eval --seed ./res/biped_champ/%s%d_fittest_final.bst --biped > out.txt"

def read_in(prefix,num,sg):
 global template
 cmd = template % (sg,prefix,num)
 #print cmd
 os.system(cmd)
 a=map(float,open("out.txt").read().split("\n")[-2].split())
 print a
 return a

fitss=[]
novss=[]
sg="biped.ne"

for k in range(7,8):
 fitss.append(read_in("fr_novbiped",k,sg))

#for k in range(1):
# novss.append(read_in("fr_fitbiped",k,sg))

f1,f2=zip(*fitss)
n1,n2=zip(*novss)

print "nov"
print sum(n1)/len(n1)
print sum(n2)/len(n2)
print "fit"
print sum(f1)/len(f1)
print sum(f2)/len(f2)
print n1
print f1
