import os
import sys
template="./mazesim --eval --seed ./bipedres/%s%d_fittest_final.bst > out.txt"

def read_in(prefix,num):
 global template
 cmd = template % (prefix,num)
 #print cmd
 os.system(cmd)
 a=map(float,open("out.txt").read().split("\n")[-2].split())
 print a
 return a

fitss=[]
novss=[]
for k in range(8):
 fitss.append(read_in("fitbiped",k))

for k in range(8):
 novss.append(read_in("novbiped",k))

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
