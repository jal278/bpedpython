import math

def scale(x,mn,mx):
 return (x-mn)/(mx-mn)

def mass_excess(x):
 return scale(x,-91656.1,119600.0)

def mse(d1,d2):
 terr=0.0
 for k in range(len(d1)):
  err=d1[k]-d2[k]
  err*=err
  terr+=err
 return terr/len(d1)

import random

"""
train=[]
test=[]
for x in range(50):
 train.append(random.uniform(-50000.0,50000.0))
 test.append(train[-1]+random.uniform(-5000,5000))

train2=[mass_excess(k) for k in train]
test2=[mass_excess(k) for k in test]
print math.sqrt(mse(train,test))/1000.0
x=mse(train2,test2)
"""

x=float(raw_input("test perf?"))

def convert_fitness(x):
 low= -91656.1 
 high =119600.0
 rng = (high-low)
 #mse = 1.0-x
 #print mse
 mse=1.0-x
 mse *= rng*rng
 return math.sqrt(mse)/1000.0
