import sys
import random


test=False

calc_evo=True
extinction=True
seed=-1
outfile="out"
nefile="biped.ne"
interval=50000

if(len(sys.argv)>1):
 seed = int(sys.argv[1])
 fname = sys.argv[2]
 nefile=sys.argv[3]

disp=False
SZX=SZY=400
screen = None
test=True

if test:
 calc_evo=False
 disp=True

if disp:
 import pygame
 from pygame.locals import *
 pygame.init()
 pygame.display.set_caption('Viz')
 screen =pygame.display.set_mode((SZX,SZY))
 
 background = pygame.Surface(screen.get_size())
 background = background.convert()
 background.fill((250, 250, 250))

def render(pop):
 global screen,background
 screen.blit(background, (0, 0))
 for robot in pop:
  x=clamp(robot.get_x(),16.0)*SZX
  y=clamp(robot.get_y(),16.0)*SZY 
  #y=mazepy.feature_detector.endy(robot)*SZY
  rect=(int(x),int(y),5,5)
  pygame.draw.rect(screen,(255,0,0),rect,0)
 pygame.display.flip()

from entropy import *
evo_fnc = calc_evolvability_cnt


if(__name__=='__main__'):
 #evo_fnc = calc_evolvability_entropy
 #initialize maze stuff with "medium maze" 

 #mazepy.mazenav.initmaze("hard_maze_list.txt")
 bipedpy.bipedsim.initmaze(nefile) #"biped.ne")
 #mazepy.mazenav.initmaze("medium_maze_list.txt")
 if(seed==-1):
  bipedpy.bipedsim.random_seed()
 else:
  random.seed(seed)
  bipedpy.bipedsim.seed(seed)

 eflag=False
 robot=None

 whole_population=[]
 psize=2000
 repop=0

 robot=bipedpy.bipedsim()
 robot.load_new(fname)

 for k in range(psize):
  robot=bipedpy.bipedsim()
  robot.init_rand()
  robot.mutate()
  robot.map()
  whole_population.append(robot)
 solved=False

 render(whole_population)
 evo=evo_fnc(robot,1000)
 print evo

 """
 robot=mazepy.mazenav()
 robot.init_rand()
 robot.mutate()
 robot.map()
 print "evolvability:", evo_fnc(robot,1000)
 
 optimize_evolvability(child)
 """
