from bipedpy import bipedsim 

#simple fitness function -- how close do we get to the goal at the end of
#a simulation?
def fitness(robot):
 return abs(robot.get_x())+abs(robot.get_y())

#initialize maze stuff with "medium maze" 
bipedsim.initmaze("biped.ne")
bipedsim.random_seed()

#create initial genome
robot=bipedsim()

#initalize it randomly and mutate it once for good measure
robot.init_rand()
robot.mutate()

#run genome in the maze simulator
robot.map()

#current fitness set to how fit our first dude is
cur_fitness=fitness(robot)

iterations=0

#loop until we discover a solution to the maze
#this is using simple hillclimbing algorithm (not neat)
while not robot.solution():

 #create new mutant from current best
 newrobot=robot.copy()
 newrobot.mutate()
 newrobot.map()
 new_fitness=fitness(newrobot)

 #new best?
 if(new_fitness>cur_fitness):
  print "New best fitness: " , new_fitness
  cur_fitness=new_fitness
  robot=newrobot

 iterations+=1

 if(iterations%1000 == 0):
  robot.save("solution_brain")
  print "iterations:%d" % iterations

print "Maze solved"
print robot.get_x(),robot.get_y(),robot.solution()

#save brain (you can visualize it with ./maze medium_maze.txt solution_brain
robot.save("solution_brain")
