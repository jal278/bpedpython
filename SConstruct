import glob
ode_path = '/home/joel/ode-0.11.1/'

#env = Environment(CCFLAGS = ' -DdTRIMESH_ENABLED -DdDOUBLE -DGRAPHICS -g -I./include')
env = Environment(CCFLAGS = ' -march=native -DdTRIMESH_ENABLED -DdDOUBLE -DGRAPHICS -O2 -g -I./include')
env.AppendENVPath('CPLUS_INCLUDE_PATH', ode_path+'include')

#current=['biped.cpp',"biped_stuff.cpp","ConfigFile.cpp"] 
#rtneat=glob.glob('rtneat/*.cpp')

current=glob.glob('*.cpp')

allsrc=current #+rtneat
allsrc.remove("mazeApp.cpp")
allsrc.remove("mazeDlg.cpp")
env.Program('mazesim', allsrc,LIBS=['pthread','tcmalloc','m','ode'],LIBPATH=['.','/usr/lib/','/usr/local/lib'])
