import glob
ode_path = '/home/jlehman/ode-0.11/'

#env = Environment(CCFLAGS = ' -DdTRIMESH_ENABLED -DdDOUBLE -DGRAPHICS -g -I./include')
env = Environment(CCFLAGS = ' -DdTRIMESH_ENABLED -DdDOUBLE -DGRAPHICS -O3 -g -I./include')
env.AppendENVPath('CPLUS_INCLUDE_PATH', ode_path+'include')

#current=['biped.cpp',"biped_stuff.cpp","ConfigFile.cpp"] 
#rtneat=glob.glob('rtneat/*.cpp')

current=glob.glob('*.cpp')

allsrc=current #+rtneat
allsrc.remove("mazeApp.cpp")
allsrc.remove("mazeDlg.cpp")
env.Program('mazesim', allsrc,LIBS=['m','ode'],LIBPATH=['.','/usr/lib/','/usr/local/lib'])
