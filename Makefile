
CC = g++ 


LIBS = -lm -lgsl -lblas
#LIBS = -L/usr/lib /usr/lib/libqthreads.so.0 -lguile -ldl -lreadline -ltermcap -lm

INCLUDES = -I ./include
#INCLUDES = -I/usr/include/g++-2 -I/usr/lib/gtkmm/include -I/usr/lib/sigc++/include -I/usr/lib/glib/include -I/usr/include/gtk-1.2 -I/usr/include/glib-1.2

#CFLAGS = -g -Wall -Wno-return-type $(INCLUDES) -DSWIG_GLOBAL
#CFLAGS = -g -Wall -Werror

CFLAGS = -g -Wall -O3 -g $(INCLUDES) 
#CFLAGS = -g -Wall -g $(INCLUDES) 

mazesim: maze.h neat.o network.o nnode.o link.o trait.o gene.o genome.o innovation.o organism.o species.o population.o experiments.o noveltyexp.o neatmain.o noveltyset.o histogram.o calc_evol.o #neatswig_wrap.o visual.o
	$(CC) $(CFLAGS) neat.o network.o nnode.o link.o trait.o gene.o genome.o innovation.o organism.o species.o population.o experiments.o neatmain.o noveltyexp.o noveltyset.o histogram.o calc_evol.o -o mazesim $(LIBS)
#	$(CC) $(CFLAGS) $(LIBS) networks.o genetics.o visual.o experiments.o neatswig_wrap.o neatmain.o -o neat `gtkmm-config --cflags --libs`

########################
calc_evol.o: calc_evol.cpp calc_evol.h
	$(CC) $(CFLAGS) -c calc_evol.cpp
histogram.o: histogram.cpp histogram.h
	$(CC) $(CFLAGS) -c histogram.cpp 
neat.o: neat.cpp neat.h
	$(CC) $(CFLAGS) -c neat.cpp -o neat.o

network.o: network.cpp network.h neat.h neat.o  
	$(CC) $(CFLAGS) -c network.cpp -o network.o

nnode.o: nnode.cpp nnode.h    
	$(CC) $(CFLAGS) -c nnode.cpp -o nnode.o

link.o: link.cpp link.h
	  $(CC) $(CFLAGS) -c link.cpp -o link.o

trait.o: trait.cpp trait.h
	  $(CC) $(CFLAGS) -c trait.cpp -o trait.o

gene.o: gene.cpp gene.h
	  $(CC) $(CFLAGS) -c gene.cpp -o gene.o

genome.o: genome.cpp genome.h
	  $(CC) $(CFLAGS) -c genome.cpp -o genome.o

innovation.o: innovation.cpp innovation.h
	  $(CC) $(CFLAGS) -c innovation.cpp -o innovation.o

organism.o: organism.cpp organism.h    
	$(CC) $(CFLAGS) -c organism.cpp -o organism.o

species.o: species.cpp species.h organism.h
	  $(CC) $(CFLAGS) -c species.cpp -o species.o

population.o: population.cpp population.h organism.h
	  $(CC) $(CFLAGS) -c population.cpp -o population.o

experiments.o: experiments.cpp experiments.h network.h species.h
	$(CC) $(CFLAGS) -c experiments.cpp -o experiments.o

neatmain.o: neatmain.cpp neatmain.h neat.h population.h
	$(CC) $(CFLAGS) -c neatmain.cpp -o neatmain.o

noveltyexp.o: noveltyexp.cpp noveltyset.h experiments.h maze.h 
	$(CC) $(CFLAGS) -c noveltyexp.cpp -o noveltyexp.o

noveltyset.o: noveltyset.cpp noveltyset.h
	$(CC) $(CFLAGS) -c noveltyset.cpp -o noveltyset.o



########################

clean:
	rm -f noveltyexp.o neat.o network.o nnode.o link.o trait.o gene.o genome.o innovation.o organism.o species.o population.o experiments.o neatmain.o rtneat
