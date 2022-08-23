CC	      = g++
CC_FLAGS        = -g3 -O3 -Wall -D_GLIBCXX_DEBUG -I  /opt/homebrew/Cellar/gsl/2.7.1/include/
LD_FLAGS        = -L/opt/homebrew/Cellar/gsl/2.7.1/lib  -lgsl -lgslcblas -lm -lstdc++ 
SI_OBJECTS      = simulation_run_infinite.o utilities_comp.o

sim_i:  $(SI_OBJECTS)
	$(CC) $(CC_FLAGS) $(SI_OBJECTS) -o comp_sim_infinite $(LD_FLAGS)
simulation_run_infinite.o: simulation_run_infinite.cpp
	$(CC) $(CC_FLAGS) -c simulation_run_infinite.cpp
utilities_comp.o: utilities_comp.cpp
	$(CC) $(CC_FLAGS) -c utilities_comp.cpp


