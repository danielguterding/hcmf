CXX      = clang++
CXXFLAGS = -Wall -O3 -I/home/guterding/local/eigen-3.2.1 -fopenmp
LDFLAGS  = -lm -lboost_system

OBJECTS = main.o hamiltonian.o
DEFINES =

all : $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(DEFINES) $(OBJECTS) $(LDFLAGS) -o heisenberg
	
main.o : main.cpp hamiltonian.hpp typedefs.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c main.cpp -o main.o
	
hamiltonian.o : hamiltonian.cpp hamiltonian.hpp typedefs.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c hamiltonian.cpp -o hamiltonian.o
	
clean:
	rm heisenberg *.o
