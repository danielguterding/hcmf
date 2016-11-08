CXX      = clang++
CXXFLAGS = -Wall -O3 -std=c++11 -fopenmp -I/home/guterding/local/eigen3
LDFLAGS  = -lm -lboost_system -lboost_filesystem

OBJECTS = main.o hamiltonian.o files.o
DEFINES =

all : $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(DEFINES) $(OBJECTS) $(LDFLAGS) -o heisenberg
	
main.o : main.cpp hamiltonian.hpp files.hpp typedefs.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c main.cpp -o main.o
	
hamiltonian.o : hamiltonian.cpp hamiltonian.hpp typedefs.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c hamiltonian.cpp -o hamiltonian.o
	
files.o : files.cpp files.hpp typedefs.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c files.cpp -o files.o	
	
clean:
	rm heisenberg *.o
