GUROBIROOT=/pathToGurobi	# path to the GUROBI root folder, example below:
# GUROBIROOT=/data/apps/gurobi1000
# if your GUROBI version is not 10.0, change 'gurobi1000' to the appropriate version

CFLAGS	= -m64 -O3 -std=c++11
INC		= ${GUROBIROOT}/linux64/include/
CPPLIB   = -L${GUROBIROOT}/linux64/lib/ ${GUROBIROOT}/linux64/lib/libgurobi_g++5.2.a -lgurobi100

all: CondDeMix

CondDeMix: objects
	g++ $(CFLAGS) -o $@ main.cpp DataMatrix.o IO.o method.o -I$(INC) $(CPPLIB) -lm

objects:
	g++ -c $(CFLAGS) DataMatrix.cpp -o DataMatrix.o
	g++ -c $(CFLAGS) IO.cpp -o IO.o
	g++ -c $(CFLAGS) method.cpp -o method.o -I$(INC) $(CPPLIB) -lm

clean:
	rm -f *.o CondDeMix
