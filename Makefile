# if MODE set to 'release', enable optimizations and hide debug info
# if COMPILER set to 'ibm', IBM compiler is used, gcc otherwise
# if OMP is set to something, we compile with OpenMP support using IBM
# compiler, COMPILER var is ignored

ifeq ($(MODE),release)
	CPPFLAGS := $(CPPFLAGS) -D"LOG_LEVEL=2" -O3
else
	CPPFLAGS :=  $(CPPFLAGS) -D"LOG_LEVEL=3" -g -O0
endif

ifdef OMP
CXX := mpixlcxx_r
CPPFLAGS := $(CPPFLAGS) -qsmp=omp
else
ifeq ($(COMPILER),ibm)
CXX := mpixlcxx
else
# gcc 4.1.2 on BlueGene
CXX := mpicxx
CPPFLAGS := $(CPPFLAGS) -std=gnu++98
endif
endif


all: main.cpp poisson.o
	$(CXX) $(CPPFLAGS) main.cpp poisson.o matrix.o -o poisson

poisson.o: matrix.o

.PHONE: clean
clean:
	- rm matrix.o poisson.o poisson
	- rm -rf res/
