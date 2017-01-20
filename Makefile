CPPFLAGS := $(CPPFLAGS)

ifeq ($(MODE),release)
   CPPFLAGS := $(CPPFLAGS) -D"LOG_LEVEL=2" -O3
else
   mode = debug
   CPPFLAGS :=  $(CPPFLAGS) -D"LOG_LEVEL=3" -g -O0
endif

ifeq ($(COMPILER),ibm)
   CXX := mpixlcxx_r
else
   # gcc 4.1.2 on BlueGene
   CXX := mpicxx
   CPPFLAGS := $(CPPFLAGS) -std=gnu++98
endif


all: main.cpp poisson.o
	$(CXX) $(CPPFLAGS) main.cpp poisson.o matrix.o -o poisson

poisson.o: matrix.o

.PHONE: clean
clean:
	- rm matrix.o poisson.o poisson
