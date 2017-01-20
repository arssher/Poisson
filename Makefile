CXX = mpicxx
CPPFLAGS := $(CPPFLAGS) -D"LOG_LEVEL=$(LOG_LEVEL)" -g -O0 -std=gnu++98
# CXX = mpixlcxx_r
# CPPFLAGS := $(CPPFLAGS) -D"LOG_LEVEL=$(LOG_LEVEL)" -g -O0

all: main.cpp poisson.o
	$(CXX) $(CPPFLAGS) main.cpp poisson.o matrix.o -o poisson

poisson.o: matrix.o

.PHONE: clean
clean:
	- rm matrix.o poisson.o poisson
