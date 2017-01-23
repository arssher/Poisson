#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <exception>

#include "lib/macrologger.h"
#include "poisson.hpp"

double FTest(double x, double y)  { return 20; }
double PhiTest(double x, double y)  { return 42.3; }

double F2(double x, double y)  { return (x*x + y*y) * sin(x * y); }
double Phi2(double x, double y)  { return 1 + sin(x * y); }

double FN(double x, double y)  { return (x*x + y*y) / ((1 + x*y) * (1 + x*y));}
double PhiN(double x, double y)  { return log(1 + x*y); }

int main(int argc, char** argv) {
	MPI_Init(&argc,&argv);

	if (argc < 6) {
		printf("Usage: './poisson x0 y0 l m sdi [dump]' where\n"
               "(x0, y0) is the left bottom corner of the square\n"
			   "l is the square size\n"
               "m*m is the number of dots, so m is grid size\n"
			   "sdi is number of steep descent iterations\n"
			   "dump the solution to [dump] directory. Dump is not performed"
			   "if the argument if absent");
		return -1;
	}

	double x0 = atof(argv[1]);
	double y0 = atof(argv[2]);
    double square_size = atof(argv[3]);
	int grid_size = atoi(argv[4]);
	int sdi_iterations = atoi(argv[5]);
	char *dump_dir = NULL;
	if (argc >= 7)
		dump_dir = argv[6];

	try {
		Poisson poiss(x0, y0, square_size, grid_size, sdi_iterations,
					  &F2, &Phi2, dump_dir);
	}
	catch (std::exception &e) {
		fprintf(stderr, "%s\n", e.what());
	}


	MPI_Finalize();
	return 0;
}
