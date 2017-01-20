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

int main(int argc, char** argv) {
	MPI_Init(&argc,&argv);

	if (argc != 6) {
		printf("Usage: './poisson x0 y0 l m sdi' where\n"
               "(x0, y0) is the left bottom corner of the square\n"
			   "l is the square size\n"
               "m*m is the number of dots, so m is grid size\n"
			   "sdi is number of steep descent iterations");
		return -1;
	}

	double x0 = atof(argv[1]);
	double y0 = atof(argv[2]);
    double square_size = atof(argv[3]);
	int grid_size = atoi(argv[4]);
	int sdi_iterations = atoi(argv[5]);

	try {
		Poisson poiss(x0, y0, square_size, grid_size, sdi_iterations,
					  &F2, &Phi2);
	}
	catch (std::exception &e) {
		fprintf(stderr, "%s\n", e.what());
	}


	MPI_Finalize();
	return 0;

    // Initialize the MPI environment
    // MPI_Init(NULL, NULL);

    // // Get the number of processes
    // int world_size;
    // MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // // Get the rank of the process
    // int world_rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // // Get the name of the processor
    // char processor_name[MPI_MAX_PROCESSOR_NAME];
    // int name_len;
    // MPI_Get_processor_name(processor_name, &name_len);

    // // Print off a hello world message
    // printf("Hello world from processor %s, rank %d"
    //        " out of %d processors\n",
    //        processor_name, world_rank, world_size);

    // // Finalize the MPI environment.
    // MPI_Finalize();
}
