#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <exception>

#include "lib/macrologger.h"
#include "poisson.hpp"

double F(double x, double y)  { return 23.2; }
double Phi(double x, double y)  { return 42.3; }

int main(int argc, char** argv) {
	MPI_Init(&argc,&argv);

	if (argc != 5) {
		printf("Usage: './poisson x0 y0 l m' where (x0, y0) is the left"
               "bottom corner of the square, l is the square size and"
               "m*m is the number of dots, so m is grid size.");
		return -1;
	}

	double x0 = atof(argv[1]);
	double y0 = atof(argv[2]);
    double square_size = atof(argv[3]);
	int grid_size = atoi(argv[4]);

	try {
		Poisson poiss(x0, y0, square_size, grid_size, &F, &Phi);
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
