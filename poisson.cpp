#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "lib/macrologger.h"
#include "poisson.hpp"


Poisson::Poisson(double x0, double y0, double square_size, int grid_size)
	: dots_per_proc(NULL) {
	/* It is pretty important to initialize dots_per_proc and dots to NULL:
     * otherwise, unallocated, they will be freed in the destructor
     */
	for (int i = 0; i < 2; i++ ) {
		dots[i] = NULL;
	}
	/* Set rank and size */
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	LOG_DEBUG_MASTER("Size is %d\n", size);

	/* define processors grid and block unneeded processors */
	if (size < 4) {
		throw PoissonException("Error: we require at least 4 processors");
	}
	proc_grid_size = static_cast<int>(floor(sqrt(size)));
	LOG_DEBUG_MASTER("Proc grid size is %d\n", proc_grid_size);

	/* we use only square grid of processors */
	if (rank >= proc_grid_size*proc_grid_size) {
		LOG_DEBUG("Process with rank %d is not used and will be blocked\n",
				  rank);
		/* Do I need this? */
		MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, 0, &comm);
		MPI_Barrier(MPI_COMM_WORLD);
		return;
	}
	MPI_Comm_split(MPI_COMM_WORLD, 1948, 0, &comm);
	ranks[0] = rank % proc_grid_size;
	ranks[1] = rank / proc_grid_size;

	if (grid_size < proc_grid_size) {
		throw PoissonException("Error: grid size must be >= proc grid size\n");
	}
	DistributeDots(grid_size);
	CalculateDots(x0, y0, square_size, grid_size);


	LOG_DEBUG_MASTER("Square starts at (%f, %f) of size %f, grid size %d\n",
		   x0, y0, square_size, grid_size);

	/* To synchronize with unused processes */
	MPI_Barrier(MPI_COMM_WORLD);
}

/* Find this processor's dots, namely, it sets dots_num[] and borders[].
 * Since everything is square, the distribution is symmetrical.
 * Returns -1, if this processor is not used.
 */
void Poisson::DistributeDots(int grid_size) {
	dots_per_proc = new int[proc_grid_size]();
	int dots_left_to_distribute = grid_size;
	int i = 0;
	while (dots_left_to_distribute > 0) {
		dots_per_proc[i]++;
		i++;
		i = i % proc_grid_size;
		dots_left_to_distribute--;
	}

	/* Learn how many dots will this processor have */
	for (int i = 0; i < 2; i++) {
		dots_num[i] = dots_per_proc[ranks[i]];
	}

	/* Does this processor touches the borders? */
	int borders_rank_check[4] =
		{0, 0, proc_grid_size - 1, proc_grid_size - 1};
	for (int i = 0; i < 4; i++) {
		borders[i] = (ranks[i % 2] == borders_rank_check[i]);
	}

	LOG_DEBUG("My rank is %d, i.e. (%d , %d), I own %d x %d dots and right touching is %d\n",
			  rank, ranks[0], ranks[1], dots_num[0], dots_num[1], borders[2]);
}

/* Calculate dots area this processor works on
 * dots_per_proc and dots_num must be set at this point
 */
void Poisson::CalculateDots(double x0, double y0, double square_size,
						   int grid_size) {
	printf("0000000000000000000000\n");
	/* Grid step. Grid is uniform and symmetrical. */
	double step = square_size / grid_size;
	LOG_DEBUG_MASTER("Step is %f", step);
	printf("AAAAAAAAAAAAAAAAAAAAAAAAA\n");
	/* bottom left point of this processor working area */
	double bottom_left[2] = {x0, y0};
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < ranks[i]; j++) {
			bottom_left[i] += step*dots_per_proc[j];
		}
	}

	printf("BBBBBBBBBBBBBBBBBBBBBBBBBBBb\n");
	/* now calculate the dots themselves */
	for (int i = 0; i < 2; i++) {
		dots[i] = new double[dots_num[i]];
		if (rank == 2) {
			printf("Allocated %p\n", (void*)dots[i]);
		}
		dots[i][0] = bottom_left[i];
		for (int j = 1; j < dots_num[i]; j++) {
			dots[i][j] += step;
		}
	}
	printf("FFFFFFFFFFFFFFFFFFFFFFFF\n");

	if (rank == 2) {
		printf("X values of proc with rank 2:\n");
		for (int i = 0; i < dots_num[0]; i++)
			printf("%f ", dots[0][i]);
		printf("\n");
		printf("Y values of proc with rank 2:\n");
		for (int i = 0; i < dots_num[1]; i++)
			printf("%f ", dots[1][i]);
	}
}


Poisson::~Poisson() {
	if (dots_per_proc)
		delete[] dots_per_proc;
    for (int i = 0; i < 2; i++) {
		if (dots[i]) {
			if (rank == 2) {
				printf("Not Freeing %p\n", (void*)dots[i]);
			}
			delete[] dots[i];
		}
	}
}
