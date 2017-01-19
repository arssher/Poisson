#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>

#include "lib/macrologger.h"
#include "poisson.hpp"


Poisson::Poisson(double x0, double y0, double square_size, int grid_size,
				 double (*_F)(double x, double y),
				 double (*_Phi)(double x, double y))
	: F(_F), Phi(_Phi), dots_per_proc(NULL) {
	/* It is pretty important to initialize ptrs to NULL;
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
	LOG_DEBUG_MASTER("Square starts at (%f, %f) of size %f, grid size %d\n",
					 x0, y0, square_size, grid_size);
	DistributeDots(grid_size);
	CalculateDots(x0, y0, square_size, grid_size);
	Solve();



	/* To synchronize with unused processes */
	MPI_Barrier(MPI_COMM_WORLD);
}

/* Find this processor's dots, namely, it sets dots_num[2], borders[4] and
 * inner_dots_range[4].
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

	/* set inner_dots_range */
	int dots_range[4] = {0, 0, dots_num[0], dots_num[1]};
	int add_or_substract[4] = {1, 1, -1, -1};
	for (int i = 0; i < 4; i++) {
	    inner_dots_range[i] = dots_range[i] + add_or_substract[i] * borders[i];
	}

	LOG_DEBUG("My rank is %d, i.e. (%d , %d), I own %d x %d dots and right touching is %d\n",
			  rank, ranks[0], ranks[1], dots_num[0], dots_num[1], borders[2]);
}

/* Calculate values of dots this processor works on
 * dots_per_proc and dots_num must be set at this point
 * sets dots[2][]
 */
void Poisson::CalculateDots(double x0, double y0, double square_size,
						   int grid_size) {
	/* Grid step. Grid is uniform and symmetrical. */
	double step = square_size / grid_size;
	LOG_DEBUG_MASTER("Step is %f", step);
	/* bottom left point of this processor working area */
	double bottom_left[2] = {x0, y0};
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < ranks[i]; j++) {
			bottom_left[i] += step*dots_per_proc[j];
		}
	}

	/* now calculate the dots themselves */
	for (int i = 0; i < 2; i++) {
		dots[i] = new double[dots_num[i]];
		dots[i][0] = bottom_left[i];
		for (int j = 1; j < dots_num[i]; j++) {
			dots[i][j] += step;
		}
	}

	/* For debugging */
	// if (rank == 0) {
	// 	printf("X values of proc with rank 0:\n");
	// 	for (int i = 0; i < dots_num[0]; i++)
	// 		printf("%f ", dots[0][i]);
	// 	printf("\n");
	// 	printf("Y values of proc with rank 0:\n");
	// 	for (int i = 0; i < dots_num[1]; i++)
	// 		printf("%f ", dots[1][i]);
	// 	printf("\n");
	// }
}

/* Now, with dots distributed and calculated, run the solver */
void Poisson::Solve() {
	resid_matr = Matrix(dots_num[0], dots_num[1]);
	sol_matr = Matrix(dots_num[0], dots_num[1]);
	InitSolMatr();
}

/* We initialize sol_matr to Phi on the borders, and to random values in the
 * inner dots.
 */
void Poisson::InitSolMatr() {
	/* The borders. Each corner will be initialized twice, but that's not a
	 * big deal and simplifies code a bit.
	 */

	int fixed_coord[4] = {0, 0, dots_num[0], dots_num[1]};
	for (int r = 0; r < 4; r++) {
		if (borders[r]) {
			if (r % 2  == 0) { /* vertical border */
				for (int j = 0; j < dots_num[1]; j++) {
					sol_matr(fixed_coord[r], j) =
						(*Phi)(dots[0][fixed_coord[r]], dots[1][j]);
				}
			}
			else { /* horizontal border */
				for (int i = 0; i < dots_num[0]; i++) {
					sol_matr(i, fixed_coord[r]) =
						(*Phi)(dots[0][i], dots[1][fixed_coord[r]]);
				}
			}
		}
	}

	/* Inner part */
	srand(time(NULL));
	for (int i = inner_dots_range[0]; i < inner_dots_range[2]; i++)
		for (int j = inner_dots_range[1]; j < inner_dots_range[3]; j++) {
			sol_matr(i, j) = rand();
		}

	/* For debugging */
	if (rank == 8) {
		printf("sol values of proc with rank 8:\n");
		sol_matr.Print();
	}
}


Poisson::~Poisson() {
	if (dots_per_proc)
		delete[] dots_per_proc;
    for (int i = 0; i < 2; i++) {
		if (dots[i]) {
			delete[] dots[i];
		}
	}
}
