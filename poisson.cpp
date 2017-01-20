#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

#include "lib/macrologger.h"
#include "poisson.hpp"

static double zero_filler(double x, double y) { return 0.0; };

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
	for (int i = 0; i < 4; i++ ) {
		send_buffers[i] = NULL;
		recv_buffers[i] = NULL;
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

	/* We demand it to avoid fuss when the same dot occurs to be several
	 * corners, when dot matrix is less than 2x2.
     */
	if (grid_size < 2*proc_grid_size) {
		throw PoissonException("Error: grid size must be >= 2*proc grid size\n");
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

	LOG_DEBUG("My rank is %d, i.e. (%d , %d), I own %d x %d dots and top touching is %d\n",
			  rank, ranks[0], ranks[1], dots_num[0], dots_num[1], borders[3]);
}

/* Calculate values of dots this processor works on
 * dots_per_proc and dots_num must be set at this point
 * sets dots[2][]
 */
void Poisson::CalculateDots(double x0, double y0, double square_size,
						   int grid_size) {
	/* Grid step. Grid is uniform and symmetrical. */
	step = square_size / (grid_size - 1);
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
			dots[i][j] = dots[i][j - 1] + step;
		}
	}

	/* For debugging */
	// int deb_rank = 4;
	// if (rank == deb_rank) {
	// 	printf("X values of proc with rank %d:\n", deb_rank);
	// 	for (int i = 0; i < dots_num[0]; i++)
	// 		printf("%f ", dots[0][i]);
	// 	printf("\n");
	// 	printf("Y values of proc with rank %d:\n", deb_rank);
	// 	for (int i = 0; i < dots_num[1]; i++)
	// 		printf("%f ", dots[1][i]);
	// 	printf("\n");
	// }
}

/* Now, with dots distributed and calculated, run the solver */
void Poisson::Solve() {
	resid_matr = Matrix(dots_num[0], dots_num[1]);
	sol_matr = Matrix(dots_num[0], dots_num[1]);

	/* allocate send & recv buffers */
	for (int i = 0; i < 4; i++) {
		send_buffers[i] = new double[dots_num[(i + 1) % 2]];
		recv_buffers[i] = new double[dots_num[(i + 1) % 2]];
	}
	IsThisProcSendsFirst();

	InitSolMatr();
	CalcResidMatr();
}

/* We initialize sol_matr to Phi on the borders, and to random values in the
 * inner dots.
 */
void Poisson::InitSolMatr() {
	FillBorders(sol_matr, Phi);

	/* Inner part */
	srand(time(NULL));
	for (int i = inner_dots_range[0]; i < inner_dots_range[2]; i++)
		for (int j = inner_dots_range[1]; j < inner_dots_range[3]; j++) {
			/* not sure which values here are better */
			sol_matr(i, j) = rand() % 4096;
		}

	/* For debugging */
	int deb_rank = 2;
	if (rank == deb_rank) {
		printf("sol values of proc with rank %d:\n", deb_rank);
		sol_matr.Print();
	}
}

/*
 * Calculate residuals matr
 */
void Poisson::CalcResidMatr() {
	FillBorders(resid_matr, &zero_filler); /* borders are zeroed */

	ApplyLaplace(sol_matr, resid_matr);

	/* For debugging */
	// int deb_rank = 3;
	// if (rank == deb_rank) {
	// 	printf("resid matr of proc with rank %d:\n", deb_rank);
	// 	resid_matr.Print();
	// }
}


/*
 * Apply discrete Laplace operator to 'matr' and put the result to 'lap_matr'.
 * Border values of 'lap_matr' not touched. Matrices must be already allocated
 * with our usual working size (dots_num[0], dots_num[1]).
 */
void Poisson::ApplyLaplace(const Matrix &matr, Matrix &lap_matr) {
	/* Inner values. I mean, locally inner; for now we doesn't care whether
	 * the borders of this processor are global borders or not
	 */
	for (int i = 1; i < dots_num[0] - 1; i++)
		for (int j = 1; j < dots_num[1] - 1; j++) {
			lap_matr(i, j) = LaplaceFormula(matr(i, j), matr(i - 1, j),
											matr(i + 1, j), matr(i, j - 1),
											matr(i, j + 1));
		}

	ExchangeData(matr);
}

/*
 * Formula of Laplace discrete operator for one dot. Again, everything is
 * square, grid is uniform, which simplifies things.
 */
double Poisson::LaplaceFormula(double center, double left, double right,
							 double bottom, double top) {
	double numerator = 4*center - left - right - bottom - top;
	return numerator / (step*step);
}

/*
 * Exchange data with 4 neighbors needed for Laplace operator.
 * send_buffer, recv_buffers must be allocated.
 * I think that we could get on with less buffers, but who cares.
 */
void Poisson::ExchangeData(const Matrix &matr) {
	/* prepare send_buffers */
	for (int j = 0; j < dots_num[1]; j++) { /* vertical */
		send_buffers[0][j] = matr(0, j);
		send_buffers[2][j] = matr(dots_num[1] - 1, j);
	}
	for (int i = 0; i < dots_num[0]; i++) { /* horizontal */
		send_buffers[1][i] = matr(i, 0);
		send_buffers[3][i] = matr(dots_num[0] - 1, i);
	}

	// I will implement left->right for now
	if (send_first) {
		if (!borders[0]) {
			MPI_Send(send_buffers[0], dots_num[1], MPI_DOUBLE, rank - 1,
					 0, comm);
			if (rank == 2) {
				printf("Successfully sent data from 2 to 1;\n");
			}
		}
	}
	else {
		if (!borders[2]) {
			MPI_Recv(recv_buffers[0], dots_num[1], MPI_DOUBLE, rank + 1,
					 MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);
			if (rank == 1) {
				printf("Successfully got data from 2 on 1:\n");
				for (int j = 0; j < dots_num[1]; j++)
					printf("%f ", recv_buffers[0][j]);
			}
		}
	}
	/* barrier and repeat the same for the other half. */
}

/*
 * Determines whether this processor sends the data first, sets send_first
 */
void Poisson::IsThisProcSendsFirst() {
	if (ranks[1] % 2 == 0)
		send_first = ranks[0] % 2 == 0;
	else
		send_first = ranks[0] % 2 == 1;
	LOG_DEBUG("Proc %d has send_first %d", rank, send_first);
}

/* Fill the (global) borders. Each corner will be set twice, but that's not a
 * big deal and simplifies code a bit.
 */
void Poisson::FillBorders(Matrix &matr, double (*filler)(double x, double y)) {
	int fixed_coord[4] = {0, 0, dots_num[0] - 1, dots_num[1] - 1};
	for (int r = 0; r < 4; r++) {
		if (borders[r]) {
			if (r % 2  == 0) { /* vertical border */
				for (int j = 0; j < dots_num[1]; j++) {
					matr(fixed_coord[r], j) =
						(*filler)(dots[0][fixed_coord[r]], dots[1][j]);
				}
			}
			else { /* horizontal border */
				for (int i = 0; i < dots_num[0]; i++) {
					matr(i, fixed_coord[r]) =
						(*filler)(dots[0][i], dots[1][fixed_coord[r]]);
				}
			}
		}
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
	for (int i = 0; i < 4; i++) {
		if (send_buffers[i])
			delete[] send_buffers[i];
		if (recv_buffers[i])
			delete[] recv_buffers[i];
	}
}
