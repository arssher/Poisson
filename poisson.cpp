#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <algorithm>
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "lib/macrologger.h"
#include "poisson.hpp"

/* To use pragma inside macros */
#define OMP_PARA_FOR _Pragma("omp parallel for schedule(static)")

static double zero_filler(double x, double y) { return 0.0; };

Poisson::Poisson(double x0, double y0, double square_size, int grid_size,
				 int sdi_it,
				 double (*_F)(double x, double y),
				 double (*_Phi)(double x, double y),
				 char *dump_d)
	: F(_F), Phi(_Phi), dots_per_proc(NULL),
	  sdi_iterations(sdi_it), eps(0.0001), dump_dir(dump_d) {
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
	LOG_INFO_MASTER("Size is %d\n", size);
#ifdef _OPENMP
	LOG_INFO_MASTER("OpenMP is ok, max num of threads: %d",
					omp_get_max_threads());
#else
	LOG_INFO_MASTER("OpenMP is not supported and will not be used. %d", 0);
#endif

	/* define processors grid and block unneeded processors */
	// if (size < 4) {
		// throw PoissonException("Error: we require at least 4 processors");
	// }
	proc_grid_size = static_cast<int>(floor(sqrt(size)));
	LOG_INFO_MASTER("Proc grid size is %d\n", proc_grid_size);

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
	LOG_INFO_MASTER("Square starts at (%f, %f) of size %f, grid size %d\n",
					 x0, y0, square_size, grid_size);

	time_t start;
	if (rank == 0) {
		start = time(NULL);
	}
	try {
		DistributeDots(grid_size);
		CalculateDots(x0, y0, square_size, grid_size);
		Solve();
	}
	catch (PoissonException &e) {
		/* To synchronize with unused processes */
		MPI_Barrier(MPI_COMM_WORLD);
		throw;
	}
	LOG_INFO_MASTER("It took me %.2f seconds to compute",
					(double)(time(NULL) - start));

	if (dump_dir)
		DumpSolution();

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
	LOG_DEBUG_MASTER("Step is %.17g", step);
	/* bottom left point of this processor working area */
	double bottom_left[2] = {x0, y0};
	for (int r = 0; r < 2; r++) {
		for (int i = 0; i < ranks[r]; i++) {
			/* faster, but differs a bit from sequential execution */
			bottom_left[r] += step*(dots_per_proc[i]);
			/* slower, but the same with sequential execution */
			// for (int j = 0; j < dots_per_proc[i]; j++)
				// bottom_left[r] += step;
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
}

/* Now, with dots distributed and calculated, run the solver */
void Poisson::Solve() {
	double global_error;
	resid_matr = Matrix(dots_num[0], dots_num[1]);
	sol_matr = Matrix(dots_num[0], dots_num[1]);
	tmp_matr = Matrix(dots_num[0], dots_num[1]);

	/* allocate send & recv buffers */
	for (int i = 0; i < 4; i++) {
		send_buffers[i] = new double[dots_num[(i + 1) % 2]];
		recv_buffers[i] = new double[dots_num[(i + 1) % 2]];
	}
	IsThisProcSendsFirst();
	InitSolMatr();

	for (int i = 0; i < sdi_iterations; i++) {
		double local_error = SteepDescentIteration();
		MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_MAX,
					  comm);
		LOG_INFO_MASTER("Steep descent iteration %d done, error %f",
						i, global_error);
	}

	/* Now run CGM until convergence */
	/* CGM initialization: g(0) = resid(0) */
	g_matr = resid_matr.DeepCopy();

	global_error = eps + 1948;
	int CGM_its = 0;
	while (global_error > eps) {
		double local_error = CGMIteration();
		MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_MAX,
					  comm);
		LOG_INFO_MASTER("CGM iteration %d done, error %f", CGM_its, global_error);
		CGM_its++;
	}
	LOG_INFO_MASTER("Done! %d", 0);
}

/* One iteration of SteepDescent method. The result is in sol_matr.
 * Returns the error corresponding to max norm.
 */
double Poisson::SteepDescentIteration() {
	double error = 0.0;
	double new_sol_val;

	CalcResidMatr();
	double tau = CalcTauSteepDescent();
	LOG_DEBUG_MASTER("SD Tau is %.17g", tau);
	/* There should be the following pragma, but BlueGene\P doesn't support it
	 * either
	 */
	/* #pragma omp parallel for schedule(static) private(new_sol_val) reduction(max:error) */
	int j;
	double thread_local_error = 0.0;
	#pragma omp parallel private(j, new_sol_val) firstprivate(thread_local_error)
	{
		# pragma omp for schedule(static)
		for (int i = 0; i < dots_num[0]; i++)
			for (j = 0; j < dots_num[1]; j++) {
				new_sol_val = sol_matr(i, j) - tau * resid_matr(i, j);
				thread_local_error =
					std::max(thread_local_error, fabs(sol_matr(i, j) - new_sol_val));
				sol_matr(i, j) = new_sol_val;
			}
		#pragma omp critical
		{
			error = std::max(error, thread_local_error);
		}
	}

	return error;
}

/* One iteration of CGM method. The result is in sol_matr.
 * Returns error corresponding to max norm.
 */
double Poisson::CGMIteration() {
	double error = 0.0;
	double new_sol_val;

	CalcResidMatr();
	double alpha = CalcAlphaCGM();
	LOG_DEBUG_MASTER("CGM alpha is %.17g", alpha);
	int j;
	#pragma omp parallel for schedule(static) private(j)
	for (int i = 0; i < dots_num[0]; i++)
		for (j = 0; j < dots_num[1]; j++) {
			g_matr(i, j) = resid_matr(i, j) - alpha * g_matr(i, j);
		}
	double tau = CalcTauCGM();
	LOG_DEBUG_MASTER("CGM tau is %.17g", tau);

	double thread_local_error = 0.0;
	#pragma omp parallel private(j, new_sol_val) firstprivate(thread_local_error)
	{
		# pragma omp for schedule(static)
		for (int i = 0; i < dots_num[0]; i++)
			for (j = 0; j < dots_num[1]; j++) {
				new_sol_val = sol_matr(i, j) - tau * g_matr(i, j);
				thread_local_error =
					std::max(thread_local_error, fabs(sol_matr(i, j) - new_sol_val));
				sol_matr(i, j) = new_sol_val;
			}
		#pragma omp critical
		{
			error = std::max(error, thread_local_error);
		}
	}

	return error;
}

/* We initialize sol_matr to Phi on the borders, and to random values in the
 * inner dots.
 */
void Poisson::InitSolMatr() {
	FillBorders(sol_matr, Phi);

	/* Inner part */
	srand(time(NULL) + rank);
	/* better to use collapse in all such places, but it is not supported by
	   BlueGene */
	/* #pragma omp parallel for collapse(2) */
	int j;
	#pragma omp parallel for schedule(static) private(j)
	for (int i = inner_dots_range[0]; i < inner_dots_range[2]; i++)
		for (j = inner_dots_range[1]; j < inner_dots_range[3]; j++) {
			/* not sure which values here are better */
			// sol_matr(i, j) = (rand() % 4096 )/ 4096.0;
			sol_matr(i, j) = 0.0;
		}
}

/*
 * Calculate residuals matr, resid_matr = (delta sol_matr) - F
 */
void Poisson::CalcResidMatr() {
	FillBorders(resid_matr, &zero_filler); /* borders are zeroed */

	ApplyLaplace(sol_matr, resid_matr);

	/* substact F(x, y) */
	int j;
	#pragma omp parallel for schedule(static) private(j)
	for (int i = inner_dots_range[0]; i < inner_dots_range[2]; i++)
		for (j = inner_dots_range[1]; j < inner_dots_range[3]; j++) {
			resid_matr(i, j) -= (*F)(dots[0][i], dots[1][j]);
		}
}

#define LAPLACE_COMPUTE_BORDER(r, max_i, center_ind, left, right, bottom, top) do { \
		if (!borders[r]) { \
			assert(recv_buffer_used[(r + 2) % 4]); \
			OMP_PARA_FOR \
			for(int i = 1; i < max_i; i++) { \
				lap_matr center_ind = \
					LaplaceFormula(matr center_ind, left, right, bottom, top); \
			} \
		} \
	} while (0)

/*
 * Apply discrete Laplace operator to 'matr' and put the result to 'lap_matr'.
 * Border values of 'lap_matr' not touched. Matrices must be already allocated
 * with our usual working size (dots_num[0], dots_num[1]).
 */
void Poisson::ApplyLaplace(const Matrix &matr, Matrix &lap_matr) {
	/* Inner values. I mean, locally inner; for now we doesn't care whether
	 * the borders of this processor are global borders or not
	 */
	int j;
	#pragma omp parallel for schedule(static) private(j)
	for (int i = 1; i < dots_num[0] - 1; i++)
		for (j = 1; j < dots_num[1] - 1; j++) {
			lap_matr(i, j) = LaplaceFormula(matr(i, j),
											matr(i - 1, j), matr(i + 1, j),
											matr(i, j - 1), matr(i, j + 1));
		}

	/* fill recv_buffers */
	ExchangeData(matr);
	int deb_rank = size > 1 ? 6 : 0;

	/* now, the borders, without corners */
	int max_x = dots_num[0] - 1;
	int max_y = dots_num[1] - 1;
	/* left border */
	LAPLACE_COMPUTE_BORDER(0, max_y, (0, i),
						   recv_buffers[2][i], matr(1, i),
						   matr(0, i - 1), matr(0, i + 1));
	/* bottom border */
	LAPLACE_COMPUTE_BORDER(1, max_x, (i, 0),
						   matr(i - 1, 0), matr(i + 1, 0),
						   recv_buffers[3][i], matr(i, 1));
	/* right border */
	LAPLACE_COMPUTE_BORDER(2, max_y, (max_x, i),
						   matr(max_x - 1, i), recv_buffers[0][i],
						   matr(max_x, i - 1), matr(max_x, i + 1));

	/* top border */
	LAPLACE_COMPUTE_BORDER(3, max_x, (i, max_y),
						   matr(i - 1, max_y), matr(i + 1, max_y),
						   matr(i, max_y - 1), recv_buffers[1][i]);

	/* and the corners */
	if (!borders[0] && !borders[1])
		/* bottom left */
		lap_matr(0, 0) = LaplaceFormula(matr(0, 0),
										recv_buffers[2][0], matr(1, 0),
										recv_buffers[3][0], matr(0, 1));
	if (!borders[0] && !borders[3])
		/* top left */
		lap_matr(0, max_y) =
			LaplaceFormula(matr(0, max_y),
						   recv_buffers[2][max_y], matr(1, max_y),
						   matr(0, max_y - 1), recv_buffers[1][0]);

	if (!borders[2] && !borders[1])
		/* bottom right */
		lap_matr(max_x, 0) =
			LaplaceFormula(matr(max_x, 0),
						   matr(max_x - 1, 0), recv_buffers[0][0],
						   recv_buffers[3][max_x], matr(max_x, 1));

	if (!borders[2] && !borders[3])
		/* top right */
		lap_matr(max_x, max_y) =
			LaplaceFormula(matr(max_x, max_y),
						   matr(max_x - 1, max_y), recv_buffers[0][max_y],
						   matr(max_x, max_y - 1), recv_buffers[1][max_x]);
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
 * It takes every local border of 'matr' and if it is not a global border,
 * sends it to neighbors; data is put to recv_buffers.
 * send_buffer, recv_buffers must be allocated.
 * I think that we could get on with less buffers, but who cares.
 */
void Poisson::ExchangeData(const Matrix &matr) {
	/* reset recv_buffer_used */
	for (int r = 0; r < 4; r++)
		recv_buffer_used[r] = false;

	/* prepare send_buffers */
	#pragma omp parallel for schedule(static)
	for (int j = 0; j < dots_num[1]; j++) { /* vertical */
		send_buffers[0][j] = matr(0, j);
		send_buffers[2][j] = matr(dots_num[0] - 1, j);
	}
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < dots_num[0]; i++) { /* horizontal */
		send_buffers[1][i] = matr(i, 0);
		send_buffers[3][i] = matr(i, dots_num[1] - 1);
	}

	int dest_ranks[4] = {rank - 1, rank - proc_grid_size, rank + 1, rank +
						 proc_grid_size};
	int src_ranks[4] = {rank + 1, rank + proc_grid_size, rank - 1, rank -
						proc_grid_size};
	for (int r = 0; r < 4; r++) {
		int buf_size = dots_num[(r + 1) % 2];
		SendDataOneDirection(r, buf_size, src_ranks, dest_ranks);
		send_first = !send_first;
		SendDataOneDirection(r, buf_size, src_ranks, dest_ranks);
		send_first = !send_first;
	}
}

/*
 * Working horse for ExchangeData, sends data in one direction from half
 * of the nodes (with those having sends_first == true) to the other half.
 * r, index in arrays such as borders[4] sets the direction.
 */
void Poisson::SendDataOneDirection(int r, int buf_size, int *src_ranks,
									int* dest_ranks) {
	if (send_first) {
		/* we send the data iff we have no global border on that side */
		if (!borders[r]) {
			MPI_Send(send_buffers[r], buf_size, MPI_DOUBLE,
					 dest_ranks[r], 0, comm);
		}
	}
	else {
		/* and receive it iff we have no global border on the receiving side */
		if (!borders[(r + 2) % 4]) {
			MPI_Status recv_status;
			recv_buffer_used[r] = true;
			MPI_Recv(recv_buffers[r], buf_size, MPI_DOUBLE,
					 src_ranks[r], MPI_ANY_TAG, comm, &recv_status);
			assert(recv_status.MPI_SOURCE == src_ranks[r]);
		}
	}
	MPI_Barrier(comm);
}

/*
 * Determines whether this processor sends the data first, sets send_first
 */
void Poisson::IsThisProcSendsFirst() {
	if (ranks[1] % 2 == 0)
		send_first = ranks[0] % 2 == 0;
	else
		send_first = ranks[0] % 2 == 1;
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


/* Calculate tau in steep descent method.
 * resid_matr must be calculated at this point and tmp_matr allocated
 */
double Poisson::CalcTauSteepDescent() {
	double numerator = resid_matr.ScalarProduct(resid_matr, step);

	FillBorders(tmp_matr, &zero_filler); /* in principle this is not necessary */
	ApplyLaplace(resid_matr, tmp_matr);

	double denominator = tmp_matr.ScalarProduct(resid_matr, step);

	SumTwoDoublesGlobally(numerator, denominator);
	LOG_DEBUG_MASTER("Collected numerator %.17g", numerator);
	LOG_DEBUG_MASTER("Collected denominator %.17g", denominator);
	if (fabs(denominator) < 10e-7)
		throw PoissonException("Error: Denominator close to zero in CalcTauSteepDescent\n");
	return numerator / denominator;
}

/* Calculate alpha in CGM method
 * resid_matr must be calculated at this point and tmp_matr allocated
 */
double Poisson::CalcAlphaCGM() {
	FillBorders(tmp_matr, &zero_filler); /* in principle this is not necessary */
	ApplyLaplace(resid_matr, tmp_matr);
	double numerator = tmp_matr.ScalarProduct(g_matr, step);

	ApplyLaplace(g_matr, tmp_matr);
	double denominator = tmp_matr.ScalarProduct(g_matr, step);

	SumTwoDoublesGlobally(numerator, denominator);
	LOG_DEBUG_MASTER("Collected alpha numerator %f", numerator);
	LOG_DEBUG_MASTER("Collected alpha denominator %f", denominator);
	if (fabs(denominator) < 10e-7)
		throw PoissonException("Error: Denominator close to zero in CalcAlphaCGM\n");
	return numerator / denominator;
}

/* Calculate tau in CGM method.
 * resid_matr and g_matr must be calculated at this point and tmp_matr allocated
 */
double Poisson::CalcTauCGM() {
	double numerator = resid_matr.ScalarProduct(g_matr, step);

	FillBorders(tmp_matr, &zero_filler); /* in principle this is not necessary */
	ApplyLaplace(g_matr, tmp_matr);
	double denominator = tmp_matr.ScalarProduct(g_matr, step);

	SumTwoDoublesGlobally(numerator, denominator);
	if (fabs(denominator) < 10e-7)
		throw PoissonException("Error: Denominator close to zero in CalcTauCGM\n");

	return numerator / denominator;
}

/* Does Allreduce for CalcTau and CalcAlpha. Receives local num and denom,
 * sums them globally and saves the result in place
 */
void Poisson::SumTwoDoublesGlobally(double &numerator, double &denominator) {
	double numer_denom_local[2] = { numerator, denominator };
	double numer_denom_global[2];

	MPI_Allreduce(numer_denom_local, numer_denom_global, 2, MPI_DOUBLE, MPI_SUM,
				  comm);
	numerator = numer_denom_global[0];
	denominator = numer_denom_global[1];
}

/* Dump the solution of this processor to 'dump_dir'/solution.'rank'.csv */
void Poisson::DumpSolution() {
	assert(dump_dir);
	LOG_INFO_MASTER("Starting dumping...%d", 0);

	if (rank == 0) {
		struct stat st = {0};
		if (stat(dump_dir, &st) == -1) {
			mkdir(dump_dir, 0700);
		}
	}
	/* Wait until directory will be created */
	MPI_Barrier(comm);

	std::stringstream ss;
	ss << dump_dir << "/solution." << rank << ".csv";
	const char *fname = ss.str().c_str();

	FILE *f = fopen(fname, "w");
	if (f == NULL) {
		throw PoissonException("Failed to open dump file");
	}
	for (int i = 0; i < dots_num[0]; i++)
		for (int j = 0; j < dots_num[1]; j++) {
			fprintf(f, "%.17g, %.17g, %.17g\n", dots[0][i], dots[1][j], sol_matr(i, j));
		}
	fclose(f);
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
