#include <exception>

#include "matrix.hpp"

class PoissonException : public std::exception {
public:
    PoissonException(const char *_msg): msg(_msg) {}

    virtual const char* what() const throw() {
        return msg;
    }

private:
	const char *msg;
};

/* Dummy solver of Poisson equation Dirichlet problem.
 * Works only with square grids and uses only square number of processors.
 * It is really not hard to implement rectanges as well, I just don't need it.
 * We index both processors (ranks[0], rank[1]) and dots indexes from the left
 * bottom corner, just like usual axis.
 * Class can be used only once, i.e. one instance -- one computation.
 * I use class just as a scope for some global variables, for automatic
 * destructor call and for virtual F and Phi.
 */
class Poisson {
private:
	/* data describing grid of processors */
	int rank; /* this processor MPI rank */
	int size; /* total number of processors */
	/* processors grid consists of proc_grid_size x proc_grid_size processors */
	int proc_grid_size;
	MPI_Comm comm; /* square processors grid communicator */
	/* this processor location on [x, y] axis, starting from 0 */
	int ranks[2];

	/* data describing grid of dots, global and for this processor */
	/* Grid step. Grid is uniform and symmetrical. */
	double step;
	/* Number of dots per process along either x or y axis,
     * array of proc_grid_size size.
     * dots_per_proc[i] = (how many x dots proc with ranks[0] = i has) =
     * (how many y dots proc with ranks[1] = i has)
     */
	int *dots_per_proc;
	/* Number of dots along [x, y] axis on this specific processor */
	int dots_num[2];
	/* Dots coordinates [x, y] themselves */
	double *dots[2];
	/* this processor touches the [left, bottom, right, top] border? */
	bool borders[4];
    /* These 4 indexes (x_min, y_min, x_max, y_max) show range of inner dots
	 * area part for this processor. Note that 'inner' here means globally
	 * inner, i.e. out of global border. In the simplest case, when
	 * the processor doesn't touch any borders, it will be (0, 0, dots_num[0],
	 * dots_num[1]).  If, for instance, it touches the right border and the
	 * bottom border, it will be (0, 1, dots_num[0] - 1, dots_num[1]). We use
	 * it to avoid ugly '-1' and '+1' indexes while looping over the inner
	 * part.
	 */
	int inner_dots_range[4];
	/* Send buffers for messages going to the left, bottom, right and top from
     * this processor.
	 */
	double *send_buffers[4];
	/* receive buffers for messages going from right, top, left, bottom to
     * to this processor.
	 */
	double *recv_buffers[4];
	/* we need to mark matrix in checkerboard order to avoid deadlocks */
	bool send_first;

	/* residuals */
	Matrix resid_matr;
	/* solution, 'p' function in the manual */
	Matrix sol_matr;
	/* Used in thau calculations. Store it here to avoid reallocating
     * every time
	 */
	Matrix tmp_matr;

	int sdi_iterations;

	/* Right hand side of Poisson equation */
	double (*F)(double x, double y);
	/* boundary function */
	double (*Phi)(double x, double y);

	void DistributeDots(int grid_size);
	void CalculateDots(double x0, double y0, double square_size, int grid_size);
	void Solve();
	void InitSolMatr();
	void CalcResidMatr();
	void ApplyLaplace(const Matrix &matr, Matrix &lap_matr);
	double LaplaceFormula(double center, double left, double right,
						  double bottom, double top);
	void ExchangeData(const Matrix &matr);
	void SendDataOneDirection(int r, int buf_size, int* src_ranks,
							  int *dest_ranks);
	void IsThisProcSendsFirst();
	void FillBorders(Matrix &matr, double (*filler)(double x, double y));
	double CalcTauSteepDescent();
	double SteepDescentIteration();
public:
   /* (x0, y0) and square_size define the square we are working on.
    * grid_size*grid_size is the total number of dots
    */
	Poisson(double x0, double y0, double square_size, int grid_size,
		    int sdi_it,
			double (*_F)(double x, double y),
			double (*_Phi)(double x, double y));
	~Poisson();
};
