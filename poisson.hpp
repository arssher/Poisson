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
 */
class Poisson {
private:
	int rank; /* this processor MPI rank */
	int size; /* total number of processors */
	/* processors grid consists of proc_grid_size x proc_grid_size processors */
	int proc_grid_size;
	MPI_Comm comm; /* square processors grid communicator */
	/* this processor location on [x, y] axis, starting from 0 */
	int ranks[2];

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

	void DistributeDots(int grid_size);
	void CalculateDots(double x0, double y0, double square_size, int grid_size);
public:
   /* (x0, y0) and square_size define the square we are working on.
    * grid_size*grid_size is the total number of dots
    */
	Poisson(double x0, double y0, double square_size, int grid_size);
	~Poisson();

	/* Right hand side of Poisson equation */
	virtual double F(double x, double y) const = 0;


};
