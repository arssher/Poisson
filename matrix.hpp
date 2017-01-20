/* Dummy matrix of doubles with m rows and n columns, growing
 * from the bottom left corner.
 */
class Matrix {
public:
	/* = 0 to make Matrix() constructor available */
	Matrix(int _n = 0, int _m = 0);
	~Matrix();
	Matrix(const Matrix &matr);
	Matrix& operator=(const Matrix &matr);
	/* return element on column i and row j */
	double operator()(int i, int j) const { return data[j*n + i]; }
	double& operator()(int i, int j) { return data[j*n + i]; }
	/* The caller must allocate memory himself */
	void GetRow(int j, double *buf);
	void GetColumn(int i, double *buf);
	/* Scalar product this*matr2 as defined in the manual */
	double ScalarProduct(const Matrix &matr2, double step) const;
	void Print() const;
private:
	int n; /* cols */
	int m; /* rows */
	double *data;
	int *refcounter;
};
