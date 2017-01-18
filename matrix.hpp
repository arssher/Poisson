/* Dummy m x n matrix of doubles */
class Matrix {
public:
	Matrix(int m, int n);
	~Matrix();
	Matrix(const Matrix &matr);
	Matrix& operator=(const Matrix &matr);
	double operator()(int i, int j) const { return data[i*n + j]; }
	double& operator()(int i, int j) { return data[i*n + j]; }
	/* The caller must allocate memory himself */
	void GetRow(int i, double *buf);
	void GetColumn(int j, double *buf);
	void Print() const;
private:
	int m;
	int n;
	double *data;
	int *refcounter;
};
