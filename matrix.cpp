#include <stdio.h>

#include "matrix.hpp"

Matrix::Matrix(int _n, int _m) : m(_m), n(_n) {
	refcounter = new int(1);
	data = new double[n*m];
}

Matrix::~Matrix() {
	(*refcounter)--;
	if (refcounter == 0) {
		delete refcounter;
		delete[] data;
	}
}

Matrix::Matrix(const Matrix &matr) : n(matr.n), m(matr.m) {
	refcounter = matr.refcounter;
	(*refcounter)++;
	data = matr.data;
}

Matrix& Matrix::operator=(const Matrix &matr) {
	n = matr.n;
	m = matr.m;
	refcounter = matr.refcounter;
	(*refcounter)++;
	data = matr.data;
	return *this;
}

void Matrix::Print() const {
	for (int j = m - 1; j >= 0; j--) {
		for (int i = 0; i < n; i++) {
			printf("%f ", (*this)(i, j));
		}
		printf("\n");
	}
}

void Matrix::GetRow(int j, double *buf) {
	for (int p = 0; p < n; p++)
		buf[p] = (*this)(p, j);
}

void Matrix::GetColumn(int i, double *buf) {
	for (int q = 0; q < m; q++)
		buf[q] = (*this)(i, q);
}
