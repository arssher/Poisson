#include <stdio.h>

#include "matrix.hpp"

Matrix::Matrix(int _m, int _n) : m(_m), n(_n) {
	refcounter = new int(1);
	data = new double[m*n];
}

Matrix::~Matrix() {
	(*refcounter)--;
	if (refcounter == 0) {
		delete refcounter;
		delete[] data;
	}
}

Matrix::Matrix(const Matrix &matr) : m(matr.m), n(matr.n) {
	refcounter = matr.refcounter;
	(*refcounter)++;
	data = matr.data;
}

Matrix& Matrix::operator=(const Matrix &matr) {
	m = matr.m;
	n = matr.n;
	refcounter = matr.refcounter;
	(*refcounter)++;
	data = matr.data;
}

void Matrix::Print() const {
	for (int i = m - 1; i >= 0; i--) {
		for (int j = 0; j < n; j++) {
			printf("%f ", (*this)(i, j));
		}
		printf("\n");
	}
}

void Matrix::GetRow(int i, double *buf) {
	for (int p = 0; p < n; p++)
		buf[p] = (*this)(i, p);
}

void Matrix::GetColumn(int j, double *buf) {
	for (int q = 0; q < m; q++)
		buf[q] = (*this)(q, j);
}
