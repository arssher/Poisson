#include <stdio.h>
#include "assert.h"

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

Matrix Matrix::DeepCopy() const {
	Matrix matr(n, m);
	#pragma omp parallel for schedule (static)
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++) {
			matr(i, j) = (*this)(i, j);
		}

	return matr;
}

void Matrix::Print() const {
	for (int j = m - 1; j >= 0; j--) {
		for (int i = 0; i < n; i++) {
			printf("%019.17g  ", (*this)(i, j));
		}
		printf("\n");
	}
}

void Matrix::GetRow(int j, double *buf) {
	#pragma omp parallel for schedule (static)
	for (int p = 0; p < n; p++)
		buf[p] = (*this)(p, j);
}

void Matrix::GetColumn(int i, double *buf) {
	#pragma omp parallel for schedule (static)
	for (int q = 0; q < m; q++)
		buf[q] = (*this)(i, q);
}

double Matrix::ScalarProduct(const Matrix &matr2, double step) const {
	assert(n == matr2.n);
	assert(m == matr2.m);

	double sp = 0.0;
	#pragma omp parallel for schedule (static) reduction(+:sp)
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++) {
			sp += (*this)(i, j) * matr2(i, j) * step * step;
		}

	return sp;
}
