#pragma once
#include <iostream>
#include <math.h>
#include <float.h>
#include <iomanip>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#define UNUSED(x) (void)(x)

bool is_double(const std::string& s);
int read_ff(const std::string& filename, double* result, int n);
void PrintDouble(double* matrix, int n, int r);
double find_residual(double *a, double *b, double *block, double *sum_array, int n, int m);
double formula(int s, int n, int i, int j);
int inverse_matrix(double *matrix, double *inverse_matrix, int *index, int n, double matrix_norm, int row_ind);
// inline void mult(double *a, double *b, double *res, int m1, int m2, int m3, int m, double norm);
void get_block(double* matrix, double* matrix_block, int n, int m, int i, int j);
void put_block(double *matrix, double *block, int n, int m, int i_block, int j_block);
double norma(double *matrix, int n);
double ravno(double x, double norma);
void ravnoBlock(double *block, int m, double norma);
void pcord(int i, int j);
void swap_rows(double* matr, int k, int l, int n, int m);
void subtraction(double *a, double *b, int i_block, int j_block, int k, int m, int l);
void initialize_matrix(double* solution, int n);
void copyMatrix(double* in, double* out, int current_size, int req_size);
void setZeroBlock(double *matrix, double *block, int i_block, int j_block, int n, int m);
int solve(int n, int m, double* matr, double* block, double* solution, double* inverse, double* tmp, double* block1, /*double* block2,*/ double matrix_norm, int *block_index);

inline void mult(double *a, double *b, double *res, int m1, int m2, int m3, int m, double norm)
{
	int t = 0, q = 0, r = 0;
	int v = m1, h = m3, ah = m2;
	int v3 = v % 3, h3 = h % 3;
	double s00 = 0, s01 = 0, s02 = 0;
	double s10 = 0, s11 = 0, s12 = 0;
	double s20 = 0, s21 = 0, s22 = 0;
    UNUSED(norm);
    int count_b = 0;
    for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < m; j++)
            {
                res[i * m + j] = 0.0;
                double local = fabs(b[i * m + j]);
                if ( 1e+250 * norm < local || local < 1e-100 * norm)
                {
                    b[i * m + j] = 0.;
                    count_b++;
                }
            }
        }
    if(count_b == m*m) {
        return;
    }
	for (r = 0; r < v3; r++)
	{
		for (t = 0; t < h3; t++)
		{
			double sum;
			sum = 0;
			for (q = 0; q < ah; q++)
				sum += a[r * m + q] * b[q * m + t];
			res[r * m + t] += sum;
		}
		for (; t < h; t += 3)
		{
			s00 = 0;
			s01 = 0;
			s02 = 0;
			for (q = 0; q < ah; q++)
			{
				s00 += a[r * m + q] * b[q * m + t];
				s01 += a[r * m + q] * b[q * m + t + 1];
				s02 += a[r * m + q] * b[q * m + t + 2];
			}
			res[r * m + t] += s00;
			res[r * m + t + 1] += s01;
			res[r * m + t + 2] += s02;
		}
	}
	for (; r < v; r += 3)
	{
		for (t = 0; t < h3; t++)
		{
			s00 = 0;
			s10 = 0;
			s20 = 0;
			for (q = 0; q < ah; q++)
			{
				s00 += a[r * m + q] * b[q * m + t];
				s10 += a[(r + 1) * m + q] * b[q * m + t];
				s20 += a[(r + 2) * m + q] * b[q * m + t];
			}
			res[r * m + t] += s00;
			res[(r + 1) * m + t] += s10;
			res[(r + 2) * m + t] += s20;
		}
		for (; t < h; t += 3)
		{
			s00 = 0;
			s01 = 0;
			s02 = 0;
			s10 = 0;
			s11 = 0;
			s12 = 0;
			s20 = 0;
			s21 = 0;
			s22 = 0;
			for (q = 0; q < ah; q++)
			{
				s00 += a[r * m + q] * b[q * m + t];
				s01 += a[r * m + q] * b[q * m + t + 1];
				s02 += a[r * m + q] * b[q * m + t + 2];
				s10 += a[(r + 1) * m + q] * b[q * m + t];
				s11 += a[(r + 1) * m + q] * b[q * m + t + 1];
				s12 += a[(r + 1) * m + q] * b[q * m + t + 2];
				s20 += a[(r + 2) * m + q] * b[q * m + t];
				s21 += a[(r + 2) * m + q] * b[q * m + t + 1];
				s22 += a[(r + 2) * m + q] * b[q * m + t + 2];
			}
			res[r * m + t] += s00;
			res[r * m + t + 1] += s01;
			res[r * m + t + 2] += s02;
			res[(r + 1) * m + t] += s10;
			res[(r + 1) * m + t + 1] += s11;
			res[(r + 1) * m + t + 2] += s12;
			res[(r + 2) * m + t] += s20;
			res[(r + 2) * m + t + 1] += s21;
			res[(r + 2) * m + t + 2] += s22;
		}
	}
}