#include "solve.hpp"

void add(double *a, double *b, int n, int m)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            a[i*m+j] = fabs(a[i*m+j]) + fabs(b[i*m+j]);
}

double norma_block (double* matr, double* block, double* block1, int n, int m) 
{
    int k = n/m;
    int bl = n - (k*m);
    int fsize;
    double max = 0.;
    if(bl != 0) {
        fsize = k+1;
    } else {
        fsize = k;
    }

    for (int i =0; i < fsize; i++) {
        get_block(matr, block, n , m, 0 , i);
        for (int j = 0; j < fsize; j++) {
            get_block(matr, block1, n ,m, j, i);
            add(block, block1, m , m);
            double tmp = norma(block, m);
            if (tmp > max) {
                max = tmp;
            }
        }
    }
    return max;
}

double residual_matrix(double* matr, double* solution, double* tmp, double* block, double* block1, int n, int m) {
    mult(solution, matr, tmp, n, n, n, n);
    for (int i =0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            tmp[i*n+j] -= (i == j ? 1 : 0);
        }
    }
    return norma_block(tmp, block, block1, n ,m);
}
double matrix_residual(double* matr, double* solution, double* tmp, double* block, double* block1, int n, int m) {
    mult(matr, solution, tmp, n, n, n, n);
    for (int i =0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            tmp[i*n+j] -= (i == j ? 1 : 0);
        }
    }
    return norma_block(tmp, block, block1, n ,m);
}
