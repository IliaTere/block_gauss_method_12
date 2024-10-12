#pragma once
#include <iostream>
#include <math.h>
#include <float.h>
#include "mult.hpp"
#include <iomanip>
#include <cstring>
//#include "reader.hpp"
#define UNUSED(x) (void)(x)
#define eps -DBL_MAX

using namespace std;

void get_block (double *matr, double *block , int n, int m, int i , int j )
{  
    for (int i=0; i < m*m; i++) 
        block[i] = 0;
    int ii=0,jj=0;
    if(m == 0 ) {m = 1;}
    int k = n / m;
    int l = n - k * m;
    int v = (j < k ? m : l), h = (i < k ? m : l);
    double *ma_bl = matr + i * n * m+ j * m;
    for (ii = 0; ii < h; ii++)
    {
        for (jj = 0; jj < v; jj++)
        {
            block[ii *m+jj]=ma_bl[ii *n+jj];
        }   
    }
}
void put_block (double *matr, double *block , int n, int m, int i , int j )
{
int ii=0,jj=0;
int k = n / m;
int l = n - k * m;
int v = (j < k ? m : l), h = (i < k ? m : l);
double *ma_bl = matr + i * n * m+ j * m;
for (ii = 0; ii < h; ii++)
{
    for (jj = 0; jj < v; jj++)
    {
        ma_bl[ii *n+jj]=block[ii *m+jj];
    }   
}
}

double norma(double* block, int m) {
    int i,j;
    double norm_m=0, temp;
    for(i=0;i<m;i++)  
    {
        temp=0;
        for(j=0;j<m;j++)
            temp+= fabs(block[j*m + i]);
        if(temp>norm_m)
            norm_m=temp;
    }
    return norm_m;
} 

int treug(double * a, double * b, int n, double norma, double* c) {
    int i;
    int j;
    int k;
    int t;
    double p;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            b[i * n + j] = 0;
        }
        for (int i = 0; i < n; i++) {
            b[i * n + i] = 1;
        }
    }

    for (i = 0; i < n; i++) {
        t = -1;
        for (j = i; j < n; j++)
            if (a[j * n + i] >  5e-15 * norma || a[j * n + i] < -5e-15 * norma) {
                //printf("%e\n", a[j * n + i]);
                t = 1;
            }
        if (t == -1) {
            return -1;
        }

        p = a[i * n + i];
        t = i;

        for (j = i; j < n; j++) {
            if (fabs(a[j * n + i]) > fabs(p)) {
                p = a[j * n + i];
                t = j;
            }
        }

        for (j = 0; j < n; j++) c[j] = a[t * n + j];
        for (j = 0; j < n; j++) a[t * n + j] = a[i * n + j];
        for (j = 0; j < n; j++) a[i * n + j] = c[j];

        for (j = 0; j < n; j++) c[j] = b[t * n + j];
        for (j = 0; j < n; j++) b[t * n + j] = b[i * n + j];
        for (j = 0; j < n; j++) b[i * n + j] = c[j];

        p = a[i * n + i];

        for (j = 0; j < n; j++) {
            a[i * n + j] = a[i * n + j] / p;
            b[i * n + j] = b[i * n + j] / p;
        }

        for (j = i + 1; j < n; j++) {
            p = a[j * n + i];
            for (k = 0; k < n; k++) {
                a[j * n + k] = a[j * n + k] - a[i * n + k] * p;
                b[j * n + k] = b[j * n + k] - b[i * n + k] * p;
            }
        }

    }
    return 1;
}

void diag(double * a, double * b, int n) {
    int i, j, k;
    double p;
    for (i = n - 1; i > 0; i--) {
        for (j = 0; j < i; j++) {
            p = a[j * n + i];
            for (k = 0; k < n; k++) {
                a[j * n + k] = a[j * n + k] - a[i * n + k] * p;
                b[j * n + k] = b[j * n + k] - b[i * n + k] * p;
            }
        }
    }
}

double ravno(double x, double norma)
{
    if(fabs(x)<1e-7*norma){
        return 0.0;}
    else {return x;}
}

void ravnoBlock(double *block, int m, double norma)
{
    int i,j;
    for(i=0;i<m;i++)
    {
        for(j=0;j<m;j++)
        {
            block[j*m + i] = ravno(block[j*m + i], norma);
        }
    }
}
void pcord(int i, int j)
{
    printf("(%d, %d)\n", i, j);
}

int findmax(double* matr, double* block, int n, int m,int l, int j, double norm, double* buffer) {
    int k = -1;
    int t;
    double tmp;
    double max_n = DBL_MAX;
    double* inverse = new double[m * m];
    for (int i = l; i < n/m; i++) {
        get_block(matr, block, n ,m , i, j);
        t = treug(block, inverse, m, norm, buffer);
        if (t==-1) {
            continue;
        }
        diag(block,inverse,m);
        tmp = norma(inverse, m);
        if (tmp < max_n) {
            max_n = tmp;
            k = i;
        }
        
    }
    delete [] inverse;
    return k;
}


void swap_rows(double* matr, int k, int l, int n, int m) {
    double* tmp = new double[m*n];
    for (int i = 0; i < m*n; i++) {
        tmp[i] = matr[k*m * n + i];
    }
    for (int i = 0; i < m*n; i++) {
        matr[m*k * n + i] = matr[m*l * n + i];
    }
    for (int i = 0; i < m*n; i++) {
        matr[m*l * n + i] = tmp[i];
    }
    delete[] tmp;
}

void subtraction(double* a, double* b, int m) {
    for (int i = 0; i < m*m; i++) {
        a[i] -= b[i];
    }    
}
void initialize_matrix(double* solution, int n) {
    std::memset(solution, 0, n * n * sizeof(double));

    for (int i = 0; i < n; ++i) {
        solution[i * n + i] = 1.0;
    }
}

int solve(int n, int m, double* matr, double* block, double* solution, double* inverse, double* tmp, double* block1, double* block2, double norma) {
    int k = n / m;
    int bl = n - (k * m);
    int l = m;
    int t = -1;
    for (int p = 0; p < k; p++) {
        t = findmax(matr, block, n, m, p, p, norma, tmp);
        if (t == -1) {
            printf("Не нашлось обратной у findmax\n");
            return 1;
        }
        if (t != p) {
            swap_rows(matr, p, t, n, m);
            swap_rows(solution, p, t, n, m);
        }

        get_block(matr, block, n, l, p, p);

        treug(block, inverse, l, norma, tmp);
        diag(block, inverse, l);

        for (int s = p + 1; s < k + 1; s++) {
            get_block(matr, block, n, m, p, s);

            mult(inverse, block, tmp, m, m, m, m, norma);
            put_block(matr, tmp, n, m, p, s);
        }

        for (int s = 0; s < k + 1; s++) {
            get_block(solution, block, n, m, p, s);
            mult(inverse, block, tmp, m, m, m, m, norma);
            put_block(solution, tmp, n, m, p, s);
        }

        for (int i = p + 1; i < k + 1; i++) { // Умножение m x m
            for (int j = p + 1; j < k + 1; j++) {
                get_block(matr, block, n, m, i, j);
                get_block(matr, block1, n, m, i, p);
                get_block(matr, block2, n, m, p, j);
                mult(block1, block2, tmp, m, m, m, m, norma);
                subtraction(block, tmp, m);
                put_block(matr, block, n, m, i, j);
            }
        }
        for (int i = p + 1; i < k + 1; i++) { // Умножение m x m
            for (int j = 0; j < k + 1; j++) {
                get_block(solution, block, n, m, i, j);

                get_block(matr, block1, n, m, i, p);
                get_block(solution, block2, n, m, p, j);
                mult(block1, block2, tmp, m, m, m, m, norma);
                subtraction(block, tmp, m);
                put_block(solution, block, n, m, i, j);
            }
        }
    }

    get_block(matr, block1, n, m, k, k);

    for (int i = 0; i < bl; i++) {
        for (int j = 0; j < bl; j++) {
            block[i * bl + j] = block1[i * m + j];
        }
    }

    if (treug(block, inverse, bl, norma, tmp) == -1) {
        printf("Метод не применим\n");
        return 1;
    } else {
        diag(block, inverse, bl);
        for (int i = 0; i < k; i++) {
            get_block(solution, block, n, m, k, i);

            mult(inverse, block, tmp, bl, bl, bl, m, norma); // умножение l x l

            put_block(solution, tmp, n, m, k, i);
        }
        get_block(solution, block1, n, m, k, k);
        for (int i = 0; i < bl; i++) {
            for (int j = 0; j < bl; j++) {
                block[i * bl + j] = block1[i * m + j];
            }
        }

        mult(inverse, block, tmp, bl, bl, bl, bl, norma); // умножение l x l

        if (bl == 1) {
            for (int i = 0; i < bl; i++) {
                for (int j = 0; j < bl; j++) {
                    block[i * bl + j] = tmp[i * bl + j];
                }
            }
        } else {
            for (int i = 0; i < bl; i++) {
                for (int j = 0; j < bl; j++) {
                    block[(i)*m + (j)] = tmp[i * bl + j];
                }
            }
        }

        put_block(solution, block, n, m, k, k);
    }

    for (int j = 0; j < k; j++) {
        get_block(solution, block, n, m, k - 1, j);
        get_block(matr, block1, n, m, k - 1, k);
        get_block(solution, block2, n, m, k, j);
        mult(block1, block2, tmp, m, m, m, m, norma);
        subtraction(block, tmp, m);
        put_block(solution, block, n, m, k - 1, j);
    }

    get_block(solution, block, n, m, k - 1, k);
    get_block(matr, block1, n, m, k - 1, k);
    get_block(solution, block2, n, m, k, k);
    mult(block1, block2, tmp, m, m, m, m, norma); // Умножение m x l
    subtraction(block, tmp, m);
    put_block(solution, block, n, m, k - 1, k);

    for (int i = k - 2; i >= 0; i--) {
        for (int j = 0; j < k; j++) {
            get_block(solution, block2, n, m, i, j);

            for (int t = i + 1; t < k; t++) {
                get_block(matr, block, n, m, i, t);
                get_block(solution, block1, n, m, t, j);
                mult(block, block1, tmp, m, m, m, m, norma);
                subtraction(block2, tmp, m);
            }
            get_block(matr, block, n, m, i, k);
            get_block(solution, block1, n, m, k, j);
            mult(block, block1, tmp, m, m, m, m, norma); // Умножение m x l
            subtraction(block2, tmp, m);
            put_block(solution, block2, n, m, i, j);
        }
    }

    if (bl != 0) {
        for (int i = k - 2; i > 0; i--) {
            get_block(solution, block2, n, m, i, k);

            for (int t = i + 1; t < k; t++) {
                get_block(matr, block, n, m, i, t);
                get_block(solution, block1, n, m, t, k);
                mult(block, block1, tmp, m, m, m, m, norma);
                subtraction(block2, tmp, m);
            }

            get_block(matr, block, n, m, i, k);
            get_block(solution, block1, n, m, k, k);
            mult(block, block1, tmp, m, m, m, m, norma);
            subtraction(block2, tmp, m);
            put_block(solution, block2, n, m, i, k);
        }
    }
    return 0;
}
