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
/**
 * @brief Извлекает блок из матрицы.
 * 
 * Эта функция извлекает блок размером `m x m` из матрицы `matr` размером `n x n`,
 * начиная с позиции `(i, j)`. Результат сохраняется в массиве `block`.
 * 
 * @param matr Указатель на матрицу, из которой извлекается блок.
 * @param block Указатель на массив, в который будет сохранен блок.
 * @param n Размер матрицы `matr` (n x n).
 * @param m Размер блока (m x m).
 * @param i Индекс строки начала блока в матрице `matr`.
 * @param j Индекс столбца начала блока в матрице `matr`.
 * 
 * @note Если `m` равно 0, то размер блока устанавливается в 1.
 * @warning Функция предполагает, что матрица `matr` и массив `block` имеют достаточный размер.
 * 
 * @author Ваше Имя
 * @date 2023-10-05
 */
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

void subtraction(double* a, double* b, int row, int col) {
    for (int i=0; i < row; i++) {
        for (int j=0; j < col; j++) {
            a[i * col + j] -= b[i * col + j];
        }
    }  
}
void initialize_matrix(double* solution, int n) {
    std::memset(solution, 0, n * n * sizeof(double));

    for (int i = 0; i < n; ++i) {
        solution[i * n + i] = 1.0;
    }
}

int solve(int n, int m, double* matr, double* block, double* solution, double* inverse, double* tmp, double* block1, double* block2, double norma) {
    int i, j;
    int k = n/m;
    int l = n%m;
    int bl = (l==0?k:k+1);
    // printf("k: %d\nl: %d\nbl: %d\n ", k, l,bl); // h
    for (int p=0; p<bl; p++) {
        if (p != k) 
        {
            int t = findmax(matr, block, n, m, p, p, norma, tmp);
            if (t == -1) {
                printf("Не нашлось обратной у findmax\n");
                return 1;
            }
            if (t != p) {
                swap_rows(matr, p, t, n, m);
                swap_rows(solution, p, t, n, m);
            }
            get_block(matr, block, n, m, p, p);
            treug(block, inverse, m, norma, tmp);
            diag(block, inverse, m);// В inverse обратная на которую надо домножить строку
        
            for (int s = p + 1; s < bl; s++) 
            {
                get_block(matr, block, n, m, p, s);
                int rowb = (p != k ? m : l);
                int colb = (s != k ? m : l);
                mult(inverse, block, tmp, m, m, rowb, colb, norma);
                put_block(matr, tmp, n, m, p, s);
                for (int ss=p+1; ss<bl; ss++) 
                {
                    get_block(matr, block, n , m, ss, s);
                    subtraction(block, tmp, ss!=k?m:l, s!=k?m:l);
                    put_block(matr, block, n ,m , ss, s);
                }
            }
            for (int s = 0; s < bl; s++) 
            {
                get_block(solution, block, n, m, p, s);
                int rowb = (p != k ? m : l);
                int colb = (s != k ? m : l);
                mult(inverse, block, tmp,m, m, rowb, colb, norma);
                put_block(solution, tmp, n, m, p, s);
                for(int ss=p+1; ss<bl; ss++) 
                {
                    get_block(solution, block, n , m, ss, s);
                    subtraction(block, tmp, ss!=k?m:l, s!=k?m:l);
                    put_block(solution, block, n ,m , ss, s);
                }
            }
        } else 
        {
            PrintDouble(solution, n , n );
            get_block(matr, block, n, m, p, p);
            treug(block, inverse, m, norma, tmp); // Размер взятого блока
            diag(block, inverse, m);
            put_block(matr, inverse, n, m, p, p);
            for (i = 0; i < bl; i++)
            {
                get_block(solution, block, n, m, p, i);
                mult(inverse, block, tmp, l, l,p!=k?m:l, i!=k?m:l, norma);
                PrintDouble(tmp, m, m);
                put_block(solution, tmp, n ,m, p, i);
            }
            continue;
        }
    }
    // Обратный ход
    
    return 0;
}
