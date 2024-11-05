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
void get_block(double *matrix, double *block, int n, int m,int i_block, int j_block)
{
	int k = n / m;
	int l = n % m;
	int h = (i_block < k ? m : l);
	int v = (j_block < k ? m : l);
	int tmp;
	for (int i = 0; i < m; i++)
	{
		tmp = i * m;
		for (int j = 0; j < m; j++)
			if (i < h && j < v)
				block[tmp + j] = matrix[n * (i_block * m + i) + j_block * m + j];
			else
				block[tmp + j] = 0;
	}
}
void put_block(double *matrix, double *block, int n, int m, int i_block, int j_block)
{
	int k = n / m;
	int l = n % m;
	int h = (i_block < k ? m : l);
	int v = (j_block < k ? m : l);
	int tmp;
	for (int i = 0; i < h; i++)
	{
		tmp = i * m;
		for (int j = 0; j < v; j++)
			matrix[n * (i_block * m + i) + j_block * m + j] = block[tmp + j];
	}
}

double norma(double *matrix, int n)
{
	double norm = 0, tmp_sum = 0;
	for (int i = 0; i < n; i++)
	{
		tmp_sum = 0;
		for (int j = 0; j < n; j++)
			tmp_sum += fabs(matrix[j * n + i]);
		if (tmp_sum > norm)
			norm = tmp_sum;
	}
	return norm;
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

bool inverse_matrix(double*a, double* b, int n, double norma, double* c)
{
    int t = treug(a, b, n, norma, c);
    if (t == -1) {
        return false;
    }
    diag(a, b, n);
    return true;
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

void subtraction(double *a, double *b, int i_block, int j_block, int k, int m, int l)
{
	int h = (i_block < k ? m : l);
	int v = (j_block < k ? m : l);
	int tmp;
	for (int i = 0; i < h; i++)
	{
		tmp = i * m;
		for (int j = 0; j < v; j++)
			a[tmp + j] = a[tmp + j] - b[tmp + j];
	}
}
void initialize_matrix(double* solution, int n) {
    std::memset(solution, 0, n * n * sizeof(double));

    for (int i = 0; i < n; ++i) {
        solution[i * n + i] = 1.0;
    }
}

void copyMatrix(double* in, double* out, int current_size, int req_size) {
    for (int i = 0; i < current_size*current_size; i++)
    {
        out[i] = 0;
    }
    for (int i = 0; i < req_size; ++i) {
        for (int j = 0; j < req_size; ++j) {
            out[i * req_size + j] = in[i * current_size + j];
        }
    }
}
void setZeroBlock(double *matrix, double *block, int i_block, int j_block, int n, int m)
{
    // Инициализация блока
    for (int row = 0; row < m; row++)
    {
        for (int col = 0; col < m; col++)
        {
            // Проверка, является ли элемент диагональным
            if (row == col && i_block == j_block)
            {
                block[row * m + col] = 1.0; // Установка диагонального элемента в 1
            }
            else
            {
                block[row * m + col] = 0.0; // Установка недиагонального элемента в 0
            }
        }
    }

    // Установка блока в матрицу
    put_block(matrix, block, n, m, i_block, j_block);
}


int solve(int n, int m, double* matr, double* block, double* solution, double* inverse, double* tmp, double* block1, double* block2, double matrix_norm, int *block_index) {
    int s, ss, i, j;
    int buff = 0;
    int k = n/m;
    int l = n%m;
    int bl = (l==0?k:k+1);
    int index_0 = 0, index_1 = 0;
    double buffer_norma, min_norma = 0;
    int min_index = 0;

    for (int p=0; p<k; p++) {
        min_norma = 0;
        min_index = -1;
        for(j=p; j < k; j++)
        {
            buffer_norma = 0;
            get_block(matr, block, n, m, j, p);
            if (gauss_classic_row(block, inverse, block_index, m, matrix_norm, m) != -2)
            {
                buffer_norma = norma(inverse, m);
                if (min_norma > buffer_norma || min_index == -1)
                {
                    min_norma = buffer_norma;
                    min_index = j;
                }
            }
        }
        if (min_index == -1) {
            printf("Метод не применим\n");
            return -1;
        }
        if (min_index != p) {
            swap_rows(matr, min_index, p, n, m);
            swap_rows(solution, min_index, p, n, m);
            // printf("После перестановки\n");
        }

        get_block(matr, block, n, m, p, p);
        // PrintDouble(block, m , m);
        gauss_classic_row(block, inverse, block_index, m, matrix_norm, m);
        // printf("Обратная\n");
        // PrintDouble(inverse, m, m);
        get_block(matr, block, n, m, p, p);
        setZeroBlock(matr, block, p, p, n , m);
        for (s = p+1; s < bl; s++) 
        {
            get_block(matr, block, n, m, p, s);

            // PrintDouble(block, m, m);
            // PrintDouble(inverse, m,m);

            mult(inverse, block, tmp, m, m, s==k?l:m, m);

            // PrintDouble(tmp, m, m);

            put_block(matr, tmp, n, m, p, s);
        }
        for (s = 0; s < bl; s++) 
        {
            get_block(solution, block, n, m, p, s);
            mult(inverse, block, tmp, m, m, m, m);
            put_block(solution, tmp, n, m, p, s);
        }
        // PrintDouble(matr, n , n);
        for (s = p+1; s < bl; s++)
        {
            get_block(matr, block, n , m, s, p);
            get_block(matr, block1, n, m, s, p);
            setZeroBlock(matr, block1, s, p, n, m);
            for (ss=0; ss<bl; ss++) // Тут с p+1
            {
                get_block(matr, block1, n , m, p, ss);
                
                // printf("Начало умножения\n");
                // PrintDouble(block1, m , m);
                // PrintDouble(block, m , m);
                
                mult(block, block1, tmp, m, m, m, m);
                
                // PrintDouble(tmp, m , m);
                // printf("Конец умножения, начало вычитания\n");

                get_block(matr, block1, n ,m, s, ss);
                
                // pcord(s, ss);
                // PrintDouble(block1, m ,m);
                // printf(" - \n");
                // PrintDouble(tmp, m, m);
                
                subtraction(block1, tmp, s, ss, k, m, l);
                
                // PrintDouble(block1, m, m);
                
                put_block(matr, block1, n, m, s, ss);
            }
            for (ss=0; ss<bl; ss++) // Тут с p+1
            {
                get_block(solution, block1, n , m, p, ss);
                mult(block, block1, tmp, m , m, m, m);
                get_block(solution, block1, n ,m, s, ss);
                subtraction(block1, tmp, s, ss, k, m, l);
                put_block(solution, block1, n, m, s, ss);
                
            }
        }
    }
    // PrintDouble(matr, n, n);
    // PrintDouble(solution, n, n);
    // printf("До if ###################\n");
    if (l!=0)
    {
        get_block(matr, block, n ,m, k, k);
        // PrintDouble(block, m, m);
        if (gauss_classic_row(block, inverse, block_index, l, matrix_norm, m) == -2)
        {
            printf("Вот по этосму\n");
            return -1;
        }
        setZeroBlock(matr, block1, k, k, n, m);
        for (j = 0; j <= k; j++) {
            get_block(solution, block, n, m, k, j);
            // PrintDouble(block, m, m);
            // PrintDouble(inverse, m, m);
            mult(inverse, block, tmp, l, l, j==k?l:m, m);
            
            put_block(solution, tmp, n ,m, k, j);
        }
    }
    // printf("После if ###################\n");
    // PrintDouble(matr, n, n);
    // PrintDouble(solution, n, n);
    for (i = n - 1; i >= 0; i--)
		for (k = i - 1; k >= 0; k--)
		{
			index_0 = k * n;
			index_1 = i * n;
			for (j = 0; j < n; j++)
			{
				solution[index_0 + j] -= solution[index_1 + j] * matr[index_0 + i];
			}
			matr[index_0 + i] = 0;
		}

	// for (int i = 0; i < n; i++)
	// {
	// 	index_0 = index_block[i] * n;
	// 	index_1 = i * n;
	// 	for (int j = 0; j < n; j++)
	// 		matr[index_0 + j] = solution[index_1 + j];
	// }

	// for (int i = 0; i < n; i++)
	// {
	// 	index_0 = i * n;
	// 	for (int j = 0; j < n; j++)
	// 		solution[index_0 + j] = matr[index_0 + j];
	// }
    printf("\n");
    return 0;
}
